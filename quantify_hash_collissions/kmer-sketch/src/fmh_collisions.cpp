
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "FastxReader.hpp"
#include "Hash.hpp"

static void usage() {
    std::cerr <<
        "Usage: fmh_collisions --input FILE [--kmer N] [--scale S] [--seed SEED]\n"
        "                      --summary OUT_SUMMARY --collisions OUT_COLLISIONS\n"
        "\n"
        "Options:\n"
        "  --input FILE        Input FASTA file\n"
        "  --kmer N            K-mer size (default: 31)\n"
        "  --scale S           FracMinHash scale threshold (default: 0.001)\n"
        "  --seed SEED         MurmurHash seed (default: 42)\n"
        "  --summary FILE      Output file for summary statistics\n"
        "  --collisions FILE   Output file for per-collision k-mer lists\n"
        "  --help, -h          Show this help\n"
        "\n"
        "Always uses canonical k-mers (lex-smallest of k-mer and its reverse complement).\n"
        "Ambiguous bases (non-ACGT) are skipped.\n";
}

static std::string get_arg(const std::vector<std::string>& args,
                            const std::string& key,
                            const std::string& def = "") {
    for (size_t i = 0; i < args.size(); ++i)
        if (args[i] == key && i + 1 < args.size()) return args[i + 1];
    return def;
}

static bool has_flag(const std::vector<std::string>& args, const std::string& key) {
    for (const auto& a : args) if (a == key) return true;
    return false;
}

// Returns the canonical form (lex-min of kmer and its reverse complement).
// Returns an empty string if any base is ambiguous (non-ACGT).
// Mirrors the logic in KmerScanner.hpp.
static std::string canonical_kmer(const char* p, size_t k) {
    // Validate forward strand and build reverse complement simultaneously.
    std::string rc(k, 'N');
    for (size_t j = 0; j < k; ++j) {
        char c = p[j];
        switch (c) {
            case 'A': rc[k - 1 - j] = 'T'; break;
            case 'C': rc[k - 1 - j] = 'G'; break;
            case 'G': rc[k - 1 - j] = 'C'; break;
            case 'T': rc[k - 1 - j] = 'A'; break;
            default:  return {};           // ambiguous — skip
        }
    }
    // Lexicographic comparison: choose the smaller of fwd and rc.
    bool take_rc = false;
    for (size_t j = 0; j < k; ++j) {
        char a = p[j], b = rc[j];
        if (a < b) { take_rc = false; break; }
        if (a > b) { take_rc = true;  break; }
    }
    return take_rc ? rc : std::string(p, k);
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }

    std::vector<std::string> args(argv + 1, argv + argc);
    if (has_flag(args, "--help") || has_flag(args, "-h")) { usage(); return 0; }

    const std::string inpath      = get_arg(args, "--input");
    const std::string summary_out = get_arg(args, "--summary");
    const std::string collis_out  = get_arg(args, "--collisions");

    if (inpath.empty() || summary_out.empty() || collis_out.empty()) {
        std::cerr << "Error: --input, --summary, and --collisions are all required.\n\n";
        usage();
        return 1;
    }

    const size_t   k     = (size_t)std::stoull(get_arg(args, "--kmer",  "31"));
    const double   scale = std::stod(get_arg(args, "--scale", "0.001"));
    const uint64_t seed  = (uint64_t)std::stoull(get_arg(args, "--seed", "42"));

    if (scale <= 0.0 || scale > 1.0) {
        std::cerr << "Error: scale must be in (0, 1].\n";
        return 1;
    }

    // Hashes <= threshold are selected into the sketch.
    const uint64_t HASH_MAX   = std::numeric_limits<uint64_t>::max();
    const uint64_t threshold  = (uint64_t)(scale * (double)HASH_MAX);

    std::cerr << "Parameters: k=" << k << ", scale=" << scale
              << ", seed=" << seed << ", threshold=" << threshold << "\n"
              << "Input: " << inpath << "\n";

    // Core data structure:
    //   hash value -> set of distinct canonical k-mer strings that hashed to it
    // Only entries with hash <= threshold are stored (i.e., the FMH sketch).
    std::unordered_map<uint64_t, std::unordered_set<std::string>> sketch;

    uint64_t total_kmer_positions = 0; // all positions scanned (counting duplicates)
    uint64_t selected_positions   = 0; // positions where hash <= threshold

    {
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            // Uppercase in-place for consistent hashing.
            for (char& c : seq)
                if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');

            const size_t n = seq.size();
            for (size_t i = 0; i + k <= n; ++i) {
                ++total_kmer_positions;
                std::string canon = canonical_kmer(&seq[i], k);
                if (canon.empty()) continue; // skip ambiguous k-mers

                const uint64_t h = hashutil::murmurhash64(canon.data(), k, seed);
                if (h <= threshold) {
                    ++selected_positions;
                    sketch[h].insert(std::move(canon));
                }
            }
        }
    }

    // Collect collision entries (hashes with >1 distinct canonical k-mer).
    std::vector<std::pair<uint64_t, std::vector<std::string>>> collisions;
    uint64_t n_collision_hashes = 0;
    uint64_t n_kmers_in_collisions = 0;
    uint64_t total_distinct_kmers_selected = 0;

    for (auto& [h, kmers] : sketch) {
        total_distinct_kmers_selected += kmers.size();
        if (kmers.size() > 1) {
            ++n_collision_hashes;
            n_kmers_in_collisions += kmers.size();
            std::vector<std::string> sorted(kmers.begin(), kmers.end());
            std::sort(sorted.begin(), sorted.end());
            collisions.emplace_back(h, std::move(sorted));
        }
    }
    // Sort collision records by hash value for deterministic output.
    std::sort(collisions.begin(), collisions.end());

    const uint64_t sketch_size = sketch.size(); // distinct hash values below threshold

    // --- Write summary file ---
    std::ofstream sum_f(summary_out);
    if (!sum_f) { std::cerr << "Cannot open summary output: " << summary_out << "\n"; return 2; }

    sum_f << "# FracMinHash collision analysis\n"
          << "# input=" << inpath << "\n"
          << "# kmer_size=" << k << "\n"
          << "# scale=" << scale << "\n"
          << "# seed=" << seed << "\n"
          << "# threshold=" << threshold << "\n"
          << "# hash_function=MurmurHash3_x64_128_low64\n"
          << "# total_kmer_positions_scanned=" << total_kmer_positions << "\n"
          << "# selected_positions (hash<=threshold, with duplicates)=" << selected_positions << "\n"
          << "# total_distinct_canonical_kmers_selected=" << total_distinct_kmers_selected << "\n"
          << "sketch_size\t" << sketch_size << "\n"
          << "n_collision_hashes\t" << n_collision_hashes << "\n"
          << "n_kmers_in_collisions\t" << n_kmers_in_collisions << "\n";

    // --- Write collision detail file ---
    std::ofstream col_f(collis_out);
    if (!col_f) { std::cerr << "Cannot open collisions output: " << collis_out << "\n"; return 2; }

    col_f << "# FracMinHash collision details\n"
          << "# Each line: a hash value with 2+ distinct canonical k-mers\n"
          << "# Columns: hash_value TAB kmer1 TAB kmer2 [TAB kmer3 ...]\n";

    for (auto& [h, kmers] : collisions) {
        col_f << h;
        for (const auto& km : kmers) col_f << '\t' << km;
        col_f << '\n';
    }

    std::cerr << "Done.\n"
              << "  Sketch size (distinct hashes):     " << sketch_size << "\n"
              << "  Collision hashes (>1 k-mer):        " << n_collision_hashes << "\n"
              << "  K-mers involved in collisions:      " << n_kmers_in_collisions << "\n";

    return 0;
}
