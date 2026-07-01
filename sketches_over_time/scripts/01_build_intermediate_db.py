#!/usr/bin/env python3
"""
Build an intermediate DuckDB database with hash birth/death statistics.

Two heavy computation stages scan the 193B-row signature table:
  Stage 1 — hash_life     : (min_hash, birth_year, death_year)
                             GROUP BY over all 193B rows.
  Stage 2 — sample_stats  : (acc, release_year, total_hashes, new_hashes,
                              fraction_new)
                             Full-scan join of signature_mins with hash_life.

Then lightweight derived tables are built from those two:
  Stage 3 — cumulative_hashes_by_year  : cumulative distinct hashes per year
           — yearly_sample_stats        : per-year aggregate stats over samples

The key intermediate table is 'sample_year', which encodes the
metadata filter (and test-mode limit) once; all subsequent queries join
against it rather than re-querying the metadata database.

Usage examples
--------------
  # Smoke test — 500 samples, fast:
  python 01_build_intermediate_db.py --test --num-test-samples 500

  # Full run (expect many hours):
  python 01_build_intermediate_db.py

  # Filter to ILLUMINA platform only:
  python 01_build_intermediate_db.py --metadata-filter "platform='ILLUMINA'"

  # Resume after interruption — skips stages whose tables already exist:
  python 01_build_intermediate_db.py --resume

  # Run only stage 1 (useful to checkpoint before stage 2):
  python 01_build_intermediate_db.py --stage 1

  # Run only stage 2 (after stage 1 completed):
  python 01_build_intermediate_db.py --stage 2 --resume

Metadata filter notes
---------------------
The --metadata-filter string is injected verbatim as an extra AND clause
into the WHERE condition of the metadata_geo_joined table.  Reference
column names without any table prefix.  Examples:
    "platform='ILLUMINA'"
    "libraryselection='RANDOM'"
    "platform='ILLUMINA' AND librarysource='METAGENOMIC'"
"""

import argparse
import duckdb
import sys
import time
from pathlib import Path

# ── Fixed paths ────────────────────────────────────────────────────────────────
SRC_DB  = '/scratch/shared_data/Logan_yacht_data/processed_data/database_all.db'
META_DB = '/scratch/shared_data/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined_5M.duckdb'

BASE_DIR = Path(__file__).resolve().parent.parent   # sketches_over_time/
TMP_DIR  = BASE_DIR / 'tmp'
DATA_DIR = BASE_DIR / 'data'
CSV_DIR  = DATA_DIR / 'csv'
LOG_DIR  = BASE_DIR / 'logs'

KSIZE = 31

# ── Logging ────────────────────────────────────────────────────────────────────
def log(msg, level='INFO'):
    ts = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f'[{ts}] [{level}] {msg}', flush=True)


class _Timer:
    def __init__(self, label):
        self.label = label
    def __enter__(self):
        self.t0 = time.time()
        log(f'START  {self.label}')
        return self
    def __exit__(self, *_):
        dt = time.time() - self.t0
        h, rem = divmod(int(dt), 3600)
        m, s   = divmod(rem, 60)
        log(f'DONE   {self.label}   [{h:02d}:{m:02d}:{s:02d}]')

def timed(label):
    return _Timer(label)


# ── Argument parsing ───────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--test', action='store_true',
                   help='Run on a small subset of samples (smoke test)')
    p.add_argument('--num-test-samples', type=int, default=500,
                   help='Number of samples in test mode (default: 500)')
    p.add_argument('--metadata-filter', type=str, default='',
                   help='Extra SQL WHERE clause on metadata columns (no table prefix). '
                        "Example: \"platform='ILLUMINA'\"")
    p.add_argument('--output-db', type=str, default=None,
                   help='Output DuckDB path (auto-named from flags if not set)')
    p.add_argument('--threads', type=int, default=768,
                   help='DuckDB thread count (default: 768)')
    p.add_argument('--memory-limit', type=str, default='2500GB',
                   help='DuckDB memory limit (default: 2500GB)')
    p.add_argument('--stage', type=int, choices=[1, 2, 3], default=None,
                   help='Run only this stage (1=hash_life, 2=sample_stats, 3=aggregates). '
                        'Default: run all three stages.')
    p.add_argument('--resume', action='store_true',
                   help='Skip stages whose output tables already exist')
    p.add_argument('--no-csv', action='store_true',
                   help='Skip CSV export of derived/aggregate tables')
    p.add_argument('--export-hash-life-csv', action='store_true',
                   help='Also export hash_life to CSV (can be very large for full runs)')
    p.add_argument('--export-sample-stats-csv', action='store_true',
                   help='Also export sample_stats to CSV (~5M rows)')
    return p.parse_args()


# ── DB helpers ─────────────────────────────────────────────────────────────────
def tables(con):
    return {r[0] for r in con.execute('SHOW TABLES').fetchall()}

def nrows(con, t):
    return con.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0]

def drop_if_exists(con, t):
    if t in tables(con):
        con.execute(f'DROP TABLE {t}')
        log(f'Dropped existing table: {t}')

def skip_or_drop(con, t, resume):
    """If resume and table exists, return True (skip).  Otherwise drop and return False."""
    if t in tables(con):
        if resume:
            log(f'{t}: already exists ({nrows(con, t):,} rows), skipping  (--resume)')
            return True
        drop_if_exists(con, t)
    return False


# ── Output-DB naming ───────────────────────────────────────────────────────────
def make_output_path(args):
    if args.output_db:
        return Path(args.output_db)
    parts = []
    if args.test:
        parts.append(f'test{args.num_test_samples}')
    if args.metadata_filter:
        slug = (args.metadata_filter
                .replace("'", '').replace('"', '')
                .replace('=', '_').replace(' ', '_')
                .replace('(', '').replace(')', '')
                .strip('_')[:40])
        parts.append(slug)
    suffix = '_'.join(parts) if parts else 'full'
    return DATA_DIR / f'intermediate_{suffix}.duckdb'


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    for d in [TMP_DIR, DATA_DIR, CSV_DIR, LOG_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    out_db = make_output_path(args)
    log(f'Output DB         : {out_db}')
    log(f'Source signatures : {SRC_DB}')
    log(f'Metadata DB       : {META_DB}')
    log(f'Test mode         : {args.test}  (n={args.num_test_samples})'
        if args.test else 'Test mode         : off (full dataset)')
    log(f'Metadata filter   : {args.metadata_filter or "(none)"}')
    log(f'Threads           : {args.threads}')
    log(f'Memory limit      : {args.memory_limit}')
    log(f'Stage(s)          : {args.stage or "all"}')
    log(f'Resume mode       : {args.resume}')

    con = duckdb.connect(str(out_db))

    # Configure parallelism and memory
    con.execute(f'SET threads = {args.threads}')
    con.execute(f"SET memory_limit = '{args.memory_limit}'")
    con.execute(f"SET temp_directory = '{TMP_DIR}'")
    con.execute('SET enable_progress_bar = true')
    con.execute('SET enable_progress_bar_print = true')

    # Attach source databases read-only
    with timed('Attach source databases'):
        con.execute(f"ATTACH '{SRC_DB}'  AS src  (READ_ONLY)")
        con.execute(f"ATTACH '{META_DB}' AS meta (READ_ONLY)")

    # ── sample_year ────────────────────────────────────────────────────────────
    # Pre-compute the (acc → release_year) mapping with filters applied once.
    # All subsequent joins use this table, keeping stage queries filter-agnostic.
    if not skip_or_drop(con, 'sample_year', args.resume):
        meta_filter = f"AND ({args.metadata_filter})" if args.metadata_filter else ''
        test_limit  = f"LIMIT {args.num_test_samples}"   if args.test           else ''
        with timed('sample_year'):
            con.execute(f"""
                CREATE TABLE sample_year AS
                SELECT acc, YEAR(releasedate) AS release_year
                FROM   meta.metadata_geo_joined
                WHERE  releasedate IS NOT NULL
                {meta_filter}
                {test_limit}
            """)
        log(f'sample_year rows  : {nrows(con, "sample_year"):,}')

    run1 = args.stage in (None, 1)
    run2 = args.stage in (None, 2)
    run3 = args.stage in (None, 3)

    # ── Stage 1: hash_life ─────────────────────────────────────────────────────
    # Full scan of 193B rows, GROUP BY min_hash.
    # Expected runtime: hours on the full dataset.
    if run1:
        if not skip_or_drop(con, 'hash_life', args.resume):
            with timed('Stage 1 — hash_life  (GROUP BY min_hash over all rows)'):
                con.execute(f"""
                    CREATE TABLE hash_life AS
                    SELECT
                        sm.min_hash,
                        MIN(sy.release_year) AS birth_year,
                        MAX(sy.release_year) AS death_year
                    FROM src.sigs_dna.signature_mins sm
                    INNER JOIN sample_year sy ON sm.sample_id = sy.acc
                    WHERE sm.ksize = {KSIZE}
                    GROUP BY sm.min_hash
                """)
            n = nrows(con, 'hash_life')
            log(f'hash_life rows    : {n:,}  (distinct hashes)')

    # ── Stage 2: sample_stats ──────────────────────────────────────────────────
    # Second full scan, now with a join against hash_life to flag new hashes.
    if run2:
        if 'hash_life' not in tables(con):
            log('hash_life missing — run stage 1 first', 'ERROR')
            sys.exit(1)
        if not skip_or_drop(con, 'sample_stats', args.resume):
            with timed('Stage 2 — sample_stats  (join 193B rows with hash_life)'):
                con.execute(f"""
                    CREATE TABLE sample_stats AS
                    SELECT
                        sm.sample_id                                                         AS acc,
                        sy.release_year,
                        COUNT(*)                                                             AS total_hashes,
                        SUM(CASE WHEN hl.birth_year = sy.release_year THEN 1 ELSE 0 END)    AS new_hashes,
                        AVG(CASE WHEN hl.birth_year = sy.release_year THEN 1.0 ELSE 0.0 END) AS fraction_new
                    FROM src.sigs_dna.signature_mins sm
                    INNER JOIN sample_year sy ON sm.sample_id  = sy.acc
                    INNER JOIN hash_life   hl ON sm.min_hash   = hl.min_hash
                    WHERE sm.ksize = {KSIZE}
                    GROUP BY sm.sample_id, sy.release_year
                """)
            n = nrows(con, 'sample_stats')
            log(f'sample_stats rows : {n:,}  (samples with release-date)')

    # ── Stage 3: lightweight derived tables ────────────────────────────────────
    if run3:
        # 3a — cumulative distinct hashes by year
        if not skip_or_drop(con, 'cumulative_hashes_by_year', args.resume):
            if 'hash_life' not in tables(con):
                log('hash_life missing — run stage 1 first', 'ERROR')
                sys.exit(1)
            with timed('Stage 3a — cumulative_hashes_by_year'):
                con.execute("""
                    CREATE TABLE cumulative_hashes_by_year AS
                    SELECT
                        birth_year                              AS year,
                        COUNT(*)                                AS new_hashes_this_year,
                        SUM(COUNT(*)) OVER (
                            ORDER BY birth_year
                            ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
                        )                                       AS cumulative_distinct_hashes
                    FROM hash_life
                    WHERE birth_year IS NOT NULL
                    GROUP BY birth_year
                    ORDER BY birth_year
                """)
            log(f'cumulative_hashes_by_year rows: {nrows(con, "cumulative_hashes_by_year"):,}')

        # 3b — per-year aggregate stats over samples
        if not skip_or_drop(con, 'yearly_sample_stats', args.resume):
            if 'sample_stats' not in tables(con):
                log('sample_stats missing — run stage 2 first', 'ERROR')
                sys.exit(1)
            with timed('Stage 3b — yearly_sample_stats'):
                con.execute("""
                    CREATE TABLE yearly_sample_stats AS
                    SELECT
                        release_year,
                        COUNT(*)                                                        AS num_samples,
                        SUM(total_hashes)                                               AS total_hashes_sum,
                        SUM(new_hashes)                                                 AS total_new_hashes,
                        AVG(total_hashes)                                               AS avg_total_hashes,
                        AVG(new_hashes)                                                 AS avg_new_hashes,
                        PERCENTILE_CONT(0.5)  WITHIN GROUP (ORDER BY total_hashes)     AS median_total_hashes,
                        PERCENTILE_CONT(0.5)  WITHIN GROUP (ORDER BY new_hashes)       AS median_new_hashes,
                        PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY new_hashes)       AS p25_new_hashes,
                        PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY new_hashes)       AS p75_new_hashes,
                        PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY fraction_new)     AS p25_fraction_new,
                        PERCENTILE_CONT(0.5)  WITHIN GROUP (ORDER BY fraction_new)     AS median_fraction_new,
                        PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY fraction_new)     AS p75_fraction_new
                    FROM sample_stats
                    GROUP BY release_year
                    ORDER BY release_year
                """)
            log(f'yearly_sample_stats rows: {nrows(con, "yearly_sample_stats"):,}')

        # 3c — birth/death aggregate for survival plots (compact: no per-hash rows)
        if not skip_or_drop(con, 'hash_life_counts', args.resume):
            if 'hash_life' not in tables(con):
                log('hash_life missing — run stage 1 first', 'ERROR')
                sys.exit(1)
            with timed('Stage 3c — hash_life_counts  (birth_year × death_year × duration)'):
                con.execute("""
                    CREATE TABLE hash_life_counts AS
                    SELECT
                        birth_year,
                        death_year,
                        (death_year - birth_year) AS duration,
                        COUNT(*)                  AS n_hashes
                    FROM hash_life
                    WHERE birth_year IS NOT NULL
                      AND death_year IS NOT NULL
                    GROUP BY birth_year, death_year
                    ORDER BY birth_year, death_year
                """)
            log(f'hash_life_counts rows: {nrows(con, "hash_life_counts"):,}')

    # ── CSV exports ────────────────────────────────────────────────────────────
    if not args.no_csv:
        stem = out_db.stem

        # Always export small aggregate tables
        for tbl in ['cumulative_hashes_by_year', 'yearly_sample_stats',
                    'hash_life_counts', 'sample_year']:
            if tbl in tables(con):
                csv = CSV_DIR / f'{stem}__{tbl}.csv'
                with timed(f'Export {tbl}'):
                    con.execute(f"COPY {tbl} TO '{csv}' (HEADER, DELIMITER ',')")
                log(f'  → {csv.name}  ({nrows(con, tbl):,} rows)')

        # Optional large exports
        if args.export_hash_life_csv and 'hash_life' in tables(con):
            csv = CSV_DIR / f'{stem}__hash_life.csv'
            with timed('Export hash_life (large)'):
                con.execute(f"COPY hash_life TO '{csv}' (HEADER, DELIMITER ',')")
            log(f'  → {csv.name}  ({nrows(con, "hash_life"):,} rows)')

        if args.export_sample_stats_csv and 'sample_stats' in tables(con):
            csv = CSV_DIR / f'{stem}__sample_stats.csv'
            with timed('Export sample_stats'):
                con.execute(f"COPY sample_stats TO '{csv}' (HEADER, DELIMITER ',')")
            log(f'  → {csv.name}  ({nrows(con, "sample_stats"):,} rows)')

    # ── Summary ────────────────────────────────────────────────────────────────
    log('\n=== Table summary ===')
    for row in con.execute('SHOW TABLES').fetchall():
        t = row[0]
        log(f'  {t:<35s}: {nrows(con, t):>15,} rows')

    con.close()
    log(f'Done.  Output DB: {out_db}')


if __name__ == '__main__':
    main()
