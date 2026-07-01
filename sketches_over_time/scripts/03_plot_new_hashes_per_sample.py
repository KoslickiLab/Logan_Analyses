#!/usr/bin/env python3
"""
Plot 2: New hashes per sample over time.

"New" = a hash not previously observed in any earlier-released sample.
At year granularity, a hash is counted as new for a sample if
  hash_life.birth_year == sample's release_year.

Produces two figures:
  <stem>__new_hashes_per_sample_median.png  — median + IQR overlay
  <stem>__new_hashes_per_sample_mean.png    — mean + IQR overlay

Each figure has two panels:
  Top    — new hash count per sample (scatter subsampled + summary stat)
  Bottom — fraction of a sample's hashes that are new

The scatter y-axis is capped at --cap (default 100 000) to suppress outliers;
a text annotation reports how many samples exceed the cap.

Usage
-----
  python 03_plot_new_hashes_per_sample.py [--db PATH] [--title STR]
                                          [--max-scatter-points N] [--cap N]
"""

import argparse
import duckdb
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

BASE_DIR  = Path(__file__).resolve().parent.parent
DATA_DIR  = BASE_DIR / 'data'
PLOTS_DIR = DATA_DIR / 'plots'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--db',     type=str, default=None)
    p.add_argument('--title',  type=str, default=None)
    p.add_argument('--max-scatter-points', type=int, default=100_000,
                   help='Max points in the scatter panel (default: 100 000)')
    p.add_argument('--cap', type=int, default=100_000,
                   help='Y-axis cap for new-hash-count panel (default: 100 000)')
    return p.parse_args()


def find_db(db_arg):
    if db_arg:
        return Path(db_arg)
    candidates = sorted(DATA_DIR.glob('intermediate_*.duckdb'))
    if not candidates:
        raise FileNotFoundError(f'No intermediate_*.duckdb found in {DATA_DIR}')
    full = [c for c in candidates if 'full' in c.name]
    return full[0] if full else candidates[-1]


def _fmt_cap(v):
    if v >= 1_000_000:
        return f'{v/1_000_000:.0f}M'
    if v >= 1_000:
        return f'{v/1_000:.0f}k'
    return str(v)


def make_figure(yearly, scatter, stat, cap, suptitle, out_path):
    """
    stat : 'median' or 'mean'
    """
    if stat == 'median':
        center_col   = 'median_new_hashes'
        center_label = 'Median'
        frac_col     = 'median_fraction_new'
        frac_label   = 'Median'
        color_line   = 'tomato'
        color_frac   = 'darkorange'
    else:
        center_col   = 'avg_new_hashes'
        center_label = 'Mean'
        frac_col     = None          # mean fraction derived below
        frac_label   = 'Mean'
        color_line   = 'mediumseagreen'
        color_frac   = 'steelblue'

    yrs = yearly['release_year'].astype(int).values
    rng = np.random.default_rng(42)
    jit = rng.uniform(-0.35, 0.35, len(scatter))

    fig, axes = plt.subplots(2, 1, figsize=(13, 11), sharex=False)

    # ── Top panel: new hash count ──────────────────────────────────────────────
    ax = axes[0]

    # Scatter — clip display at cap, mark clipped points
    new_h  = scatter['new_hashes'].values.copy()
    clipped = new_h > cap
    n_above = int(clipped.sum())
    new_h_disp = np.where(clipped, cap, new_h)

    valid = new_h_disp > 0
    ax.scatter(scatter['release_year'].values[valid] + jit[valid],
               new_h_disp[valid],
               alpha=0.04, s=2, color='steelblue', rasterized=True, linewidths=0)

    # Clipped points shown as triangles at the cap line
    if n_above > 0:
        ax.scatter(scatter['release_year'].values[clipped] + jit[clipped],
                   np.full(n_above, cap),
                   marker='^', s=8, alpha=0.25, color='grey',
                   rasterized=True, linewidths=0, label='_nolegend_')

    # Summary line
    ax.plot(yrs, yearly[center_col],
            'o-', color=color_line, lw=2.5, markersize=7, zorder=5,
            label=center_label)
    if stat == 'median':
        ax.fill_between(yrs, yearly['p25_new_hashes'], yearly['p75_new_hashes'],
                        alpha=0.25, color=color_line, label='IQR (25–75%)')

    ax.set_yscale('log')
    ax.set_ylim(bottom=0.9, top=cap * 1.05)
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: (f'{v/1e6:.0f}M' if v >= 1e6
                                           else f'{v/1e3:.0f}k' if v >= 1e3
                                           else f'{v:.0f}')))
    cap_str = _fmt_cap(cap)
    ax.set_ylabel(f'New Hashes Per Sample  (log scale, capped at {cap_str})',
                  fontsize=12)
    ax.set_title('New Distinct Hashes Per Sample Over Time', fontsize=13,
                 fontweight='bold')
    ax.set_xticks(sorted(yearly['release_year'].astype(int).unique()))
    ax.tick_params(axis='x', rotation=45)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    if n_above > 0:
        ax.text(0.99, 0.97, f'{n_above:,} samples > {cap_str}  (shown at cap)',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=8, color='grey', style='italic')

    # ── Bottom panel: fraction new ─────────────────────────────────────────────
    ax2 = axes[1]
    ax2.scatter(scatter['release_year'].values + jit,
                scatter['fraction_new'].values.clip(0),
                alpha=0.04, s=2, color='seagreen', rasterized=True, linewidths=0)

    if stat == 'median':
        center_frac = yearly['median_fraction_new']
    else:
        # Derive per-year mean fraction from the aggregate table columns
        # avg_new_hashes / avg_total_hashes approximates mean fraction
        center_frac = yearly['avg_new_hashes'] / yearly['avg_total_hashes']

    ax2.plot(yrs, center_frac,
             'o-', color=color_frac, lw=2.5, markersize=7, zorder=5,
             label=frac_label)
    if stat == 'median':
        ax2.fill_between(yrs, yearly['p25_fraction_new'], yearly['p75_fraction_new'],
                         alpha=0.25, color=color_frac, label='IQR (25–75%)')

    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel('Release Year', fontsize=12)
    ax2.set_ylabel('Fraction of Hashes That Are New', fontsize=12)
    ax2.set_title('Fraction of Novel Hashes Per Sample Over Time', fontsize=13,
                  fontweight='bold')
    ax2.set_xticks(sorted(yearly['release_year'].astype(int).unique()))
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    for x, y in zip(yrs, center_frac):
        ax2.annotate(f'{y:.2f}', (x, y), textcoords='offset points',
                     xytext=(0, 8), ha='center', fontsize=7, color=color_frac)

    fig.suptitle(suptitle, fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def main():
    args  = parse_args()
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    db_path = find_db(args.db)
    print(f'Reading: {db_path}')

    con = duckdb.connect(str(db_path), read_only=True)

    yearly = con.execute("""
        SELECT release_year,
               num_samples,
               avg_total_hashes,
               avg_new_hashes,
               median_total_hashes,
               median_new_hashes,
               p25_new_hashes,
               p75_new_hashes,
               median_fraction_new,
               p25_fraction_new,
               p75_fraction_new
        FROM   yearly_sample_stats
        ORDER  BY release_year
    """).df()

    n_total = con.execute('SELECT COUNT(*) FROM sample_stats').fetchone()[0]
    pct     = min(100.0, 100.0 * args.max_scatter_points / max(n_total, 1))
    scatter = con.execute(f"""
        SELECT release_year, new_hashes, fraction_new
        FROM   sample_stats
        USING  SAMPLE {pct:.6f} PERCENT (bernoulli)
    """).df()

    con.close()

    print(f'Yearly stats ({len(yearly)} years):\n{yearly.to_string(index=False)}')
    print(f'Scatter: {len(scatter):,} / {n_total:,} rows  (cap={args.cap:,})')

    stem     = db_path.stem
    base_ttl = args.title or f'New FracMinHash Hashes Per Sample (k=31)\nsource: {stem}'

    make_figure(
        yearly, scatter, stat='median', cap=args.cap,
        suptitle=base_ttl,
        out_path=PLOTS_DIR / f'{stem}__new_hashes_per_sample_median.png',
    )
    make_figure(
        yearly, scatter, stat='mean', cap=args.cap,
        suptitle=base_ttl,
        out_path=PLOTS_DIR / f'{stem}__new_hashes_per_sample_mean.png',
    )


if __name__ == '__main__':
    main()
