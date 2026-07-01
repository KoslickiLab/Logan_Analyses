#!/usr/bin/env python3
"""
Plot 1: Cumulative distinct FracMinHash hashes over time.

Reads 'cumulative_hashes_by_year' from an intermediate DuckDB database and
produces a two-panel figure:
  Left  — cumulative distinct hashes (log y-axis)
  Right — new distinct hashes introduced each year (log y-axis, bar chart)

Usage
-----
  python 02_plot_cumulative_hashes.py [--db PATH] [--output PNG] [--title STR]
"""

import argparse
import duckdb
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
    p.add_argument('--db',     type=str, default=None,
                   help='Path to intermediate DuckDB (auto-detected if omitted)')
    p.add_argument('--output', type=str, default=None,
                   help='Output PNG path (auto-named if omitted)')
    p.add_argument('--title',  type=str, default=None,
                   help='Override the plot suptitle')
    return p.parse_args()


def find_db(db_arg):
    if db_arg:
        return Path(db_arg)
    candidates = sorted(DATA_DIR.glob('intermediate_*.duckdb'))
    if not candidates:
        raise FileNotFoundError(f'No intermediate_*.duckdb found in {DATA_DIR}')
    # prefer 'full' variant
    full = [c for c in candidates if 'full' in c.name]
    return full[0] if full else candidates[-1]


def main():
    args  = parse_args()
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    db_path = find_db(args.db)
    print(f'Reading: {db_path}')

    con = duckdb.connect(str(db_path), read_only=True)
    df  = con.execute("""
        SELECT year, new_hashes_this_year, cumulative_distinct_hashes
        FROM   cumulative_hashes_by_year
        ORDER  BY year
    """).df()
    con.close()

    print(df.to_string(index=False))

    out_path = (Path(args.output) if args.output
                else PLOTS_DIR / f'{db_path.stem}__cumulative_hashes.png')
    suptitle = (args.title or
                f'Cumulative Distinct FracMinHashes Over Time')

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    years = df['year'].astype(int)
    year_ticks = sorted(years.unique())

    # ── Left: cumulative distinct hashes ──────────────────────────────────────
    ax = axes[0]
    ax.plot(years, df['cumulative_distinct_hashes'],
            'o-', color='steelblue', linewidth=2.5, markersize=7, zorder=3)
    for x, y in zip(years, df['cumulative_distinct_hashes']):
        ax.annotate(f'{y/1e9:.1f}B' if y >= 1e9 else f'{y/1e6:.0f}M',
                    (x, y), textcoords='offset points', xytext=(0, 8),
                    ha='center', fontsize=8, color='steelblue')
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: (f'{v/1e9:.1f}B' if v >= 1e9
                                           else f'{v/1e6:.0f}M' if v >= 1e6
                                           else f'{v:.0f}')))
    ax.set_xlabel('Release Year', fontsize=12)
    ax.set_ylabel('Cumulative Distinct Hashes  (log scale)', fontsize=12)
    ax.set_title('Cumulative Distinct Hashes', fontsize=13, fontweight='bold')
    ax.set_xticks(year_ticks)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.35, which='both')

    # ── Right: new hashes per year ─────────────────────────────────────────────
    ax2 = axes[1]
    bars = ax2.bar(years, df['new_hashes_this_year'],
                   color='steelblue', alpha=0.75, edgecolor='steelblue')
    ax2.set_yscale('log')
    ax2.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: (f'{v/1e9:.1f}B' if v >= 1e9
                                           else f'{v/1e6:.0f}M' if v >= 1e6
                                           else f'{v:.0f}')))
    ax2.set_xlabel('Release Year', fontsize=12)
    ax2.set_ylabel('New Distinct Hashes This Year  (log scale)', fontsize=12)
    ax2.set_title('New Distinct Hashes Per Year', fontsize=13, fontweight='bold')
    ax2.set_xticks(year_ticks)
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(True, alpha=0.35, axis='y', which='both')

    fig.suptitle(suptitle, fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    main()
