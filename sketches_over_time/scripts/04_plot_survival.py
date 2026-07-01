#!/usr/bin/env python3
"""
Plot 3: Kaplan-Meier survival curves and birth/death heatmap.

A hash is "born" in its first observed year (birth_year) and "dies" in its
last observed year (death_year).  Hashes whose death_year equals the dataset
maximum year are right-censored (we cannot tell whether they will reappear).

Survival interpretation:
  "Given a hash was just born, what is the probability it will still be
   observed T years later?"

Produces three output files:
  <db_stem>__survival_km_overall.png  — overall KM curve with censor marks
  <db_stem>__survival_km_by_birth.png — KM curves stratified by birth year
  <db_stem>__birth_death_heatmap.png  — 2-D heatmap of birth vs. death year

KM is implemented from scratch using the 'hash_life_counts' aggregate table
(counts of hashes by (birth_year, death_year)) to avoid loading billions of
individual hash rows.

Usage
-----
  python 04_plot_survival.py [--db PATH] [--title STR]
"""

import argparse
import duckdb
import numpy as np
import pandas as pd
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
    p.add_argument('--db',    type=str, default=None)
    p.add_argument('--title', type=str, default=None)
    return p.parse_args()


def find_db(db_arg):
    if db_arg:
        return Path(db_arg)
    candidates = sorted(DATA_DIR.glob('intermediate_*.duckdb'))
    if not candidates:
        raise FileNotFoundError(f'No intermediate_*.duckdb found in {DATA_DIR}')
    full = [c for c in candidates if 'full' in c.name]
    return full[0] if full else candidates[-1]


# ── Kaplan-Meier from aggregated (duration, event, count) data ─────────────────
def km_from_counts(df_counts):
    """
    Compute KM survival from aggregated (duration, event, n_hashes) data.

    Parameters
    ----------
    df_counts : DataFrame with columns  duration (int), event (0/1), n_hashes (int)
        event=1 → hash stopped being observed before max_year  (death)
        event=0 → hash was still observed in max_year          (censored)

    Returns
    -------
    DataFrame: time, n_at_risk, n_events, n_censored, survival
    """
    agg = (df_counts
           .groupby(['duration', 'event'])['n_hashes']
           .sum()
           .reset_index())
    pivot = (agg
             .pivot_table(index='duration', columns='event',
                          values='n_hashes', aggfunc='sum', fill_value=0)
             .reset_index()
             .sort_values('duration'))
    pivot.columns.name = None
    if 0 not in pivot.columns:
        pivot[0] = 0
    if 1 not in pivot.columns:
        pivot[1] = 0

    at_risk = int(pivot[[0, 1]].values.sum())
    S = 1.0
    rows = []
    for _, row in pivot.iterrows():
        n_ev  = int(row[1])
        n_cen = int(row[0])
        n_exit = n_ev + n_cen
        if n_ev > 0 and at_risk > 0:
            S *= (1.0 - n_ev / at_risk)
        rows.append({'time': int(row['duration']),
                     'n_at_risk': at_risk,
                     'n_events': n_ev,
                     'n_censored': n_cen,
                     'survival': S})
        at_risk -= n_exit
    return pd.DataFrame(rows)


def step_coords(km_df, prepend_t0=True):
    """Return (times, survivals) arrays suitable for a step plot."""
    t = km_df['time'].tolist()
    s = km_df['survival'].tolist()
    if prepend_t0:
        t = [0] + t
        s = [1.0] + s
    return np.array(t), np.array(s)


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    db_path = find_db(args.db)
    print(f'Reading: {db_path}')

    con = duckdb.connect(str(db_path), read_only=True)

    max_death = con.execute(
        'SELECT MAX(death_year) FROM hash_life_counts WHERE death_year IS NOT NULL'
    ).fetchone()[0]
    min_birth = con.execute(
        'SELECT MIN(birth_year) FROM hash_life_counts WHERE birth_year IS NOT NULL'
    ).fetchone()[0]
    max_birth = con.execute(
        'SELECT MAX(birth_year) FROM hash_life_counts WHERE birth_year IS NOT NULL'
    ).fetchone()[0]
    print(f'Birth range: {min_birth}–{max_birth}   Max death year: {max_death}')

    # Load aggregate counts
    counts = con.execute("""
        SELECT birth_year, death_year, duration, n_hashes
        FROM   hash_life_counts
        ORDER  BY birth_year, death_year
    """).df()
    con.close()

    print(f'hash_life_counts rows: {len(counts):,}')
    print(counts.to_string(index=False))

    # Derive event column: event=1 if died before max_death, 0 if censored
    counts['event'] = (counts['death_year'] < max_death).astype(int)

    stem    = db_path.stem
    title_s = args.title or f'Hash Survival Analysis (k=31) — {stem}'

    # ── Figure 1: overall KM ──────────────────────────────────────────────────
    km_all = km_from_counts(counts[['duration', 'event', 'n_hashes']].copy())
    print('\nOverall KM table:')
    print(km_all.to_string(index=False))

    fig, ax = plt.subplots(figsize=(10, 7))
    t, s = step_coords(km_all)
    ax.step(t, s, where='post', color='steelblue', lw=2.5,
            label=f'All hashes  (N={counts["n_hashes"].sum():,})')

    # Censor tick marks at each time with censored events
    cens = km_all[km_all['n_censored'] > 0]
    if not cens.empty:
        # Find survival at each censor time
        for _, crow in cens.iterrows():
            mask = km_all['time'] <= crow['time']
            s_at_t = km_all.loc[mask, 'survival'].iloc[-1] if mask.any() else 1.0
            ax.scatter([crow['time']], [s_at_t],
                       marker='+', s=80, color='steelblue', zorder=6,
                       linewidths=1.5, label='_nolegend_')

    ax.set_xlabel('Years Since First Observation  (birth_year = 0)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.set_title('Kaplan-Meier Survival Curve\n'
                 'P(hash still observed T years after first appearance)', fontsize=13,
                 fontweight='bold')
    ax.set_xlim(left=-0.2)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(range(0, int(max_death - min_birth) + 2))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)

    # Annotate final survival value
    s_final = km_all['survival'].iloc[-1]
    ax.annotate(f'S(end) = {s_final:.3f}',
                xy=(t[-1], s[-1]),
                xytext=(t[-1] - 1.5, s[-1] + 0.06),
                fontsize=10, color='steelblue',
                arrowprops=dict(arrowstyle='->', color='steelblue', lw=1.2))

    plt.tight_layout()
    out_km = PLOTS_DIR / f'{stem}__survival_km_overall.png'
    plt.savefig(out_km, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_km}')

    # ── Figure 2: KM stratified by birth year ─────────────────────────────────
    birth_years = sorted(counts['birth_year'].dropna().unique().astype(int))
    n_by        = len(birth_years)
    cmap        = plt.cm.viridis

    fig2, ax2 = plt.subplots(figsize=(12, 8))
    for i, by in enumerate(birth_years):
        grp  = counts[counts['birth_year'] == by][['duration', 'event', 'n_hashes']].copy()
        if grp.empty or grp['n_hashes'].sum() == 0:
            continue
        km_g = km_from_counts(grp)
        tg, sg = step_coords(km_g)
        color  = cmap(i / max(n_by - 1, 1))
        # Max possible duration for this birth cohort
        max_dur = int(max_death - by)
        # Extend the KM line to max_dur so all curves end at the same x
        if len(tg) > 0 and tg[-1] < max_dur:
            tg = np.append(tg, max_dur)
            sg = np.append(sg, sg[-1])
        ax2.step(tg, sg, where='post', color=color, lw=1.8, alpha=0.85,
                 label=str(by))

    ax2.set_xlabel('Years Since First Observation  (birth_year = 0)', fontsize=12)
    ax2.set_ylabel('Survival Probability', fontsize=12)
    ax2.set_title('Kaplan-Meier Curves Stratified by Birth Year', fontsize=13,
                  fontweight='bold')
    ax2.set_xlim(left=-0.2)
    ax2.set_ylim(0, 1.05)
    ax2.set_xticks(range(0, int(max_death - min_birth) + 2))
    ax2.grid(True, alpha=0.3)

    if n_by <= 16:
        ax2.legend(title='Birth Year', fontsize=8, ncol=2,
                   loc='upper right', framealpha=0.8)
    else:
        sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=plt.Normalize(vmin=min(birth_years), vmax=max(birth_years)))
        sm.set_array([])
        plt.colorbar(sm, ax=ax2, label='Birth Year', shrink=0.8)

    fig2.suptitle(title_s, fontsize=13, fontweight='bold')
    plt.tight_layout()
    out_km2 = PLOTS_DIR / f'{stem}__survival_km_by_birth.png'
    plt.savefig(out_km2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_km2}')

    # ── Figure 3: birth × death heatmap ───────────────────────────────────────
    all_years = list(range(int(min_birth), int(max_death) + 1))
    n_yr      = len(all_years)
    yr_idx    = {y: i for i, y in enumerate(all_years)}

    mat = np.zeros((n_yr, n_yr), dtype=np.float64)
    for _, row in counts.iterrows():
        bi = yr_idx.get(int(row['birth_year']), -1)
        di = yr_idx.get(int(row['death_year']), -1)
        if bi >= 0 and di >= 0:
            mat[bi, di] += row['n_hashes']

    # Only show upper-left triangle (birth_year <= death_year)
    mat_log = np.log10(mat + 1)

    fig3, ax3 = plt.subplots(figsize=(10, 9))
    im = ax3.imshow(
        mat_log,
        aspect='auto',
        origin='lower',
        cmap='YlOrRd',
        extent=[all_years[0] - 0.5, all_years[-1] + 0.5,
                all_years[0] - 0.5, all_years[-1] + 0.5],
    )
    cbar = plt.colorbar(im, ax=ax3, shrink=0.85)
    cbar.set_label('Number of hashes', fontsize=11)
    # Fixed round-number ticks; positions are log10(count + 1)
    round_counts = [0, 100, 10_000, 1_000_000, 100_000_000, 1_000_000_000, 10_000_000_000]
    round_labels = ['0', '100', '10K', '1M', '100M', '1B', '10B']
    tick_pos = [np.log10(c + 1) for c in round_counts]
    vmin, vmax = im.norm.vmin, im.norm.vmax
    valid = [(p, l) for p, l in zip(tick_pos, round_labels) if vmin <= p <= vmax + 0.01]
    cbar.set_ticks([p for p, _ in valid])
    cbar.set_ticklabels([l for _, l in valid])

    # Diagonal: birth == death (duration = 0)
    ax3.plot([all_years[0], all_years[-1]],
             [all_years[0], all_years[-1]],
             'b--', alpha=0.5, lw=1.2, label='birth = death')

    ax3.set_xlabel('Death Year  (last observed year)', fontsize=12)
    ax3.set_ylabel('Birth Year  (first observed year)', fontsize=12)
    ax3.set_title('Hash Birth Year vs Death Year', fontsize=13,
                  fontweight='bold')
    ax3.set_xticks(all_years)
    ax3.set_yticks(all_years)
    ax3.tick_params(axis='both', labelsize=9)
    ax3.tick_params(axis='x', rotation=45)
    ax3.legend(fontsize=9, loc='upper left')

    plt.tight_layout()
    out_hm = PLOTS_DIR / f'{stem}__birth_death_heatmap.png'
    plt.savefig(out_hm, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_hm}')


if __name__ == '__main__':
    main()
