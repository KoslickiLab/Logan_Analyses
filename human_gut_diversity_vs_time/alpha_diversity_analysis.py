#!/usr/bin/env python3
"""
Alpha Diversity (Observed species/organisms) Yearly Histograms, Violin, and Time Series
-----------------------------------------------------------------------------------------
Assumes files named: alpha_diversity_YYYY.csv for YYYY in [2012..2023]
with columns like:
accession,taxa_coverage_1,taxa_coverage_0_5,...,taxa_coverage_0

Outputs:
- plots/histograms/hist_{coverage}_{year}.png              (per-year histograms per coverage)
- plots/histograms/hist_{coverage}_ALL.png                 (all-years combined histograms per coverage)
- plots/violin/violin_{coverage}.png                       (violin over time for a chosen coverage)
- plots/averages/line_{coverage}_{agg}_{err}.png           (yearly average line plots with error bars)
- summaries/yearly_stats_{coverage}.csv                    (table of yearly stats per coverage and aggregator)
- summaries/available_coverages.csv                        (coverage columns discovered across files)

Edit the CONFIG section below to try different coverages, aggregations, and error metrics.
"""

from __future__ import annotations

import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import gmean, trim_mean, mstats

# -----------------------------
# CONFIG — tweak as you like
# -----------------------------

@dataclass
class Config:
    # Directory containing alpha_diversity_YYYY.csv files
    data_dir: Path = Path("/scratch/dmk333_new/Logan/Logan_Analyses/human_gut_diversity_vs_time/results")

    # Years to include (inclusive)
    years: List[int] = field(default_factory=lambda: list(range(2012, 2024)))

    # Coverage columns to consider.
    # If None, will auto-detect from the union of columns present in data files.
    coverage_columns: List[str] | None = None

    # Coverage to use for violin/time-series plots by default
    default_timeseries_coverage: str | None = "taxa_coverage_0"

    # VIOLIN PLOT SETTINGS
    # Which coverages to draw violins for:
    #   - None  => make a violin figure for EVERY coverage column
    #   - ["taxa_coverage_0", ...] => only those
    violin_coverages: ["taxa_coverage_1", "taxa_coverage_0_5", "taxa_coverage_0_25", "taxa_coverage_0_125", "taxa_coverage_0_0625", "taxa_coverage_0_03125", "taxa_coverage_0_015625", "taxa_coverage_0"]
    # Turn off extrema whiskers to avoid "error-bar" look
    violin_show_extrema: bool = False
    # Keep median line
    violin_show_median: bool = True

    # Aggregations to compute for time series plots
    # Built-ins: "mean", "median", "geometric_mean", "trimmed_mean_pXX", "winsorized_mean_pXX"
    #   where pXX is the proportion to cut/winsorize on each tail (e.g., p10 => 0.10).
    aggregators: List[str] = (
        "mean",
        "median",
        "geometric_mean",
        "trimmed_mean_p10",
        "winsorized_mean_p10",
    )

    # Error bar metric for time series: one of "std", "sem", "iqr", "ci95"
    error_metric: str = "std"

    # Histogram settings
    histogram_bins: str | int = "auto"
    histogram_density: bool = False
    histogram_alpha: float = 0.8

    # Figure size (inches) for all plots
    fig_size: Tuple[float, float] = (8.0, 5.0)

    # DPI for saved figures
    dpi: int = 200

    # Output directories
    out_dir_plots: Path = Path("plots")
    out_dir_summaries: Path = Path("summaries")


CFG = Config()


# -----------------------------
# Helper utilities
# -----------------------------

YEAR_RE = re.compile(r"alpha_diversity_(\d{4})\.csv$")


def discover_years_and_files(data_dir: Path) -> Dict[int, Path]:
    files = {}
    for p in sorted(data_dir.glob("alpha_diversity_*.csv")):
        m = YEAR_RE.search(p.name)
        if m:
            yr = int(m.group(1))
            files[yr] = p
    return files


def read_year_df(path: Path, year: int) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.insert(0, "year", year)
    return df


def discover_coverage_columns(files_by_year: Dict[int, Path]) -> List[str]:
    cols = set()
    for yr, p in files_by_year.items():
        df = pd.read_csv(p, nrows=1)
        for c in df.columns:
            if c.startswith("taxa_coverage_"):
                cols.add(c)
    return sorted(cols, key=lambda c: (
        # sort by numeric coverage descending (1, 0.5, ..., 0), if possible
        -float(c.split("taxa_coverage_")[1].replace("_", "."))
        if c.split("taxa_coverage_")[1].replace("_", ".").replace(".", "", 1).isdigit() else 0.0
    ))


def melt_long(df: pd.DataFrame, coverage_cols: List[str]) -> pd.DataFrame:
    long_df = df.melt(
        id_vars=["year", "accession"],
        value_vars=coverage_cols,
        var_name="coverage",
        value_name="alpha_diversity",
    )
    # ensure numeric
    long_df["alpha_diversity"] = pd.to_numeric(long_df["alpha_diversity"], errors="coerce")
    long_df = long_df.dropna(subset=["alpha_diversity"])
    return long_df


# -----------------------------
# Aggregation and error metrics
# -----------------------------

def parse_trim_or_winsor(name: str) -> float:
    # "trimmed_mean_p10" -> 0.10 ; "winsorized_mean_p05" -> 0.05
    m = re.search(r"_p(\d+)$", name)
    if not m:
        raise ValueError(f"Invalid aggregator name: {name}")
    p = int(m.group(1)) / 100.0
    if not (0 <= p < 0.5):
        raise ValueError("Proportion p should be in [0, 0.5).")
    return p


def agg_fn(name: str) -> Callable[[np.ndarray], float]:
    name = name.lower()
    if name == "mean":
        return lambda x: float(np.mean(x)) if len(x) else float("nan")
    if name == "median":
        return lambda x: float(np.median(x)) if len(x) else float("nan")
    if name == "geometric_mean":
        # Use scipy.stats.gmean on strictly positive values; fall back to NaN if none.
        def gm(x: np.ndarray) -> float:
            x = np.asarray(x, dtype=float)
            x = x[x > 0]
            return float(gmean(x)) if len(x) else float("nan")
        return gm
    if name.startswith("trimmed_mean_p"):
        p = parse_trim_or_winsor(name)
        return lambda x: float(trim_mean(x, proportiontocut=p)) if len(x) else float("nan")
    if name.startswith("winsorized_mean_p"):
        p = parse_trim_or_winsor(name)
        return lambda x: float(np.mean(mstats.winsorize(x, limits=p))) if len(x) else float("nan")
    raise ValueError(f"Unknown aggregator: {name}")


def error_fn(name: str) -> Callable[[np.ndarray], float]:
    name = name.lower()
    if name == "std":
        return lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else float("nan")
    if name == "sem":
        return lambda x: float(np.std(x, ddof=1) / math.sqrt(len(x))) if len(x) > 1 else float("nan")
    if name == "iqr":
        def iqr(x: np.ndarray) -> float:
            if len(x) < 2:
                return float("nan")
            q75, q25 = np.percentile(x, [75, 25])
            return float(q75 - q25)
        return iqr
    if name == "ci95":
        # Symmetric half-width using normal approximation: 1.96 * std / sqrt(n)
        return lambda x: float(1.96 * np.std(x, ddof=1) / math.sqrt(len(x))) if len(x) > 1 else float("nan")
    raise ValueError(f"Unknown error metric: {name}")


# -----------------------------
# Plot helpers
# -----------------------------

def prep_matplotlib(fig_size=(8, 5)):
    # Basic, clean settings (no external styles; avoid specifying colors)
    plt.rcParams.update({
        "figure.figsize": fig_size,
        "axes.grid": True,
        "grid.alpha": 0.4,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.dpi": 200,
    })


def savefig(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    plt.close()


# -----------------------------
# Core analysis
# -----------------------------

def main(cfg: Config = CFG):
    prep_matplotlib(cfg.fig_size)

    files_by_year = discover_years_and_files(cfg.data_dir)
    if not files_by_year:
        print("No alpha_diversity_YYYY.csv files found in", cfg.data_dir, file=sys.stderr)
        sys.exit(1)

    years = cfg.years or sorted(files_by_year.keys())
    # read and concat
    dfs = []
    for yr in years:
        if yr not in files_by_year:
            print(f"Warning: missing file for year {yr}", file=sys.stderr)
            continue
        dfs.append(read_year_df(files_by_year[yr], yr))
    if not dfs:
        print("No data ingested.", file=sys.stderr)
        sys.exit(1)
    df_all = pd.concat(dfs, ignore_index=True)

    # coverage columns
    coverage_cols = cfg.coverage_columns or discover_coverage_columns(files_by_year)
    if not coverage_cols:
        print("No coverage columns found.", file=sys.stderr)
        sys.exit(1)

    # Persist discovered coverage columns for reference
    cfg.out_dir_summaries.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"coverage": coverage_cols}).to_csv(cfg.out_dir_summaries / "available_coverages.csv", index=False)

    # Long form
    long_all = melt_long(df_all, coverage_cols)

    # 1) HISTOGRAMS — per year per coverage + all-years per coverage
    for cov in coverage_cols:
        # Per-year histograms
        for yr in years:
            subset = long_all[(long_all["coverage"] == cov) & (long_all["year"] == yr)]
            if subset.empty:
                continue
            plt.figure()
            plt.hist(
                subset["alpha_diversity"].values,
                bins=cfg.histogram_bins,
                density=cfg.histogram_density,
                alpha=cfg.histogram_alpha,
            )
            plt.title(f"Histogram of {cov} — {yr}")
            plt.xlabel("Observed alpha diversity")
            plt.ylabel("Frequency" if not cfg.histogram_density else "Density")
            plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
            savefig(cfg.out_dir_plots / "histograms" / f"hist_{cov}_{yr}.png")

        # All-years combined histogram
        subset_all = long_all[long_all["coverage"] == cov]
        if not subset_all.empty:
            plt.figure()
            plt.hist(
                subset_all["alpha_diversity"].values,
                bins=cfg.histogram_bins,
                density=cfg.histogram_density,
                alpha=cfg.histogram_alpha,
            )
            plt.title(f"Histogram of {cov} — ALL YEARS")
            plt.xlabel("Observed alpha diversity")
            plt.ylabel("Frequency" if not cfg.histogram_density else "Density")
            plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
            savefig(cfg.out_dir_plots / "histograms" / f"hist_{cov}_ALL.png")

        # 2) VIOLIN PLOTS — alpha diversity over years for selected coverages
        if cfg.violin_coverages is None:
            violin_covs = coverage_cols[:]  # all coverages
        else:
            violin_covs = [c for c in cfg.violin_coverages if c in coverage_cols]
            if not violin_covs:
                print("No valid violin_coverages found; defaulting to all coverage columns.", file=sys.stderr)
                violin_covs = coverage_cols[:]

        for violin_cov in violin_covs:
            viol_df = long_all[long_all["coverage"] == violin_cov]
            if viol_df.empty:
                continue
            data_by_year = [viol_df[viol_df["year"] == yr]["alpha_diversity"].values for yr in years]
            plt.figure()
            plt.violinplot(
                data_by_year,
                showmeans=False,
                showmedians=cfg.violin_show_median,
                showextrema=cfg.violin_show_extrema,
                widths=0.9,
            )
            plt.title(f"Alpha Diversity over Time — Violin ({violin_cov})")
            plt.xlabel("Year")
            plt.ylabel("Observed alpha diversity")
            plt.xticks(ticks=range(1, len(years) + 1), labels=[str(y) for y in years], rotation=45, ha="right")
            savefig(cfg.out_dir_plots / "violin" / f"violin_{violin_cov}.png")

    # 3) TIME SERIES — average over years with error bars per aggregator
    err = error_fn(cfg.error_metric)

    # compute summary table: for each coverage, for each year, for each aggregator
    for cov in coverage_cols:
        cov_df = long_all[long_all["coverage"] == cov].copy()
        if cov_df.empty:
            continue
        rows = []
        for yr, grp in cov_df.groupby("year"):
            x = grp["alpha_diversity"].to_numpy(dtype=float)
            n = len(x)
            for agg_name in cfg.aggregators:
                a_fn = agg_fn(agg_name)
                avg = a_fn(x)
                eb = err(x)
                rows.append({"year": int(yr), "coverage": cov, "n": n, "aggregator": agg_name, "value": avg, "error": eb})
        if rows:
            summary = pd.DataFrame(rows).sort_values(["aggregator", "year"])
            # Save table
            out_csv = cfg.out_dir_summaries / f"yearly_stats_{cov}.csv"
            summary.to_csv(out_csv, index=False)

            # Plot each aggregator as a separate line plot (one figure per aggregator to avoid clutter)
            for agg_name, grp in summary.groupby("aggregator"):
                plt.figure()
                # Keep year order consistent
                grp = grp.sort_values("year")
                plt.errorbar(grp["year"], grp["value"], yerr=grp["error"], fmt="-o", capsize=3)
                plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
                plt.title(f"Yearly {agg_name.replace('_', ' ').title()} of {cov}\n(Error: {CFG.error_metric.upper()})")
                plt.xlabel("Year")
                plt.ylabel("Alpha diversity")
                savefig(cfg.out_dir_plots / "averages" / f"line_{cov}_{agg_name}_{CFG.error_metric}.png")

    print("Done. See ./plots/ and ./summaries/")


if __name__ == "__main__":
    main()
