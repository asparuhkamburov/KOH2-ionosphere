from __future__ import annotations

r"""
FIGURE SUITE 2 — SOLAR-CYCLE BACKGROUND AND REFERENCE DIVERGENCE
===============================================================

Uses real harmonized common-hour data and an existing Solar Cycle daily
master table supplied at runtime.

Creates:
  1) Monthly equal-day inter-product differences + F10.7
  2) Observed daily VTEC vs fitted F10.7+season background for
     PyTECGg, IGS and Madrigal
  3) Monthly residual statistics CSV

No missing months are interpolated.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

COMMON = Path("KOH2_2019_2026_common_hour_values.csv")
MASTER = Path("KOH2_2019_2026_solar_geomagnetic_daily_master.csv")
MODEL_SUMMARY = Path("KOH2_2019_2026_background_regression_F107_season.csv")
OUT = Path("02_SOLAR_REFERENCE_DIVERGENCE")
DPI = 400


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create KOH2 solar-cycle reference-divergence figures."
    )
    parser.add_argument("--common-hour-file", required=True, type=Path,
                        help="Validated KOH2 common-hour values CSV.")
    parser.add_argument("--master-file", required=True, type=Path,
                        help="Solar/geomagnetic daily master CSV.")
    parser.add_argument("--model-summary-file", type=Path, default=None,
                        help="Optional F10.7+season regression summary CSV used for reported R^2 values.")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Directory for figures and derived CSV files.")
    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global COMMON, MASTER, MODEL_SUMMARY, OUT
    COMMON = args.common_hour_file.expanduser().resolve()
    MASTER = args.master_file.expanduser().resolve()
    MODEL_SUMMARY = (
        args.model_summary_file.expanduser().resolve()
        if args.model_summary_file is not None
        else MASTER.parent / "KOH2_2019_2026_background_regression_F107_season.csv"
    )
    OUT = args.output_dir.expanduser().resolve()

SERIES = [
    "pytecgg_veq_tecu",
    "igs_vtec_tecu",
    "madrigal_vtec_tecu",
    "iri_vtec_tecu",
]


def monthly_common_differences():
    df = pd.read_csv(COMMON)
    df["epoch"] = pd.to_datetime(df["epoch"], utc=True, errors="coerce")
    for c in ["f107_sfu", *SERIES]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    strict = df.dropna(subset=SERIES).copy()
    strict["date"] = strict["epoch"].dt.floor("D")

    daily = (
        strict.groupby("date")
        .agg(
            pytecgg=("pytecgg_veq_tecu", "mean"),
            igs=("igs_vtec_tecu", "mean"),
            madrigal=("madrigal_vtec_tecu", "mean"),
            pyiri=("iri_vtec_tecu", "mean"),
            f107=("f107_sfu", "median"),
            n_hours=("epoch", "size"),
        )
        .reset_index()
    )

    daily["pytecgg_minus_igs"] = daily["pytecgg"] - daily["igs"]
    daily["pytecgg_minus_madrigal"] = daily["pytecgg"] - daily["madrigal"]
    daily["igs_minus_madrigal"] = daily["igs"] - daily["madrigal"]
    daily["month"] = daily["date"].dt.to_period("M").dt.to_timestamp()

    monthly = (
        daily.groupby("month")
        .agg(
            pytecgg_minus_igs=("pytecgg_minus_igs", "mean"),
            pytecgg_minus_madrigal=("pytecgg_minus_madrigal", "mean"),
            igs_minus_madrigal=("igs_minus_madrigal", "mean"),
            f107=("f107", "mean"),
            n_days=("date", "size"),
            n_common_hours=("n_hours", "sum"),
        )
        .reset_index()
    )

    full = pd.DataFrame(
        {"month": pd.date_range("2019-01-01", "2026-12-01", freq="MS")}
    )
    monthly = full.merge(monthly, on="month", how="left")
    return daily, monthly


def plot_monthly_differences(monthly):
    fig = plt.figure(figsize=(14, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1.1])

    ax = fig.add_subplot(gs[0])
    ax.plot(
        monthly["month"], monthly["pytecgg_minus_igs"],
        marker="o", markersize=3.5, linewidth=1.5,
        label="PyTECGg VEq − IGS"
    )
    ax.plot(
        monthly["month"], monthly["pytecgg_minus_madrigal"],
        marker="o", markersize=3.5, linewidth=1.5,
        label="PyTECGg VEq − Madrigal"
    )
    ax.plot(
        monthly["month"], monthly["igs_minus_madrigal"],
        marker="o", markersize=3.5, linewidth=1.5,
        label="IGS − Madrigal"
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_ylabel("Monthly mean difference (TECU)")
    ax.set_title(
        "KOH2 monthly inter-product VTEC differences during Solar Cycle 25\n"
        "Strict common-hour data; monthly values are equal-weighted across observation days"
    )
    ax.grid(True, alpha=0.25)

    ax2 = ax.twinx()
    ax2.plot(
        monthly["month"], monthly["f107"],
        linestyle="--", linewidth=1.5, label="F10.7"
    )
    ax2.set_ylabel("Monthly mean F10.7 (sfu)")

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8)

    ax_cov = fig.add_subplot(gs[1], sharex=ax)
    ax_cov.bar(monthly["month"], monthly["n_days"].fillna(0), width=22)
    ax_cov.set_ylabel("Days")
    ax_cov.set_xlabel("Year")
    ax_cov.grid(True, alpha=0.2)

    for a in [ax, ax_cov]:
        a.xaxis.set_major_locator(mdates.YearLocator())
        a.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_monthly_reference_divergence_F107{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_background_fit():
    master = pd.read_csv(MASTER)

    pairs = [
        (
            "pytecgg_veq_mean_tecu",
            "pytecgg_veq_mean_tecu_background_fitted",
            "PyTECGg VEq",
        ),
        (
            "igs_vtec_mean_tecu",
            "igs_vtec_mean_tecu_background_fitted",
            "IGS Final GIM",
        ),
        (
            "madrigal_vtec_mean_tecu",
            "madrigal_vtec_mean_tecu_background_fitted",
            "Madrigal GNSS TEC",
        ),
    ]

    model_r2 = {}
    if MODEL_SUMMARY.is_file():
        ms = pd.read_csv(MODEL_SUMMARY)
        if {"dependent", "r2"}.issubset(ms.columns):
            model_r2 = dict(zip(ms["dependent"], ms["r2"]))

    fig, axes = plt.subplots(
        1, 3, figsize=(14, 4.8), constrained_layout=True
    )

    for ax, (obs_col, fit_col, label) in zip(axes, pairs):
        if obs_col not in master.columns or fit_col not in master.columns:
            ax.text(0.5, 0.5, f"Missing columns for {label}",
                    ha="center", va="center")
            ax.set_axis_off()
            continue

        x = pd.to_numeric(master[fit_col], errors="coerce")
        y = pd.to_numeric(master[obs_col], errors="coerce")
        good = x.notna() & y.notna()

        xg = x[good].to_numpy(float)
        yg = y[good].to_numpy(float)

        ax.scatter(xg, yg, s=20, alpha=0.65)

        if len(xg):
            lo = min(np.min(xg), np.min(yg))
            hi = max(np.max(xg), np.max(yg))
            ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)

        r2 = model_r2.get(obs_col, np.nan)
        if not np.isfinite(r2) and len(xg) >= 2:
            ss_res = np.sum((yg - xg) ** 2)
            ss_tot = np.sum((yg - np.mean(yg)) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

        ax.set_title(f"{label}\n$R^2$ = {r2:.3f}" if np.isfinite(r2) else label)
        ax.set_xlabel("F10.7 + seasonal background fitted VTEC (TECU)")
        ax.grid(True, alpha=0.2)

    axes[0].set_ylabel("Observed daily mean VTEC (TECU)")
    fig.suptitle(
        "Observed KOH2 VTEC versus fitted F10.7 + seasonal background",
        fontsize=14,
    )

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_observed_vs_F107_season_background_fit{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def main():
    configure_runtime(parse_args())

    OUT.mkdir(parents=True, exist_ok=True)

    if not COMMON.is_file():
        raise FileNotFoundError(COMMON)
    if not MASTER.is_file():
        raise FileNotFoundError(MASTER)

    daily, monthly = monthly_common_differences()
    daily.to_csv(OUT / "KOH2_daily_common_reference_differences.csv", index=False)
    monthly.to_csv(OUT / "KOH2_monthly_common_reference_differences.csv", index=False)

    plot_monthly_differences(monthly)
    plot_background_fit()

    print("Output:", OUT)


if __name__ == "__main__":
    main()
