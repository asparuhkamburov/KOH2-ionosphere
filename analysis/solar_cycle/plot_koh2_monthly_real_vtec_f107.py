from __future__ import annotations

r"""
KOH2 MONTHLY MEAN VTEC + F10.7 FROM REAL PROCESSED DATA
========================================================

Recommended monograph figure.

INPUT
-----
KOH2_2019_2026_common_hour_values.csv

Supply the validated common-hour CSV with --input-file.

This harmonized file contains actual processed/common-hour values for:
    - PyTECGg VEq
    - IGS Final GIM VTEC at KOH2
    - Madrigal GNSS TEC near KOH2
    - PyIRI climatological VTEC
    - F10.7

SCIENTIFIC AGGREGATION
----------------------
1. Keep only STRICT COMMON HOURS where all four VTEC series are available.
2. Compute one DAILY mean for each VTEC series from those same common hours.
3. Compute the MONTHLY mean from the daily means.
   => every available observation day receives equal weight.
4. F10.7 is reduced to one value per day before monthly averaging.
5. Missing months are NOT interpolated. Lines break across data gaps.
6. Coverage is shown separately:
       - number of common observation days per month
       - number of strict common hours per month

OUTPUTS
-------
Select the output directory with --output-dir. If omitted, MONTHLY_FIGURES
is created beside the input CSV.

    KOH2_monthly_mean_VTEC_F107_strict_common_hours.png
    KOH2_monthly_mean_VTEC_F107_strict_common_hours.pdf
    KOH2_monthly_mean_VTEC_F107_strict_common_hours.svg
    KOH2_monthly_strict_common_hour_statistics.csv

Also creates an optional second figure based on each product's available data:
    KOH2_monthly_mean_VTEC_F107_available_data.png
    KOH2_monthly_available_data_statistics.csv

IMPORTANT
---------
The STRICT COMMON-HOUR figure is recommended for direct inter-product
comparison in the monograph because all VTEC products use identical epochs.

No synthetic/interpolated VTEC values are created.
"""

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


# =============================================================================
# SETTINGS
# =============================================================================

INPUT_FILE = Path("KOH2_2019_2026_common_hour_values.csv")
OUTPUT_DIR = Path("MONTHLY_FIGURES")

START_DATE = "2019-01-01"
END_DATE = "2026-12-31"

# Recommended plot uses equal-day monthly means.
# Optional uncertainty shading = +/- one monthly standard error of daily means.
SHOW_SEM_BANDS = True

DPI = 400

# Solar Cycle 25 maximum period shown only as broad contextual shading.
# Remove or edit these dates if you prefer no highlighted interval.
SHADE_HIGH_ACTIVITY = True
HIGH_ACTIVITY_START = pd.Timestamp("2024-09-01")
HIGH_ACTIVITY_END = pd.Timestamp("2025-12-31")


# =============================================================================
# COLUMN DEFINITIONS
# =============================================================================

SERIES = {
    "pytecgg_veq_tecu": "PyTECGg VEq",
    "igs_vtec_tecu": "IGS Final GIM",
    "madrigal_vtec_tecu": "Madrigal GNSS TEC",
    "iri_vtec_tecu": "PyIRI climatology",
}

REQUIRED = {
    "epoch",
    "f107_sfu",
    *SERIES.keys(),
}


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create KOH2 monthly VTEC + F10.7 products from the validated "
            "common-hour table."
        )
    )

    parser.add_argument(
        "--input-file",
        required=True,
        type=Path,
        help="Validated KOH2 common-hour values CSV.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: MONTHLY_FIGURES beside the input file."
        ),
    )

    parser.add_argument(
        "--start-date",
        default="2019-01-01",
        help="Start date, YYYY-MM-DD. Default: 2019-01-01.",
    )

    parser.add_argument(
        "--end-date",
        default="2026-12-31",
        help="End date, YYYY-MM-DD. Default: 2026-12-31.",
    )

    parser.add_argument(
        "--no-sem-bands",
        action="store_true",
        help="Disable ±1 SEM shading in the strict-common-hour figure.",
    )

    parser.add_argument(
        "--no-high-activity-shading",
        action="store_true",
        help="Disable the contextual Solar Cycle 25 high-activity shading.",
    )

    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global INPUT_FILE, OUTPUT_DIR
    global START_DATE, END_DATE
    global SHOW_SEM_BANDS, SHADE_HIGH_ACTIVITY

    INPUT_FILE = args.input_file.expanduser().resolve()

    OUTPUT_DIR = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else INPUT_FILE.parent / "MONTHLY_FIGURES"
    )

    # Validate date syntax while retaining the original inclusive date-window logic.
    pd.Timestamp(args.start_date)
    pd.Timestamp(args.end_date)

    START_DATE = args.start_date
    END_DATE = args.end_date

    if args.no_sem_bands:
        SHOW_SEM_BANDS = False

    if args.no_high_activity_shading:
        SHADE_HIGH_ACTIVITY = False


# =============================================================================
# LOAD
# =============================================================================

def load_common_hour_data() -> pd.DataFrame:
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            "\nInput file not found:\n"
            f"    {INPUT_FILE}\n\n"
            "Check that the common-hour comparison has already been run."
        )

    df = pd.read_csv(INPUT_FILE)

    missing = REQUIRED - set(df.columns)

    if missing:
        raise RuntimeError(
            "Input file is missing required columns:\n"
            + "\n".join(f"    {c}" for c in sorted(missing))
        )

    df = df.copy()

    df["epoch"] = pd.to_datetime(
        df["epoch"],
        utc=True,
        errors="coerce",
    )

    for col in [
        "f107_sfu",
        *SERIES.keys(),
    ]:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce",
        )

    df = df.dropna(
        subset=["epoch"]
    ).copy()

    df = df[
        (df["epoch"] >= pd.Timestamp(START_DATE, tz="UTC"))
        & (df["epoch"] <= pd.Timestamp(END_DATE, tz="UTC") + pd.Timedelta(days=1))
    ].copy()

    df["date"] = df["epoch"].dt.floor("D")
    df["month"] = df["epoch"].dt.to_period("M").dt.to_timestamp()

    return df.sort_values("epoch").reset_index(drop=True)


# =============================================================================
# STRICT COMMON-HOUR MONTHLY STATISTICS
# =============================================================================

def strict_common_monthly(df: pd.DataFrame):
    value_cols = list(SERIES.keys())

    strict = df.dropna(
        subset=value_cols
    ).copy()

    if strict.empty:
        raise RuntimeError(
            "No strict common hours found where all four VTEC series are valid."
        )

    # Daily means from identical strict common hours.
    daily = (
        strict.groupby("date")
        .agg(
            pytecgg_veq_tecu=("pytecgg_veq_tecu", "mean"),
            igs_vtec_tecu=("igs_vtec_tecu", "mean"),
            madrigal_vtec_tecu=("madrigal_vtec_tecu", "mean"),
            iri_vtec_tecu=("iri_vtec_tecu", "mean"),
            f107_sfu=("f107_sfu", "median"),
            n_common_hours=("epoch", "size"),
        )
        .reset_index()
    )

    daily["month"] = (
        pd.to_datetime(daily["date"])
        .dt.to_period("M")
        .dt.to_timestamp()
    )

    # Equal-day monthly means.
    grouped = daily.groupby("month")

    rows = []

    for month, g in grouped:
        row = {
            "month": month,
            "n_days": int(len(g)),
            "n_common_hours": int(g["n_common_hours"].sum()),
            "median_common_hours_per_day": float(g["n_common_hours"].median()),
            "f107_mean_sfu": float(g["f107_sfu"].mean()),
            "f107_median_sfu": float(g["f107_sfu"].median()),
        }

        for col, label in SERIES.items():
            x = pd.to_numeric(g[col], errors="coerce").dropna()

            row[f"{col}_monthly_mean"] = (
                float(x.mean()) if len(x) else np.nan
            )

            row[f"{col}_monthly_median"] = (
                float(x.median()) if len(x) else np.nan
            )

            row[f"{col}_monthly_std"] = (
                float(x.std(ddof=1)) if len(x) >= 2 else np.nan
            )

            row[f"{col}_monthly_sem"] = (
                float(x.std(ddof=1) / np.sqrt(len(x)))
                if len(x) >= 2
                else np.nan
            )

        rows.append(row)

    monthly = pd.DataFrame(rows)

    # Reindex to every calendar month so missing months become real gaps.
    full_months = pd.date_range(
        pd.Timestamp(START_DATE),
        pd.Timestamp(END_DATE),
        freq="MS",
    )

    monthly = (
        monthly.set_index("month")
        .reindex(full_months)
        .rename_axis("month")
        .reset_index()
    )

    # Keep counts as 0 in months with no observations.
    monthly["n_days"] = monthly["n_days"].fillna(0).astype(int)
    monthly["n_common_hours"] = (
        monthly["n_common_hours"]
        .fillna(0)
        .astype(int)
    )

    return daily, monthly


# =============================================================================
# AVAILABLE-DATA MONTHLY STATISTICS
# =============================================================================

def available_data_monthly(df: pd.DataFrame):
    """
    Alternative descriptive series.

    Each product is averaged from its own valid hours:
      hourly -> daily mean -> monthly equal-day mean.

    Unlike the strict figure, sample epochs can differ among products.
    """
    months = pd.date_range(
        pd.Timestamp(START_DATE),
        pd.Timestamp(END_DATE),
        freq="MS",
    )

    result = pd.DataFrame({"month": months})

    for col in SERIES:
        sub = df.dropna(
            subset=[col]
        )[
            ["date", col]
        ].copy()

        daily = (
            sub.groupby("date")[col]
            .mean()
            .rename("daily_mean")
            .reset_index()
        )

        daily["month"] = (
            pd.to_datetime(daily["date"])
            .dt.to_period("M")
            .dt.to_timestamp()
        )

        m = (
            daily.groupby("month")["daily_mean"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )

        m[f"{col}_monthly_mean"] = m["mean"]

        m[f"{col}_monthly_sem"] = np.where(
            m["count"] >= 2,
            m["std"] / np.sqrt(m["count"]),
            np.nan,
        )

        m[f"{col}_n_days"] = m["count"].astype(int)

        m = m[
            [
                "month",
                f"{col}_monthly_mean",
                f"{col}_monthly_sem",
                f"{col}_n_days",
            ]
        ]

        result = result.merge(
            m,
            on="month",
            how="left",
        )

    # F10.7: exactly one daily value, then equal-day monthly mean.
    f107_daily = (
        df[["date", "f107_sfu"]]
        .dropna()
        .groupby("date")["f107_sfu"]
        .median()
        .reset_index()
    )

    f107_daily["month"] = (
        pd.to_datetime(f107_daily["date"])
        .dt.to_period("M")
        .dt.to_timestamp()
    )

    f107_monthly = (
        f107_daily.groupby("month")["f107_sfu"]
        .agg(["mean", "count"])
        .reset_index()
        .rename(
            columns={
                "mean": "f107_mean_sfu",
                "count": "f107_n_days",
            }
        )
    )

    result = result.merge(
        f107_monthly,
        on="month",
        how="left",
    )

    return result


# =============================================================================
# PLOTTING
# =============================================================================

def format_time_axis(ax):
    ax.set_xlim(
        pd.Timestamp(START_DATE),
        pd.Timestamp(END_DATE),
    )

    ax.xaxis.set_major_locator(
        mdates.YearLocator()
    )

    ax.xaxis.set_major_formatter(
        mdates.DateFormatter("%Y")
    )

    ax.xaxis.set_minor_locator(
        mdates.MonthLocator(bymonth=[4, 7, 10])
    )

    ax.grid(
        True,
        which="major",
        alpha=0.25,
    )


def add_high_activity_shading(ax):
    if not SHADE_HIGH_ACTIVITY:
        return

    ax.axvspan(
        HIGH_ACTIVITY_START,
        HIGH_ACTIVITY_END,
        alpha=0.08,
    )

    ymax = ax.get_ylim()[1]

    ax.text(
        HIGH_ACTIVITY_START
        + (HIGH_ACTIVITY_END - HIGH_ACTIVITY_START) / 2,
        ymax * 0.97,
        "High-activity / Solar Cycle 25 maximum period",
        ha="center",
        va="top",
        fontsize=9,
    )


def plot_strict(monthly: pd.DataFrame):
    fig = plt.figure(
        figsize=(14, 9),
        constrained_layout=True,
    )

    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[4.0, 1.2],
    )

    # -------------------------------------------------------------------------
    # Top: VTEC + F10.7
    # -------------------------------------------------------------------------
    ax = fig.add_subplot(gs[0])

    for col, label in SERIES.items():
        ycol = f"{col}_monthly_mean"

        line = ax.plot(
            monthly["month"],
            monthly[ycol],
            marker="o",
            markersize=4,
            linewidth=1.7,
            label=label,
        )[0]

        if SHOW_SEM_BANDS:
            sem_col = f"{col}_monthly_sem"

            y = pd.to_numeric(
                monthly[ycol],
                errors="coerce",
            ).to_numpy(dtype=float)

            sem = pd.to_numeric(
                monthly[sem_col],
                errors="coerce",
            ).to_numpy(dtype=float)

            good = (
                np.isfinite(y)
                & np.isfinite(sem)
            )

            ax.fill_between(
                monthly["month"],
                y - sem,
                y + sem,
                where=good,
                alpha=0.10,
                color=line.get_color(),
                linewidth=0,
            )

    ax.set_ylabel(
        "Monthly mean VTEC (TECU)"
    )

    ax.set_title(
        "KOH2 monthly mean VTEC and F10.7 during Solar Cycle 25 (2019–2026)\n"
        "Strict common-hour comparison; monthly means are equal-weighted across observation days"
    )

    format_time_axis(ax)
    add_high_activity_shading(ax)

    ax_f = ax.twinx()

    fline = ax_f.plot(
        monthly["month"],
        monthly["f107_mean_sfu"],
        linestyle="--",
        marker="s",
        markersize=4,
        linewidth=1.5,
        label="F10.7",
    )[0]

    ax_f.set_ylabel(
        "Monthly mean F10.7 (sfu)"
    )

    # Combined legend.
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax_f.get_legend_handles_labels()

    ax.legend(
        handles1 + handles2,
        labels1 + labels2,
        loc="upper left",
        ncol=2,
        frameon=True,
    )

    # -------------------------------------------------------------------------
    # Bottom: coverage
    # -------------------------------------------------------------------------
    ax_cov = fig.add_subplot(gs[1], sharex=ax)

    ax_cov.bar(
        monthly["month"],
        monthly["n_days"],
        width=22,
        label="Observation days",
    )

    ax_cov.set_ylabel(
        "Days / month"
    )

    ax_cov.set_xlabel(
        "Year"
    )

    format_time_axis(ax_cov)

    ax_hours = ax_cov.twinx()

    ax_hours.plot(
        monthly["month"],
        monthly["n_common_hours"],
        marker=".",
        linewidth=1.0,
        label="Strict common hours",
    )

    ax_hours.set_ylabel(
        "Common hours / month"
    )

    h1, l1 = ax_cov.get_legend_handles_labels()
    h2, l2 = ax_hours.get_legend_handles_labels()

    ax_cov.legend(
        h1 + h2,
        l1 + l2,
        loc="upper left",
        ncol=2,
        frameon=True,
    )

    fig.text(
        0.01,
        0.005,
        "Note: VTEC monthly values are calculated from actual processed common-hour data only. "
        "No missing months are interpolated. Shaded uncertainty is ±1 SEM of daily means.",
        fontsize=9,
    )

    base = (
        OUTPUT_DIR
        / "KOH2_monthly_mean_VTEC_F107_strict_common_hours"
    )

    fig.savefig(
        base.with_suffix(".png"),
        dpi=DPI,
        bbox_inches="tight",
    )

    fig.savefig(
        base.with_suffix(".pdf"),
        bbox_inches="tight",
    )

    fig.savefig(
        base.with_suffix(".svg"),
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_available(monthly: pd.DataFrame):
    fig, ax = plt.subplots(
        figsize=(14, 7)
    )

    for col, label in SERIES.items():
        ax.plot(
            monthly["month"],
            monthly[f"{col}_monthly_mean"],
            marker="o",
            markersize=4,
            linewidth=1.6,
            label=label,
        )

    ax.set_title(
        "KOH2 monthly mean VTEC and F10.7 (2019–2026)\n"
        "Available-data version: each product uses its own valid observation hours"
    )

    ax.set_ylabel(
        "Monthly mean VTEC (TECU)"
    )

    ax.set_xlabel(
        "Year"
    )

    format_time_axis(ax)
    add_high_activity_shading(ax)

    ax_f = ax.twinx()

    ax_f.plot(
        monthly["month"],
        monthly["f107_mean_sfu"],
        linestyle="--",
        marker="s",
        markersize=4,
        linewidth=1.5,
        label="F10.7",
    )

    ax_f.set_ylabel(
        "Monthly mean F10.7 (sfu)"
    )

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax_f.get_legend_handles_labels()

    ax.legend(
        h1 + h2,
        l1 + l2,
        loc="upper left",
        ncol=2,
    )

    fig.text(
        0.01,
        0.01,
        "Descriptive alternative only: product sampling is not identical. "
        "Use the strict common-hour figure for direct inter-product comparison.",
        fontsize=9,
    )

    fig.tight_layout(
        rect=[0, 0.03, 1, 1]
    )

    fig.savefig(
        OUTPUT_DIR
        / "KOH2_monthly_mean_VTEC_F107_available_data.png",
        dpi=DPI,
        bbox_inches="tight",
    )

    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()
    configure_runtime(args)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "=" * 100
    )

    print(
        "KOH2 MONTHLY MEAN VTEC + F10.7 FROM REAL COMMON-HOUR DATA"
    )

    print(
        "=" * 100
    )

    print(
        "Input:",
        INPUT_FILE,
    )

    print(
        "Output:",
        OUTPUT_DIR,
    )

    print(
        "Date range:",
        START_DATE,
        "to",
        END_DATE,
    )

    print(
        "SEM bands:",
        SHOW_SEM_BANDS,
    )

    print(
        "High-activity shading:",
        SHADE_HIGH_ACTIVITY,
    )

    df = load_common_hour_data()

    print(
        f"Input rows: {len(df):,}"
    )

    strict_daily, strict_monthly = strict_common_monthly(
        df
    )

    available_monthly = available_data_monthly(
        df
    )

    strict_csv = (
        OUTPUT_DIR
        / "KOH2_monthly_strict_common_hour_statistics.csv"
    )

    strict_daily_csv = (
        OUTPUT_DIR
        / "KOH2_daily_strict_common_hour_means.csv"
    )

    available_csv = (
        OUTPUT_DIR
        / "KOH2_monthly_available_data_statistics.csv"
    )

    strict_monthly.to_csv(
        strict_csv,
        index=False,
    )

    strict_daily.to_csv(
        strict_daily_csv,
        index=False,
    )

    available_monthly.to_csv(
        available_csv,
        index=False,
    )

    plot_strict(
        strict_monthly
    )

    plot_available(
        available_monthly
    )

    months_with_data = int(
        (
            strict_monthly["n_days"]
            > 0
        ).sum()
    )

    total_common_days = int(
        (
            strict_daily["date"]
            .nunique()
        )
    )

    total_common_hours = int(
        strict_daily[
            "n_common_hours"
        ].sum()
    )

    print()
    print(
        "STRICT COMMON-HOUR COVERAGE"
    )

    print(
        "-" * 100
    )

    print(
        f"Common observation days: {total_common_days:,}"
    )

    print(
        f"Strict common hours:      {total_common_hours:,}"
    )

    print(
        f"Months with data:         {months_with_data}"
    )

    print()
    print(
        "OUTPUTS"
    )

    print(
        "-" * 100
    )

    for path in [
        strict_csv,
        strict_daily_csv,
        OUTPUT_DIR
        / "KOH2_monthly_mean_VTEC_F107_strict_common_hours.png",
        OUTPUT_DIR
        / "KOH2_monthly_mean_VTEC_F107_strict_common_hours.pdf",
        OUTPUT_DIR
        / "KOH2_monthly_mean_VTEC_F107_strict_common_hours.svg",
        available_csv,
        OUTPUT_DIR
        / "KOH2_monthly_mean_VTEC_F107_available_data.png",
    ]:
        print(
            path
        )

    print()
    print(
        "Recommended monograph figure:"
    )

    print(
        OUTPUT_DIR
        / "KOH2_monthly_mean_VTEC_F107_strict_common_hours.png"
    )


if __name__ == "__main__":
    main()
