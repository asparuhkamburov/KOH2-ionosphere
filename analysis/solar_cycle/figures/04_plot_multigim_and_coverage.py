from __future__ import annotations

r"""
FIGURE SUITE 4 — MULTI-GIM ROBUSTNESS AND DATA COVERAGE
=======================================================

Uses:
  - strict common-day CODE/ESA/JPL/UPC multi-GIM yearly summary
  - KOH2 Solar Cycle daily master

Creates:
  1) PyTECGg and pyOASIS bias against CODE/ESA/JPL/UPC by year
  2) PyTECGg-minus-pyOASIS bias difference by GIM
  3) Year × month heatmap of KOH2 observation-day coverage

The multi-GIM reader is intentionally tolerant of different column names.
It can infer the GIM and method from text columns if necessary.
"""

from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MULTI_ROOT = Path("TEC_VALIDATION_MULTI_GIM_2019_2026")
MASTER = Path("KOH2_2019_2026_solar_geomagnetic_daily_master.csv")
OUT = Path("04_MULTIGIM_AND_COVERAGE")
DPI = 400


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create KOH2 multi-GIM robustness and temporal-coverage figures."
    )
    parser.add_argument("--multigim-root", required=True, type=Path,
                        help="Directory containing the strict common-day multi-GIM yearly summary CSV.")
    parser.add_argument("--master-file", required=True, type=Path,
                        help="Solar/geomagnetic daily master CSV used for coverage.")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Directory for figures and derived CSV files.")
    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global MULTI_ROOT, MASTER, OUT
    MULTI_ROOT = args.multigim_root.expanduser().resolve()
    MASTER = args.master_file.expanduser().resolve()
    OUT = args.output_dir.expanduser().resolve()

GIMS = ["CODE", "ESA", "JPL", "UPC"]


def find_multigim_file():
    candidates = [
        MULTI_ROOT / "KOH2_2019_2026_multiGIM_strict_common_day_yearly_summary.csv",
        *sorted(MULTI_ROOT.glob("*strict*common*day*yearly*summary*.csv")),
    ]
    for p in candidates:
        if p.is_file():
            return p
    raise FileNotFoundError(
        f"No strict common-day multi-GIM summary found in {MULTI_ROOT}"
    )


def infer_text(df, row, needles):
    for col in df.columns:
        val = str(row[col])
        u = val.upper()
        for needle in needles:
            if needle.upper() in u:
                return needle
    return None


def normalize_multigim(df):
    cols_lower = {c.lower(): c for c in df.columns}

    year_col = next(
        (c for c in df.columns if c.lower() == "year"),
        None,
    )

    bias_candidates = [
        "mean_daily_bias_tecu",
        "bias_tecu",
        "median_daily_bias_tecu",
    ]
    bias_col = next(
        (cols_lower[x] for x in bias_candidates if x in cols_lower),
        None,
    )

    if year_col is None or bias_col is None:
        raise RuntimeError(
            "Could not identify year/bias columns.\n"
            f"Columns: {list(df.columns)}"
        )

    ref_col = next(
        (
            c for c in df.columns
            if c.lower() in {"gim", "reference", "reference_product", "product"}
        ),
        None,
    )

    method_col = next(
        (
            c for c in df.columns
            if c.lower() in {"method", "comparison", "tec_method"}
        ),
        None,
    )

    rows = []

    for _, row in df.iterrows():
        year = pd.to_numeric(row[year_col], errors="coerce")
        bias = pd.to_numeric(row[bias_col], errors="coerce")

        if not np.isfinite(year) or not np.isfinite(bias):
            continue

        if ref_col is not None:
            ref_text = str(row[ref_col]).upper()
            gim = next((g for g in GIMS if g in ref_text), None)
        else:
            gim = infer_text(df, row, GIMS)

        if method_col is not None:
            method_text = str(row[method_col]).lower()
        else:
            method_text = " ".join(str(v).lower() for v in row.values)

        if "pytecgg" in method_text:
            method = "PyTECGg VTEC"
        elif "pyoasis" in method_text:
            method = "pyOASIS VTEC"
        else:
            method = infer_text(df, row, ["PyTECGg", "pyOASIS"])
            if method == "PyTECGg":
                method = "PyTECGg VTEC"
            elif method == "pyOASIS":
                method = "pyOASIS VTEC"

        if gim and method:
            rows.append(
                {
                    "year": int(year),
                    "gim": gim,
                    "method": method,
                    "bias_tecu": float(bias),
                }
            )

    out = pd.DataFrame(rows)

    if out.empty:
        raise RuntimeError(
            "Could not infer multi-GIM rows from the CSV.\n"
            f"Columns: {list(df.columns)}"
        )

    return out


def plot_multigim(long):
    fig, axes = plt.subplots(
        1, 2, figsize=(13, 5.2), sharey=True, constrained_layout=True
    )

    for ax, method in zip(axes, ["PyTECGg VTEC", "pyOASIS VTEC"]):
        sub = long[long["method"] == method]
        for gim in GIMS:
            g = sub[sub["gim"] == gim].sort_values("year")
            if not g.empty:
                ax.plot(
                    g["year"], g["bias_tecu"],
                    marker="o", linewidth=1.7, label=gim
                )
        ax.axhline(0, linewidth=0.8)
        ax.set_title(method)
        ax.set_xlabel("Year")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("Mean daily bias vs GIM (TECU)")
    fig.suptitle(
        "KOH2 strict common-day robustness across CODE, ESA, JPL and UPC GIMs",
        fontsize=14,
    )

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_multiGIM_bias_robustness{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_method_offset(long):
    pivot = (
        long.pivot_table(
            index=["year", "gim"],
            columns="method",
            values="bias_tecu",
            aggfunc="mean",
        )
        .reset_index()
    )

    if not {"PyTECGg VTEC", "pyOASIS VTEC"}.issubset(pivot.columns):
        print("Skipping method-offset plot: both methods not available.")
        return

    pivot["PyTECGg_minus_pyOASIS"] = (
        pivot["PyTECGg VTEC"] - pivot["pyOASIS VTEC"]
    )

    pivot.to_csv(
        OUT / "KOH2_multiGIM_PyTECGg_minus_pyOASIS.csv",
        index=False,
    )

    fig, ax = plt.subplots(figsize=(10, 5.8), constrained_layout=True)

    for gim in GIMS:
        g = pivot[pivot["gim"] == gim].sort_values("year")
        if not g.empty:
            ax.plot(
                g["year"], g["PyTECGg_minus_pyOASIS"],
                marker="o", linewidth=1.7, label=gim
            )

    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("PyTECGg bias − pyOASIS bias (TECU)")
    ax.set_title(
        "Inter-method VTEC offset across selected GIM products"
    )
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_multiGIM_inter_method_offset{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_coverage():
    if not MASTER.is_file():
        raise FileNotFoundError(MASTER)

    df = pd.read_csv(MASTER)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.dropna(subset=["date"]).copy()
    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month

    counts = (
        df.groupby(["year", "month"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=range(2019, 2027), columns=range(1, 13), fill_value=0)
    )

    counts.to_csv(OUT / "KOH2_observation_days_by_year_month.csv")

    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)

    im = ax.imshow(counts.to_numpy(), aspect="auto")

    ax.set_yticks(np.arange(len(counts.index)))
    ax.set_yticklabels(counts.index.astype(str))

    month_labels = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ]
    ax.set_xticks(np.arange(12))
    ax.set_xticklabels(month_labels)

    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            value = int(counts.iloc[i, j])
            ax.text(j, i, str(value), ha="center", va="center", fontsize=9)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Observation days")

    ax.set_title(
        "KOH2 temporal coverage of the Solar Cycle 25 analysis\n"
        "Number of available observation days by month and year"
    )
    ax.set_xlabel("Month")
    ax.set_ylabel("Year")

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_observation_coverage_heatmap{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def main():
    configure_runtime(parse_args())

    OUT.mkdir(parents=True, exist_ok=True)

    path = find_multigim_file()
    raw = pd.read_csv(path)
    long = normalize_multigim(raw)

    long.to_csv(OUT / "KOH2_multiGIM_normalized_plot_data.csv", index=False)

    print("Multi-GIM input:", path)
    print("Normalized rows:", len(long))

    plot_multigim(long)
    plot_method_offset(long)
    plot_coverage()

    print("Output:", OUT)


if __name__ == "__main__":
    main()
