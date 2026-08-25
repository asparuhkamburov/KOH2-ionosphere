from __future__ import annotations

r"""
FIGURE SUITE 3 — GEOMAGNETIC DISTURBANCE RESPONSE
=================================================

Uses the existing Solar Cycle daily master and the lagged-bootstrap output.

Creates:
  1) Quiet vs storm boxplots for key validation metrics
  2) Within-year storm-minus-quiet RMSE differences
  3) Lagged Spearman-rho heatmap with bootstrap-robust results marked '*'

No new downloads.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MASTER = Path("KOH2_2019_2026_solar_geomagnetic_daily_master.csv")
LAGGED = Path("KOH2_lagged_geomagnetic_all_tests.csv")
OUT = Path("03_GEOMAGNETIC_RESPONSE")
DPI = 400


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create KOH2 geomagnetic-response figures from validated daily analyses."
    )
    parser.add_argument("--master-file", required=True, type=Path,
                        help="Solar/geomagnetic daily master CSV.")
    parser.add_argument("--lagged-file", required=True, type=Path,
                        help="Lagged-bootstrap all-tests CSV.")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Directory for figures and derived CSV files.")
    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global MASTER, LAGGED, OUT
    MASTER = args.master_file.expanduser().resolve()
    LAGGED = args.lagged_file.expanduser().resolve()
    OUT = args.output_dir.expanduser().resolve()

KEY_METRICS = [
    ("pytecgg_vtec_vs_igs_rmse_tecu", "PyTECGg VTEC vs IGS RMSE"),
    ("pytecgg_veq_vs_igs_rmse_tecu", "PyTECGg VEq vs IGS RMSE"),
    ("pyoasis_vtec_vs_igs_rmse_tecu", "pyOASIS vs IGS RMSE"),
    ("pyoasis_vtec_vs_igs_bias_tecu", "pyOASIS vs IGS bias"),
    ("pytecgg_vtec_vs_madrigal_rmse_tecu", "PyTECGg VTEC vs Madrigal RMSE"),
    ("igs_minus_madrigal_mean_tecu", "IGS − Madrigal"),
]


def plot_storm_quiet(master):
    available = [(c, lab) for c, lab in KEY_METRICS if c in master.columns]
    if not available:
        raise RuntimeError("None of the selected storm/quiet metrics exist.")

    n = len(available)
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(13, 4.2 * nrows), constrained_layout=True
    )
    axes = np.atleast_1d(axes).ravel()

    for ax, (col, label) in zip(axes, available):
        q = pd.to_numeric(
            master.loc[master["activity_class"] == "quiet", col],
            errors="coerce",
        ).dropna()
        s = pd.to_numeric(
            master.loc[master["activity_class"] == "storm", col],
            errors="coerce",
        ).dropna()

        ax.boxplot(
            [q.to_numpy(), s.to_numpy()],
            tick_labels=[f"Quiet\nn={len(q)}", f"Storm\nn={len(s)}"],
            showmeans=True,
        )
        ax.set_title(label)
        ax.set_ylabel("TECU")
        ax.grid(True, axis="y", alpha=0.2)

        if len(q) and len(s):
            delta = s.mean() - q.mean()
            ax.text(
                0.5, 0.97,
                f"Δ mean = {delta:+.2f} TECU",
                transform=ax.transAxes,
                ha="center", va="top", fontsize=9,
            )

    for ax in axes[len(available):]:
        ax.set_axis_off()

    fig.suptitle(
        "KOH2 validation metrics on quiet and geomagnetically disturbed days",
        fontsize=15,
    )

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_quiet_vs_storm_validation_metrics{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_within_year(master):
    metrics = [
        ("pytecgg_vtec_vs_igs_rmse_tecu", "PyTECGg VTEC vs IGS"),
        ("pytecgg_veq_vs_igs_rmse_tecu", "PyTECGg VEq vs IGS"),
        ("pyoasis_vtec_vs_igs_rmse_tecu", "pyOASIS vs IGS"),
        ("pytecgg_vtec_vs_madrigal_rmse_tecu", "PyTECGg VTEC vs Madrigal"),
    ]

    rows = []
    for year, yg in master.groupby("year"):
        for col, label in metrics:
            if col not in yg.columns:
                continue
            q = pd.to_numeric(
                yg.loc[yg["activity_class"] == "quiet", col],
                errors="coerce",
            ).dropna()
            s = pd.to_numeric(
                yg.loc[yg["activity_class"] == "storm", col],
                errors="coerce",
            ).dropna()
            if len(q) and len(s):
                rows.append(
                    {
                        "year": int(year),
                        "metric": label,
                        "delta": float(s.mean() - q.mean()),
                        "n_quiet": len(q),
                        "n_storm": len(s),
                    }
                )

    contrasts = pd.DataFrame(rows)
    contrasts.to_csv(OUT / "KOH2_within_year_storm_minus_quiet_RMSE.csv", index=False)

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)

    for metric, g in contrasts.groupby("metric"):
        g = g.sort_values("year")
        ax.plot(g["year"], g["delta"], marker="o", linewidth=1.6, label=metric)

    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Storm − quiet mean RMSE (TECU)")
    ax.set_title(
        "Within-year storm–quiet contrast in validation RMSE\n"
        "Positive values indicate larger mean RMSE on storm-classified days"
    )
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_within_year_storm_minus_quiet_RMSE{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_lag_heatmap():
    if not LAGGED.is_file():
        raise FileNotFoundError(LAGGED)

    df = pd.read_csv(LAGGED)

    required = {
        "outcome", "predictor", "lag_days_after_activity",
        "spearman_rho", "bootstrap_robust"
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Lagged file missing: {sorted(missing)}")

    wanted = [
        "pytecgg_vtec_vs_igs_rmse_tecu",
        "pytecgg_veq_vs_igs_rmse_tecu",
        "pyoasis_vtec_vs_igs_bias_tecu",
        "pytecgg_vtec_vs_madrigal_rmse_tecu",
        "igs_minus_madrigal_mean_tecu_background_residual",
    ]

    sub = df[df["outcome"].isin(wanted)].copy()
    if sub.empty:
        raise RuntimeError(
            "None of the requested lagged outcomes found.\n"
            "Available examples:\n"
            + "\n".join(df["outcome"].dropna().astype(str).unique()[:30])
        )

    predictor_order = ["kp_max", "Ap_daily", "dst_intensity", "symh_intensity"]
    cols = [(p, lag) for p in predictor_order for lag in [0, 1, 2]]

    row_order = [x for x in wanted if x in set(sub["outcome"])]

    matrix = np.full((len(row_order), len(cols)), np.nan)
    robust = np.zeros_like(matrix, dtype=bool)

    for i, outcome in enumerate(row_order):
        for j, (pred, lag) in enumerate(cols):
            hit = sub[
                (sub["outcome"] == outcome)
                & (sub["predictor"] == pred)
                & (pd.to_numeric(sub["lag_days_after_activity"], errors="coerce") == lag)
            ]
            if len(hit):
                matrix[i, j] = pd.to_numeric(
                    hit.iloc[0]["spearman_rho"], errors="coerce"
                )
                robust[i, j] = bool(hit.iloc[0]["bootstrap_robust"])

    labels = {
        "pytecgg_vtec_vs_igs_rmse_tecu": "PyTECGg VTEC–IGS RMSE",
        "pytecgg_veq_vs_igs_rmse_tecu": "PyTECGg VEq–IGS RMSE",
        "pyoasis_vtec_vs_igs_bias_tecu": "pyOASIS–IGS bias",
        "pytecgg_vtec_vs_madrigal_rmse_tecu": "PyTECGg VTEC–Madrigal RMSE",
        "igs_minus_madrigal_mean_tecu_background_residual":
            "IGS–Madrigal residual",
    }

    fig, ax = plt.subplots(
        figsize=(13, 5.8), constrained_layout=True
    )

    im = ax.imshow(
        matrix,
        aspect="auto",
        vmin=-0.5,
        vmax=0.5,
        cmap="RdBu_r",
    )

    ax.set_yticks(np.arange(len(row_order)))
    ax.set_yticklabels([labels.get(x, x) for x in row_order])

    xlabels = []
    for p, lag in cols:
        name = {
            "kp_max": "Kp",
            "Ap_daily": "Ap",
            "dst_intensity": "−Dst",
            "symh_intensity": "−SYM-H",
        }[p]
        xlabels.append(f"{name}\nlag {lag}")

    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(xlabels)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                mark = "*" if robust[i, j] else ""
                ax.text(
                    j, i, f"{matrix[i, j]:+.2f}{mark}",
                    ha="center", va="center", fontsize=8,
                )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Spearman ρ")

    ax.set_title(
        "Lagged geomagnetic associations after F10.7 + seasonal background removal\n"
        "* = both year-cluster and 7-day block-bootstrap 95% CIs exclude zero"
    )

    for ext in [".png", ".pdf", ".svg"]:
        fig.savefig(
            OUT / f"KOH2_lagged_geomagnetic_robustness_heatmap{ext}",
            dpi=DPI if ext == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def main():
    configure_runtime(parse_args())

    OUT.mkdir(parents=True, exist_ok=True)

    if not MASTER.is_file():
        raise FileNotFoundError(MASTER)

    master = pd.read_csv(MASTER)
    master["year"] = pd.to_numeric(master["year"], errors="coerce")

    plot_storm_quiet(master)
    plot_within_year(master)
    plot_lag_heatmap()

    print("Output:", OUT)


if __name__ == "__main__":
    main()
