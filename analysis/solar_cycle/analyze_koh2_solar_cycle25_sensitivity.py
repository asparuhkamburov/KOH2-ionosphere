from __future__ import annotations

r"""
KOH2 SOLAR CYCLE 25 SENSITIVITY ANALYSIS
Equal-year weighting + storm/quiet validation-error analysis
==========================================================================

PURPOSE
-------
Refine the already completed Solar Cycle 25 analysis without downloading
any additional solar/geomagnetic data.

INPUT
-----

    KOH2_2019_2026_solar_geomagnetic_daily_master.csv

PART A — EQUAL-YEAR WEIGHTED BACKGROUND REGRESSION
--------------------------------------------------
The prior ordinary daily OLS gives every day equal weight. Because 2024 has
many more KOH2 observation days than other years, it can dominate the fit.

For every dependent variable, this script assigns:

    weight_i = 1 / N_valid_days_in_that_year

so every year represented in that model contributes the same total weight.

Model:
    y = b0
        + b1 * F10.7
        + annual sin/cos
        + semiannual sin/cos

The script reports:
    - ordinary OLS result
    - equal-year weighted WLS result
    - F10.7 slope
    - standard error / p-value
    - R^2
    - RMSE
    - change in F10.7 slope and R^2

Weighted R^2 and weighted RMSE are used for WLS.

PART B — STORM VS QUIET VALIDATION-ERROR ANALYSIS
--------------------------------------------------
Uses the activity classification already stored in the master file:

    quiet
    active_nonstorm
    storm
    unknown

For validation bias/RMSE metrics and key common-hour differences, calculate:

POOLED:
    - quiet n, mean, median
    - storm n, mean, median
    - storm - quiet mean difference
    - storm - quiet median difference
    - Welch t-test
    - Mann-Whitney U test
    - BH-FDR q-values

WITHIN-YEAR:
    For each year having both quiet and storm observations:
        mean(metric | storm) - mean(metric | quiet)
        median(metric | storm) - median(metric | quiet)

    Then summarize those annual contrasts with equal weight per year.

This within-year comparison is important because the distribution of storm
days differs strongly among years and the TEC background itself changes with
Solar Cycle 25.

No empirical TEC correction is applied.

OUTPUT
------
SENSITIVITY

    KOH2_equal_year_weighted_background_models.csv
    KOH2_equal_year_weighted_background_coefficients.csv
    KOH2_OLS_vs_equal_year_WLS.csv

    KOH2_storm_vs_quiet_pooled.csv
    KOH2_storm_vs_quiet_within_year.csv
    KOH2_storm_vs_quiet_equal_year_summary.csv

    KOH2_sensitivity_report.txt
"""

from pathlib import Path
import argparse
import math

import numpy as np
import pandas as pd

try:
    from scipy import stats as scipy_stats
except ImportError as exc:
    raise RuntimeError(
        "\nSciPy is required.\n"
        "Activate pytecgg_env and run:\n\n"
        "    python -m pip install scipy\n"
    ) from exc


# =============================================================================
# SETTINGS
# =============================================================================

MASTER_FILE = Path("KOH2_2019_2026_solar_geomagnetic_daily_master.csv")
OUTPUT_ROOT = Path("SENSITIVITY")

MIN_MODEL_N = 12
MIN_GROUP_N = 2

PREDICTORS = [
    "omni_f107_sfu",
    "season_sin1",
    "season_cos1",
    "season_sin2",
    "season_cos2",
]

BACKGROUND_DEPENDENTS = [
    "pytecgg_veq_mean_tecu",
    "igs_vtec_mean_tecu",
    "madrigal_vtec_mean_tecu",
    "iri_vtec_mean_tecu",
    "pytecgg_minus_igs_mean_tecu",
    "pytecgg_minus_madrigal_mean_tecu",
    "igs_minus_madrigal_mean_tecu",
    "pytecgg_minus_pyiri_mean_tecu",
    "igs_minus_pyiri_mean_tecu",
    "madrigal_minus_pyiri_mean_tecu",
]


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run KOH2 Solar Cycle 25 sensitivity analyses from the daily master table."
    )
    parser.add_argument("--master-file", required=True, type=Path,
                        help="KOH2 solar/geomagnetic daily master CSV.")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Output directory. Default: SENSITIVITY beside the master file.")
    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global MASTER_FILE, OUTPUT_ROOT
    MASTER_FILE = args.master_file.expanduser().resolve()
    OUTPUT_ROOT = (args.output_dir.expanduser().resolve()
                   if args.output_dir is not None
                   else MASTER_FILE.parent / "SENSITIVITY")


# =============================================================================
# HELPERS
# =============================================================================

def add_bh_fdr(
    table: pd.DataFrame,
    p_col: str,
    q_col: str,
):
    out = table.copy()

    out[
        q_col
    ] = np.nan

    if (
        out.empty
        or p_col not in out.columns
    ):
        return out

    p = pd.to_numeric(
        out[
            p_col
        ],
        errors="coerce",
    )

    good = p.notna()

    if not good.any():
        return out

    vals = p[
        good
    ].to_numpy(
        dtype=float
    )

    order = np.argsort(
        vals
    )

    ranked = vals[
        order
    ]

    m = len(
        ranked
    )

    q_ranked = (
        ranked
        * m
        / np.arange(
            1,
            m + 1,
            dtype=float,
        )
    )

    q_ranked = np.minimum.accumulate(
        q_ranked[
            ::-1
        ]
    )[
        ::-1
    ]

    q_ranked = np.clip(
        q_ranked,
        0.0,
        1.0,
    )

    q_vals = np.empty(
        m,
        dtype=float,
    )

    q_vals[
        order
    ] = q_ranked

    out.loc[
        good,
        q_col,
    ] = q_vals

    return out


def weighted_mean(
    x,
    w,
):
    x = np.asarray(
        x,
        dtype=float,
    )

    w = np.asarray(
        w,
        dtype=float,
    )

    good = (
        np.isfinite(
            x
        )
        & np.isfinite(
            w
        )
        & (
            w
            > 0
        )
    )

    if not np.any(
        good
    ):
        return np.nan

    return float(
        np.sum(
            w[
                good
            ]
            * x[
                good
            ]
        )
        / np.sum(
            w[
                good
            ]
        )
    )


# =============================================================================
# REGRESSION
# =============================================================================

def fit_linear_model(
    df: pd.DataFrame,
    y_col: str,
    weighted: bool,
):
    required = [
        "year",
        y_col,
    ] + PREDICTORS

    g = df[
        required
    ].copy()

    for col in [
        y_col,
    ] + PREDICTORS:
        g[
            col
        ] = pd.to_numeric(
            g[
                col
            ],
            errors="coerce",
        )

    g[
        "year"
    ] = pd.to_numeric(
        g[
            "year"
        ],
        errors="coerce",
    )

    g = g.dropna(
        subset=required
    )

    n = len(
        g
    )

    p = (
        len(
            PREDICTORS
        )
        + 1
    )

    if n < max(
        MIN_MODEL_N,
        p + 2,
    ):
        return None

    y = g[
        y_col
    ].to_numpy(
        dtype=float
    )

    X_raw = g[
        PREDICTORS
    ].to_numpy(
        dtype=float
    )

    X = np.column_stack([
        np.ones(
            n
        ),
        X_raw,
    ])

    if weighted:
        counts = g.groupby(
            "year"
        )[
            y_col
        ].transform(
            "size"
        ).to_numpy(
            dtype=float
        )

        w = (
            1.0
            / counts
        )

        # Normalize weights to mean 1. This does not alter beta estimates,
        # but keeps residual variance scaling interpretable.
        w = (
            w
            * n
            / np.sum(
                w
            )
        )
    else:
        w = np.ones(
            n,
            dtype=float,
        )

    sqrt_w = np.sqrt(
        w
    )

    Xw = (
        X
        * sqrt_w[
            :,
            None
        ]
    )

    yw = (
        y
        * sqrt_w
    )

    beta, _, _, _ = np.linalg.lstsq(
        Xw,
        yw,
        rcond=None,
    )

    fitted = (
        X
        @ beta
    )

    residual = (
        y
        - fitted
    )

    ybar_w = weighted_mean(
        y,
        w,
    )

    sse_w = float(
        np.sum(
            w
            * residual
            ** 2
        )
    )

    sst_w = float(
        np.sum(
            w
            * (
                y
                - ybar_w
            )
            ** 2
        )
    )

    r2_w = (
        1.0
        - sse_w
        / sst_w
        if sst_w > 0
        else np.nan
    )

    rmse_w = float(
        np.sqrt(
            np.sum(
                w
                * residual
                ** 2
            )
            / np.sum(
                w
            )
        )
    )

    dof = (
        n
        - p
    )

    sigma2 = (
        sse_w
        / dof
        if dof > 0
        else np.nan
    )

    try:
        xtwx_inv = np.linalg.inv(
            X.T
            @ (
                w[
                    :,
                    None
                ]
                * X
            )
        )

        cov = (
            sigma2
            * xtwx_inv
        )

        se = np.sqrt(
            np.diag(
                cov
            )
        )

    except np.linalg.LinAlgError:
        se = np.full(
            p,
            np.nan,
        )

    with np.errstate(
        divide="ignore",
        invalid="ignore",
    ):
        t_stat = (
            beta
            / se
        )

    if dof > 0:
        p_vals = (
            2.0
            * scipy_stats.t.sf(
                np.abs(
                    t_stat
                ),
                df=dof,
            )
        )
    else:
        p_vals = np.full(
            p,
            np.nan,
        )

    names = [
        "intercept",
    ] + PREDICTORS

    coefficient_rows = []

    for i, name in enumerate(
        names
    ):
        coefficient_rows.append({
            "term":
                name,
            "coefficient":
                float(
                    beta[
                        i
                    ]
                ),
            "std_error":
                float(
                    se[
                        i
                    ]
                )
                if np.isfinite(
                    se[
                        i
                    ]
                )
                else np.nan,
            "t_stat":
                float(
                    t_stat[
                        i
                    ]
                )
                if np.isfinite(
                    t_stat[
                        i
                    ]
                )
                else np.nan,
            "p_value":
                float(
                    p_vals[
                        i
                    ]
                )
                if np.isfinite(
                    p_vals[
                        i
                    ]
                )
                else np.nan,
        })

    return {
        "n":
            n,
        "n_years":
            int(
                g[
                    "year"
                ].nunique()
            ),
        "r2":
            r2_w,
        "rmse":
            rmse_w,
        "coefficients":
            coefficient_rows,
    }


def run_regression_sensitivity(
    master: pd.DataFrame,
):
    model_rows = []
    coefficient_rows = []
    comparison_rows = []

    for y_col in BACKGROUND_DEPENDENTS:
        if y_col not in master.columns:
            continue

        ols = fit_linear_model(
            master,
            y_col,
            weighted=False,
        )

        wls = fit_linear_model(
            master,
            y_col,
            weighted=True,
        )

        if ols is None or wls is None:
            continue

        for label, fit in [
            (
                "ordinary_daily_OLS",
                ols,
            ),
            (
                "equal_year_WLS",
                wls,
            ),
        ]:
            model_rows.append({
                "dependent":
                    y_col,
                "model_type":
                    label,
                "n":
                    fit[
                        "n"
                    ],
                "n_years":
                    fit[
                        "n_years"
                    ],
                "r2":
                    fit[
                        "r2"
                    ],
                "rmse_tecu":
                    fit[
                        "rmse"
                    ],
            })

            for row in fit[
                "coefficients"
            ]:
                coefficient_rows.append({
                    "dependent":
                        y_col,
                    "model_type":
                        label,
                    **row,
                })

        ols_coef = {
            row[
                "term"
            ]:
                row
            for row in ols[
                "coefficients"
            ]
        }

        wls_coef = {
            row[
                "term"
            ]:
                row
            for row in wls[
                "coefficients"
            ]
        }

        b_ols = ols_coef[
            "omni_f107_sfu"
        ][
            "coefficient"
        ]

        b_wls = wls_coef[
            "omni_f107_sfu"
        ][
            "coefficient"
        ]

        comparison_rows.append({
            "dependent":
                y_col,
            "n":
                ols[
                    "n"
                ],
            "n_years":
                ols[
                    "n_years"
                ],
            "ols_r2":
                ols[
                    "r2"
                ],
            "equal_year_wls_r2":
                wls[
                    "r2"
                ],
            "delta_r2_wls_minus_ols":
                (
                    wls[
                        "r2"
                    ]
                    - ols[
                        "r2"
                    ]
                ),
            "ols_rmse_tecu":
                ols[
                    "rmse"
                ],
            "equal_year_wls_rmse_tecu":
                wls[
                    "rmse"
                ],
            "ols_f107_slope_tecu_per_sfu":
                b_ols,
            "equal_year_wls_f107_slope_tecu_per_sfu":
                b_wls,
            "f107_slope_ratio_wls_to_ols":
                (
                    b_wls
                    / b_ols
                    if b_ols != 0
                    else np.nan
                ),
            "ols_f107_p":
                ols_coef[
                    "omni_f107_sfu"
                ][
                    "p_value"
                ],
            "equal_year_wls_f107_p":
                wls_coef[
                    "omni_f107_sfu"
                ][
                    "p_value"
                ],
        })

    return (
        pd.DataFrame(
            model_rows
        ),
        pd.DataFrame(
            coefficient_rows
        ),
        pd.DataFrame(
            comparison_rows
        ),
    )


# =============================================================================
# STORM / QUIET METRIC SELECTION
# =============================================================================

def select_storm_quiet_metrics(
    master: pd.DataFrame,
):
    metrics = []

    # Existing validation error fields from IGS and Madrigal comparisons.
    for col in master.columns:
        lower = col.lower()

        if (
            (
                lower.endswith(
                    "_bias_tecu"
                )
                or lower.endswith(
                    "_rmse_tecu"
                )
            )
            and (
                "vs_igs"
                in lower
                or "vs_madrigal"
                in lower
            )
        ):
            metrics.append(
                col
            )

    # Key harmonized common-hour inter-product differences.
    preferred = [
        "pytecgg_minus_igs_mean_tecu",
        "pytecgg_minus_madrigal_mean_tecu",
        "igs_minus_madrigal_mean_tecu",
        "pytecgg_minus_pyiri_mean_tecu",
        "igs_minus_pyiri_mean_tecu",
        "madrigal_minus_pyiri_mean_tecu",
    ]

    for col in preferred:
        if (
            col in master.columns
            and col not in metrics
        ):
            metrics.append(
                col
            )

    return metrics


# =============================================================================
# POOLED STORM VS QUIET
# =============================================================================

def pooled_storm_quiet(
    master: pd.DataFrame,
    metrics,
):
    rows = []

    quiet = master[
        master[
            "activity_class"
        ]
        == "quiet"
    ]

    storm = master[
        master[
            "activity_class"
        ]
        == "storm"
    ]

    for metric in metrics:
        q = pd.to_numeric(
            quiet[
                metric
            ],
            errors="coerce",
        ).dropna()

        s = pd.to_numeric(
            storm[
                metric
            ],
            errors="coerce",
        ).dropna()

        row = {
            "metric":
                metric,
            "quiet_n":
                len(
                    q
                ),
            "storm_n":
                len(
                    s
                ),
            "quiet_mean":
                q.mean(),
            "storm_mean":
                s.mean(),
            "storm_minus_quiet_mean":
                (
                    s.mean()
                    - q.mean()
                    if (
                        len(
                            q
                        )
                        and len(
                            s
                        )
                    )
                    else np.nan
                ),
            "quiet_median":
                q.median(),
            "storm_median":
                s.median(),
            "storm_minus_quiet_median":
                (
                    s.median()
                    - q.median()
                    if (
                        len(
                            q
                        )
                        and len(
                            s
                        )
                    )
                    else np.nan
                ),
            "quiet_std":
                q.std(
                    ddof=0
                ),
            "storm_std":
                s.std(
                    ddof=0
                ),
            "welch_t_p":
                np.nan,
            "mannwhitney_p":
                np.nan,
        }

        if (
            len(
                q
            )
            >= MIN_GROUP_N
            and len(
                s
            )
            >= MIN_GROUP_N
        ):
            try:
                t = scipy_stats.ttest_ind(
                    s,
                    q,
                    equal_var=False,
                    nan_policy="omit",
                )

                row[
                    "welch_t_p"
                ] = float(
                    t.pvalue
                )
            except Exception:
                pass

            try:
                mw = scipy_stats.mannwhitneyu(
                    s,
                    q,
                    alternative="two-sided",
                )

                row[
                    "mannwhitney_p"
                ] = float(
                    mw.pvalue
                )
            except Exception:
                pass

        rows.append(
            row
        )

    out = pd.DataFrame(
        rows
    )

    out = add_bh_fdr(
        out,
        "welch_t_p",
        "welch_t_q_bh",
    )

    out = add_bh_fdr(
        out,
        "mannwhitney_p",
        "mannwhitney_q_bh",
    )

    return out


# =============================================================================
# WITHIN-YEAR STORM VS QUIET
# =============================================================================

def within_year_storm_quiet(
    master: pd.DataFrame,
    metrics,
):
    rows = []

    for year, g in master.groupby(
        "year"
    ):
        qg = g[
            g[
                "activity_class"
            ]
            == "quiet"
        ]

        sg = g[
            g[
                "activity_class"
            ]
            == "storm"
        ]

        for metric in metrics:
            q = pd.to_numeric(
                qg[
                    metric
                ],
                errors="coerce",
            ).dropna()

            s = pd.to_numeric(
                sg[
                    metric
                ],
                errors="coerce",
            ).dropna()

            if (
                len(
                    q
                )
                < 1
                or len(
                    s
                )
                < 1
            ):
                continue

            rows.append({
                "year":
                    int(
                        year
                    ),
                "metric":
                    metric,
                "quiet_n":
                    len(
                        q
                    ),
                "storm_n":
                    len(
                        s
                    ),
                "quiet_mean":
                    q.mean(),
                "storm_mean":
                    s.mean(),
                "storm_minus_quiet_mean":
                    (
                        s.mean()
                        - q.mean()
                    ),
                "quiet_median":
                    q.median(),
                "storm_median":
                    s.median(),
                "storm_minus_quiet_median":
                    (
                        s.median()
                        - q.median()
                    ),
            })

    return pd.DataFrame(
        rows
    )


def summarize_equal_year_contrasts(
    within_year: pd.DataFrame,
):
    if within_year.empty:
        return pd.DataFrame()

    rows = []

    for metric, g in within_year.groupby(
        "metric"
    ):
        mean_diff = pd.to_numeric(
            g[
                "storm_minus_quiet_mean"
            ],
            errors="coerce",
        ).dropna()

        median_diff = pd.to_numeric(
            g[
                "storm_minus_quiet_median"
            ],
            errors="coerce",
        ).dropna()

        row = {
            "metric":
                metric,
            "n_years_with_both_classes":
                int(
                    g[
                        "year"
                    ].nunique()
                ),
            "years":
                ",".join(
                    str(
                        int(
                            y
                        )
                    )
                    for y in sorted(
                        g[
                            "year"
                        ].unique()
                    )
                ),
            "equal_year_mean_of_mean_differences":
                mean_diff.mean(),
            "median_of_yearly_mean_differences":
                mean_diff.median(),
            "equal_year_mean_of_median_differences":
                median_diff.mean(),
            "median_of_yearly_median_differences":
                median_diff.median(),
            "n_years_positive_mean_difference":
                int(
                    (
                        mean_diff
                        > 0
                    ).sum()
                ),
            "n_years_negative_mean_difference":
                int(
                    (
                        mean_diff
                        < 0
                    ).sum()
                ),
            "sign_test_p":
                np.nan,
        }

        # Exact two-sided binomial sign test across independent yearly contrasts.
        nonzero = mean_diff[
            mean_diff
            != 0
        ]

        if len(
            nonzero
        ) > 0:
            positives = int(
                (
                    nonzero
                    > 0
                ).sum()
            )

            try:
                bt = scipy_stats.binomtest(
                    positives,
                    n=len(
                        nonzero
                    ),
                    p=0.5,
                    alternative="two-sided",
                )

                row[
                    "sign_test_p"
                ] = float(
                    bt.pvalue
                )
            except Exception:
                pass

        rows.append(
            row
        )

    out = pd.DataFrame(
        rows
    )

    out = add_bh_fdr(
        out,
        "sign_test_p",
        "sign_test_q_bh",
    )

    return out


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()
    configure_runtime(args)

    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    pd.set_option(
        "display.max_columns",
        None,
    )

    pd.set_option(
        "display.width",
        300,
    )

    if not MASTER_FILE.is_file():
        raise FileNotFoundError(
            f"Missing daily master file:\n{MASTER_FILE}"
        )

    print(
        "=" * 120
    )

    print(
        "KOH2 SOLAR CYCLE 25 SENSITIVITY ANALYSIS"
    )

    print(
        "=" * 120
    )

    print(
        "Input:",
        MASTER_FILE,
    )

    master = pd.read_csv(
        MASTER_FILE
    )

    print(
        "Daily rows:",
        len(
            master
        ),
    )

    if "activity_class" not in master.columns:
        raise RuntimeError(
            "Master file does not contain activity_class. "
            "Run the previous Solar Cycle 25 analysis first."
        )

    # -------------------------------------------------------------------------
    # A. Equal-year weighted regression
    # -------------------------------------------------------------------------

    (
        model_table,
        coefficient_table,
        comparison_table,
    ) = run_regression_sensitivity(
        master
    )

    # -------------------------------------------------------------------------
    # B. Storm vs quiet
    # -------------------------------------------------------------------------

    metrics = select_storm_quiet_metrics(
        master
    )

    print(
        "Storm/quiet metrics found:",
        len(
            metrics
        ),
    )

    pooled = pooled_storm_quiet(
        master,
        metrics,
    )

    within_year = within_year_storm_quiet(
        master,
        metrics,
    )

    equal_year = summarize_equal_year_contrasts(
        within_year
    )

    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------

    model_file = (
        OUTPUT_ROOT
        / "KOH2_equal_year_weighted_background_models.csv"
    )

    coefficient_file = (
        OUTPUT_ROOT
        / "KOH2_equal_year_weighted_background_coefficients.csv"
    )

    comparison_file = (
        OUTPUT_ROOT
        / "KOH2_OLS_vs_equal_year_WLS.csv"
    )

    pooled_file = (
        OUTPUT_ROOT
        / "KOH2_storm_vs_quiet_pooled.csv"
    )

    within_year_file = (
        OUTPUT_ROOT
        / "KOH2_storm_vs_quiet_within_year.csv"
    )

    equal_year_file = (
        OUTPUT_ROOT
        / "KOH2_storm_vs_quiet_equal_year_summary.csv"
    )

    report_file = (
        OUTPUT_ROOT
        / "KOH2_sensitivity_report.txt"
    )

    model_table.to_csv(
        model_file,
        index=False,
    )

    coefficient_table.to_csv(
        coefficient_file,
        index=False,
    )

    comparison_table.to_csv(
        comparison_file,
        index=False,
    )

    pooled.to_csv(
        pooled_file,
        index=False,
    )

    within_year.to_csv(
        within_year_file,
        index=False,
    )

    equal_year.to_csv(
        equal_year_file,
        index=False,
    )

    with open(
        report_file,
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            "KOH2 SOLAR CYCLE 25 SENSITIVITY ANALYSIS\n"
        )

        f.write(
            "=" * 100
            + "\n\n"
        )

        f.write(
            "PART A — OLS VS EQUAL-YEAR WLS\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            comparison_table.to_string(
                index=False
            )
        )

        f.write(
            "\n\nPART B — POOLED STORM VS QUIET\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            pooled.to_string(
                index=False
            )
        )

        f.write(
            "\n\nPART C — EQUAL-YEAR STORM VS QUIET CONTRASTS\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            equal_year.to_string(
                index=False
            )
        )

        f.write(
            "\n"
        )

    # -------------------------------------------------------------------------
    # Console
    # -------------------------------------------------------------------------

    print()
    print(
        "=" * 120
    )

    print(
        "SENSITIVITY ANALYSIS COMPLETE"
    )

    print(
        "=" * 120
    )

    print(
        "OLS vs equal-year WLS:",
        comparison_file,
    )

    print(
        "Weighted coefficients:",
        coefficient_file,
    )

    print(
        "Pooled storm vs quiet:",
        pooled_file,
    )

    print(
        "Within-year storm vs quiet:",
        within_year_file,
    )

    print(
        "Equal-year storm vs quiet:",
        equal_year_file,
    )

    print(
        "Report:",
        report_file,
    )

    print()
    print(
        "OLS VS EQUAL-YEAR WLS"
    )

    print(
        "-" * 120
    )

    print(
        comparison_table.to_string(
            index=False
        )
    )

    print()
    print(
        "POOLED STORM VS QUIET — strongest absolute mean differences"
    )

    print(
        "-" * 120
    )

    if not pooled.empty:
        temp = pooled.copy()

        temp[
            "_absdiff"
        ] = np.abs(
            pd.to_numeric(
                temp[
                    "storm_minus_quiet_mean"
                ],
                errors="coerce",
            )
        )

        print(
            temp.sort_values(
                "_absdiff",
                ascending=False,
            ).drop(
                columns=[
                    "_absdiff",
                ]
            ).head(
                20
            ).to_string(
                index=False
            )
        )

    print()
    print(
        "EQUAL-YEAR STORM VS QUIET CONTRASTS"
    )

    print(
        "-" * 120
    )

    print(
        equal_year.to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
