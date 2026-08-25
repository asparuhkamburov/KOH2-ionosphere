from __future__ import annotations

r"""
KOH2 SOLAR CYCLE 25 — LAGGED GEOMAGNETIC ROBUSTNESS ANALYSIS
Lag 0/+1/+2 days + year-cluster bootstrap + 7-day moving-block bootstrap
=============================================================================

PURPOSE
-------
Final robustness test for the KOH2 Solar Cycle 25 analysis.

This script asks:

    After the F10.7 + seasonal background has already been removed,
    are the remaining VTEC anomalies / inter-product differences related
    to geomagnetic activity on:

        lag 0 : the SAME calendar day
        lag +1: ONE DAY EARLIER
        lag +2: TWO DAYS EARLIER

The same lag analysis is also applied to the daily validation bias/RMSE
metrics from PyTECGg and pyOASIS.

IMPORTANT LAG CONVENTION
------------------------
For a TEC/error response dated D:

    lag 0 predictor = geomagnetic index on D
    lag 1 predictor = geomagnetic index on D - 1 calendar day
    lag 2 predictor = geomagnetic index on D - 2 calendar days

This uses true CALENDAR-DAY lags, not previous available KOH2 observation rows.

UNCERTAINTY
-----------
Two bootstrap confidence intervals are calculated for the equal-year
weighted simple slope:

A. YEAR-CLUSTER BOOTSTRAP
   - resample complete years with replacement
   - every selected year-cluster receives equal total weight
   - preserves all within-year dependence

B. 7-DAY CALENDAR MOVING-BLOCK BOOTSTRAP
   - resample circular 7-day calendar blocks separately within each year
   - preserves short-term autocorrelation / storm sequences
   - KOH2 observation gaps are preserved because blocks are sampled in
     calendar time, not in observation-row order
   - every year receives equal total weight in each replicate

An association is labelled:

    bootstrap_robust = True

only if BOTH 95% bootstrap CIs exclude zero and have the same sign.

INPUTS
------
Main daily analysis table:

    
        KOH2_2019_2026_solar_geomagnetic_daily_master.csv

Existing cached geomagnetic data:

    ...\_INDEX_CACHE\OMNI2\omni2_YYYY.dat
    ...\_INDEX_CACHE\OMNI_HRO_SYMH\symh_5min_YYYY.csv

NO INTERNET DOWNLOADS ARE PERFORMED.

OUTCOMES
--------
1. Background-removed VTEC / inter-product residuals already present in the
   daily master as columns ending:

       _background_residual

2. Validation daily bias and RMSE columns for comparisons against IGS and
   Madrigal.

GEOMAGNETIC PREDICTORS
----------------------
    Kp_max
    Ap_daily_from_ap
    Dst intensity  = -Dst_min
    SYM-H intensity = -SYM-H_min

The intensity sign convention means larger positive values correspond to
stronger negative Dst/SYM-H storm excursions.

STATISTICS
----------
For every outcome x predictor x lag:

    n
    Pearson r and p
    Spearman rho and p
    BH-FDR q values
    equal-year WLS slope
    95% year-cluster bootstrap CI
    95% 7-day moving-block bootstrap CI
    bootstrap_robust flag

OUTPUT
------
LAGGED_BOOTSTRAP

    KOH2_lagged_geomagnetic_all_tests.csv
    KOH2_lagged_geomagnetic_robust_associations.csv
    KOH2_lagged_geomagnetic_best_lag_per_outcome_predictor.csv
    KOH2_lagged_geomagnetic_key_validation_summary.csv
    KOH2_lagged_bootstrap_report.txt

SETTINGS
--------
Default:
    lags = 0, 1, 2 calendar days
    block length = 7 calendar days
    bootstrap replicates = 1000
    random seed = 20260822

No empirical TEC bias correction is applied.
"""

from pathlib import Path
import argparse
import calendar
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
INDEX_CACHE = Path("_INDEX_CACHE")
OUTPUT_ROOT = Path("LAGGED_BOOTSTRAP")

YEARS = list(
    range(
        2019,
        2027,
    )
)

LAGS = [
    0,
    1,
    2,
]

BLOCK_LENGTH_DAYS = 7
N_BOOTSTRAP = 1000
RANDOM_SEED = 20260822

MIN_N_CORRELATION = 8
MIN_N_BOOTSTRAP = 8
CI_LOW = 2.5
CI_HIGH = 97.5


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Run lagged geomagnetic robustness analysis with year-cluster "
                     "and calendar moving-block bootstrap.")
    )
    parser.add_argument("--master-file", required=True, type=Path,
                        help="KOH2 solar/geomagnetic daily master CSV.")
    parser.add_argument("--index-cache", type=Path, default=None,
                        help="Existing index cache. Default: _INDEX_CACHE beside the master file.")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Output directory. Default: LAGGED_BOOTSTRAP beside the master file.")
    parser.add_argument("--bootstrap-replicates", type=int, default=N_BOOTSTRAP,
                        help=f"Bootstrap replicates. Default: {N_BOOTSTRAP}.")
    parser.add_argument("--random-seed", type=int, default=RANDOM_SEED,
                        help=f"Random seed. Default: {RANDOM_SEED}.")
    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global MASTER_FILE, INDEX_CACHE, OUTPUT_ROOT, N_BOOTSTRAP, RANDOM_SEED
    MASTER_FILE = args.master_file.expanduser().resolve()
    INDEX_CACHE = (args.index_cache.expanduser().resolve()
                   if args.index_cache is not None
                   else MASTER_FILE.parent / "_INDEX_CACHE")
    OUTPUT_ROOT = (args.output_dir.expanduser().resolve()
                   if args.output_dir is not None
                   else MASTER_FILE.parent / "LAGGED_BOOTSTRAP")
    if args.bootstrap_replicates < 1:
        raise ValueError("--bootstrap-replicates must be >= 1")
    N_BOOTSTRAP = args.bootstrap_replicates
    RANDOM_SEED = args.random_seed


# =============================================================================
# OMNI CACHE FORMAT
# =============================================================================

OMNI2_COL_YEAR = 0
OMNI2_COL_DOY = 1
OMNI2_COL_HOUR = 2
OMNI2_COL_KP = 38
OMNI2_COL_DST = 40
OMNI2_COL_AP = 49

PREDICTORS = {
    "kp_max":
        "geomag_kp_max",
    "Ap_daily":
        "geomag_Ap_daily_from_ap",
    "dst_intensity":
        "geomag_dst_intensity_nt",
    "symh_intensity":
        "geomag_symh_intensity_nt",
}


# =============================================================================
# GENERAL HELPERS
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

    good = (
        p.notna()
        & np.isfinite(
            p
        )
    )

    if not good.any():
        return out

    values = p[
        good
    ].to_numpy(
        dtype=float
    )

    order = np.argsort(
        values
    )

    ranked = values[
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

    q_values = np.empty(
        m,
        dtype=float,
    )

    q_values[
        order
    ] = q_ranked

    out.loc[
        good,
        q_col,
    ] = q_values

    return out


def percentile_ci(
    values,
):
    x = np.asarray(
        values,
        dtype=float,
    )

    x = x[
        np.isfinite(
            x
        )
    ]

    if len(
        x
    ) == 0:
        return (
            np.nan,
            np.nan,
            0,
        )

    return (
        float(
            np.percentile(
                x,
                CI_LOW,
            )
        ),
        float(
            np.percentile(
                x,
                CI_HIGH,
            )
        ),
        int(
            len(
                x
            )
        ),
    )


def ci_excludes_zero(
    low,
    high,
):
    if (
        not np.isfinite(
            low
        )
        or not np.isfinite(
            high
        )
    ):
        return False

    return (
        low
        > 0
        or high
        < 0
    )


def ci_sign(
    low,
    high,
):
    if (
        not np.isfinite(
            low
        )
        or not np.isfinite(
            high
        )
    ):
        return 0

    if low > 0:
        return 1

    if high < 0:
        return -1

    return 0


# =============================================================================
# BUILD COMPLETE DAILY GEOMAGNETIC CALENDAR FROM EXISTING CACHE
# =============================================================================

def parse_omni2_cached_year(
    year: int,
):
    path = (
        INDEX_CACHE
        / "OMNI2"
        / f"omni2_{year}.dat"
    )

    if not path.is_file():
        raise FileNotFoundError(
            "\nMissing cached OMNI2 file:\n"
            f"    {path}\n\n"
            "Run analyze_koh2_solar_cycle25_2019_2026.py first."
        )

    arr = np.loadtxt(
        path,
        usecols=[
            OMNI2_COL_YEAR,
            OMNI2_COL_DOY,
            OMNI2_COL_HOUR,
            OMNI2_COL_KP,
            OMNI2_COL_DST,
            OMNI2_COL_AP,
        ],
        dtype=float,
    )

    if arr.ndim == 1:
        arr = arr.reshape(
            1,
            -1,
        )

    df = pd.DataFrame({
        "year":
            arr[
                :,
                0
            ].astype(
                int
            ),
        "doy":
            arr[
                :,
                1
            ].astype(
                int
            ),
        "hour":
            arr[
                :,
                2
            ].astype(
                int
            ),
        "kp_raw":
            arr[
                :,
                3
            ],
        "dst_nt":
            arr[
                :,
                4
            ],
        "ap":
            arr[
                :,
                5
            ],
    })

    df.loc[
        df[
            "kp_raw"
        ]
        >= 99,
        "kp_raw",
    ] = np.nan

    df.loc[
        np.abs(
            df[
                "dst_nt"
            ]
        )
        >= 99999,
        "dst_nt",
    ] = np.nan

    df.loc[
        df[
            "ap"
        ]
        >= 999,
        "ap",
    ] = np.nan

    df[
        "kp"
    ] = (
        df[
            "kp_raw"
        ]
        / 10.0
    )

    jan1 = pd.Timestamp(
        year=year,
        month=1,
        day=1,
    )

    df[
        "date"
    ] = (
        jan1
        + pd.to_timedelta(
            df[
                "doy"
            ]
            - 1,
            unit="D",
        )
    ).dt.normalize()

    rows = []

    for date, g in df.groupby(
        "date"
    ):
        rows.append({
            "date":
                date,
            "geomag_kp_max":
                pd.to_numeric(
                    g[
                        "kp"
                    ],
                    errors="coerce",
                ).max(),
            "geomag_Ap_daily_from_ap":
                pd.to_numeric(
                    g[
                        "ap"
                    ],
                    errors="coerce",
                ).mean(),
            "geomag_dst_min_nt":
                pd.to_numeric(
                    g[
                        "dst_nt"
                    ],
                    errors="coerce",
                ).min(),
        })

    return pd.DataFrame(
        rows
    )


def parse_symh_cached_year(
    year: int,
):
    path = (
        INDEX_CACHE
        / "OMNI_HRO_SYMH"
        / f"symh_5min_{year}.csv"
    )

    if not path.is_file():
        raise FileNotFoundError(
            "\nMissing cached reduced SYM-H file:\n"
            f"    {path}\n\n"
            "Run analyze_koh2_solar_cycle25_2019_2026.py first."
        )

    df = pd.read_csv(
        path
    )

    required = {
        "date",
        "symh_nt",
    }

    missing = (
        required
        - set(
            df.columns
        )
    )

    if missing:
        raise RuntimeError(
            f"{path} missing columns: {sorted(missing)}"
        )

    df[
        "date"
    ] = pd.to_datetime(
        df[
            "date"
        ],
        errors="coerce",
    ).dt.normalize()

    df[
        "symh_nt"
    ] = pd.to_numeric(
        df[
            "symh_nt"
        ],
        errors="coerce",
    )

    return (
        df.groupby(
            "date",
            as_index=False,
        )[
            "symh_nt"
        ]
        .min()
        .rename(
            columns={
                "symh_nt":
                    "geomag_symh_min_nt",
            }
        )
    )


def build_geomagnetic_calendar():
    pieces = []

    for year in YEARS:
        omni = parse_omni2_cached_year(
            year
        )

        symh = parse_symh_cached_year(
            year
        )

        g = omni.merge(
            symh,
            on="date",
            how="outer",
        )

        pieces.append(
            g
        )

    full = pd.concat(
        pieces,
        ignore_index=True,
    )

    full[
        "date"
    ] = pd.to_datetime(
        full[
            "date"
        ],
        errors="coerce",
    ).dt.normalize()

    full = (
        full.dropna(
            subset=[
                "date",
            ]
        )
        .drop_duplicates(
            subset=[
                "date",
            ],
            keep="last",
        )
        .sort_values(
            "date"
        )
        .reset_index(
            drop=True
        )
    )

    full[
        "geomag_dst_intensity_nt"
    ] = -pd.to_numeric(
        full[
            "geomag_dst_min_nt"
        ],
        errors="coerce",
    )

    full[
        "geomag_symh_intensity_nt"
    ] = -pd.to_numeric(
        full[
            "geomag_symh_min_nt"
        ],
        errors="coerce",
    )

    return full


# =============================================================================
# MASTER + TRUE CALENDAR LAGS
# =============================================================================

def load_master_with_lags():
    if not MASTER_FILE.is_file():
        raise FileNotFoundError(
            f"Missing daily master file:\n{MASTER_FILE}"
        )

    master = pd.read_csv(
        MASTER_FILE
    )

    if "date" not in master.columns:
        raise RuntimeError(
            "Daily master does not contain date."
        )

    master[
        "date"
    ] = pd.to_datetime(
        master[
            "date"
        ],
        errors="coerce",
    ).dt.normalize()

    master = master.dropna(
        subset=[
            "date",
        ]
    ).copy()

    geomag = build_geomagnetic_calendar()

    for lag in LAGS:
        shifted = geomag.copy()

        # Activity on date A is assigned to response date A + lag.
        shifted[
            "date"
        ] = (
            shifted[
                "date"
            ]
            + pd.to_timedelta(
                lag,
                unit="D",
            )
        )

        rename = {}

        for original in PREDICTORS.values():
            rename[
                original
            ] = (
                f"{original}_lag{lag}"
            )

        keep = [
            "date",
        ] + list(
            PREDICTORS.values()
        )

        shifted = shifted[
            keep
        ].rename(
            columns=rename
        )

        master = master.merge(
            shifted,
            on="date",
            how="left",
        )

    return master


# =============================================================================
# SELECT OUTCOMES
# =============================================================================

def select_outcomes(
    master: pd.DataFrame,
):
    background = []

    validation = []

    for col in master.columns:
        lower = col.lower()

        if lower.endswith(
            "_background_residual"
        ):
            # Focus on observed/reference VTEC and inter-product differences.
            if (
                "pytecgg"
                in lower
                or "igs"
                in lower
                or "madrigal"
                in lower
            ):
                background.append(
                    col
                )

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
            validation.append(
                col
            )

    background = sorted(
        set(
            background
        )
    )

    validation = sorted(
        set(
            validation
        )
    )

    return (
        background,
        validation,
    )


# =============================================================================
# OBSERVED CORRELATIONS AND EQUAL-YEAR SLOPE
# =============================================================================

def valid_pair(
    df: pd.DataFrame,
    y_col: str,
    x_col: str,
):
    g = df[
        [
            "date",
            "year",
            y_col,
            x_col,
        ]
    ].copy()

    g[
        y_col
    ] = pd.to_numeric(
        g[
            y_col
        ],
        errors="coerce",
    )

    g[
        x_col
    ] = pd.to_numeric(
        g[
            x_col
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
        subset=[
            "date",
            "year",
            y_col,
            x_col,
        ]
    )

    return g


def correlations(
    g: pd.DataFrame,
    y_col: str,
    x_col: str,
):
    n = len(
        g
    )

    result = {
        "n":
            n,
        "pearson_r":
            np.nan,
        "pearson_p":
            np.nan,
        "spearman_rho":
            np.nan,
        "spearman_p":
            np.nan,
    }

    if n < MIN_N_CORRELATION:
        return result

    x = g[
        x_col
    ].to_numpy(
        dtype=float
    )

    y = g[
        y_col
    ].to_numpy(
        dtype=float
    )

    if (
        np.std(
            x
        )
        == 0
        or np.std(
            y
        )
        == 0
    ):
        return result

    p = scipy_stats.pearsonr(
        x,
        y,
    )

    s = scipy_stats.spearmanr(
        x,
        y,
    )

    result.update({
        "pearson_r":
            float(
                p.statistic
            ),
        "pearson_p":
            float(
                p.pvalue
            ),
        "spearman_rho":
            float(
                s.statistic
            ),
        "spearman_p":
            float(
                s.pvalue
            ),
    })

    return result


def equal_year_weighted_slope_arrays(
    x,
    y,
    cluster_labels,
):
    x = np.asarray(
        x,
        dtype=float,
    )

    y = np.asarray(
        y,
        dtype=float,
    )

    labels = np.asarray(
        cluster_labels
    )

    good = (
        np.isfinite(
            x
        )
        & np.isfinite(
            y
        )
    )

    x = x[
        good
    ]

    y = y[
        good
    ]

    labels = labels[
        good
    ]

    if len(
        x
    ) < MIN_N_BOOTSTRAP:
        return np.nan

    unique_labels, counts = np.unique(
        labels,
        return_counts=True,
    )

    if len(
        unique_labels
    ) < 2:
        return np.nan

    count_map = {
        label:
            count
        for label, count in zip(
            unique_labels,
            counts,
        )
    }

    w = np.array(
        [
            1.0
            / count_map[
                label
            ]
            for label in labels
        ],
        dtype=float,
    )

    sw = np.sum(
        w
    )

    xbar = np.sum(
        w
        * x
    ) / sw

    ybar = np.sum(
        w
        * y
    ) / sw

    denom = np.sum(
        w
        * (
            x
            - xbar
        )
        ** 2
    )

    if denom <= 0:
        return np.nan

    slope = (
        np.sum(
            w
            * (
                x
                - xbar
            )
            * (
                y
                - ybar
            )
        )
        / denom
    )

    return float(
        slope
    )


def observed_equal_year_slope(
    g: pd.DataFrame,
    y_col: str,
    x_col: str,
):
    return equal_year_weighted_slope_arrays(
        g[
            x_col
        ].to_numpy(
            dtype=float
        ),
        g[
            y_col
        ].to_numpy(
            dtype=float
        ),
        g[
            "year"
        ].to_numpy(),
    )


# =============================================================================
# YEAR-CLUSTER BOOTSTRAP
# =============================================================================

def year_cluster_bootstrap_slopes(
    g: pd.DataFrame,
    y_col: str,
    x_col: str,
    rng: np.random.Generator,
):
    years = sorted(
        int(
            y
        )
        for y in g[
            "year"
        ].unique()
    )

    if len(
        years
    ) < 2:
        return []

    by_year = {
        year:
            g[
                g[
                    "year"
                ]
                == year
            ]
        for year in years
    }

    slopes = []

    for _ in range(
        N_BOOTSTRAP
    ):
        selected = rng.choice(
            years,
            size=len(
                years
            ),
            replace=True,
        )

        x_parts = []

        y_parts = []

        labels = []

        for cluster_instance, year in enumerate(
            selected
        ):
            part = by_year[
                int(
                    year
                )
            ]

            if part.empty:
                continue

            x_part = part[
                x_col
            ].to_numpy(
                dtype=float
            )

            y_part = part[
                y_col
            ].to_numpy(
                dtype=float
            )

            x_parts.append(
                x_part
            )

            y_parts.append(
                y_part
            )

            labels.append(
                np.full(
                    len(
                        part
                    ),
                    cluster_instance,
                    dtype=int,
                )
            )

        if not x_parts:
            continue

        slope = equal_year_weighted_slope_arrays(
            np.concatenate(
                x_parts
            ),
            np.concatenate(
                y_parts
            ),
            np.concatenate(
                labels
            ),
        )

        if np.isfinite(
            slope
        ):
            slopes.append(
                slope
            )

    return slopes


# =============================================================================
# 7-DAY CALENDAR MOVING-BLOCK BOOTSTRAP
# =============================================================================

def calendar_metadata(
    master: pd.DataFrame,
):
    metadata = {}

    date_to_rows = {}

    for idx, date in master[
        "date"
    ].items():
        key = pd.Timestamp(
            date
        ).normalize()

        date_to_rows.setdefault(
            key,
            []
        ).append(
            idx
        )

    for year in YEARS:
        n_days = (
            366
            if calendar.isleap(
                year
            )
            else 365
        )

        jan1 = pd.Timestamp(
            year=year,
            month=1,
            day=1,
        )

        dates = pd.date_range(
            jan1,
            periods=n_days,
            freq="D",
        )

        metadata[
            year
        ] = {
            "dates":
                dates,
            "n_days":
                n_days,
        }

    return (
        metadata,
        date_to_rows,
    )


def generate_block_bootstrap_indices(
    master: pd.DataFrame,
    rng: np.random.Generator,
):
    metadata, date_to_rows = calendar_metadata(
        master
    )

    samples = []

    for _ in range(
        N_BOOTSTRAP
    ):
        sampled_rows = []

        for year in YEARS:
            dates = metadata[
                year
            ][
                "dates"
            ]

            n_days = metadata[
                year
            ][
                "n_days"
            ]

            n_blocks = int(
                math.ceil(
                    n_days
                    / BLOCK_LENGTH_DAYS
                )
            )

            starts = rng.integers(
                0,
                n_days,
                size=n_blocks,
            )

            sampled_offsets = []

            for start in starts:
                block = (
                    start
                    + np.arange(
                        BLOCK_LENGTH_DAYS
                    )
                ) % n_days

                sampled_offsets.extend(
                    block.tolist()
                )

            sampled_offsets = sampled_offsets[
                :n_days
            ]

            for offset in sampled_offsets:
                date = dates[
                    int(
                        offset
                    )
                ]

                rows = date_to_rows.get(
                    date,
                )

                if rows:
                    sampled_rows.extend(
                        rows
                    )

        samples.append(
            np.asarray(
                sampled_rows,
                dtype=int,
            )
        )

    return samples


def block_bootstrap_slopes(
    master: pd.DataFrame,
    y_col: str,
    x_col: str,
    block_samples,
):
    slopes = []

    for indices in block_samples:
        if len(
            indices
        ) == 0:
            continue

        sample = master.loc[
            indices,
            [
                "year",
                y_col,
                x_col,
            ]
        ].copy()

        sample[
            y_col
        ] = pd.to_numeric(
            sample[
                y_col
            ],
            errors="coerce",
        )

        sample[
            x_col
        ] = pd.to_numeric(
            sample[
                x_col
            ],
            errors="coerce",
        )

        sample = sample.dropna(
            subset=[
                "year",
                y_col,
                x_col,
            ]
        )

        if len(
            sample
        ) < MIN_N_BOOTSTRAP:
            continue

        slope = equal_year_weighted_slope_arrays(
            sample[
                x_col
            ].to_numpy(
                dtype=float
            ),
            sample[
                y_col
            ].to_numpy(
                dtype=float
            ),
            sample[
                "year"
            ].to_numpy(),
        )

        if np.isfinite(
            slope
        ):
            slopes.append(
                slope
            )

    return slopes


# =============================================================================
# ANALYSIS
# =============================================================================

def run_all_tests(
    master: pd.DataFrame,
    outcomes,
):
    rng_blocks = np.random.default_rng(
        RANDOM_SEED
    )

    print(
        f"Generating {N_BOOTSTRAP} seven-day calendar block-bootstrap samples ..."
    )

    block_samples = generate_block_bootstrap_indices(
        master,
        rng_blocks,
    )

    rows = []

    total_tests = (
        len(
            outcomes
        )
        * len(
            PREDICTORS
        )
        * len(
            LAGS
        )
    )

    test_number = 0

    for outcome_type, y_col in outcomes:
        for predictor_label, predictor_base in PREDICTORS.items():
            for lag in LAGS:
                test_number += 1

                x_col = (
                    f"{predictor_base}_lag{lag}"
                )

                print(
                    f"{test_number:3d}/{total_tests:3d}  "
                    f"{outcome_type} | {y_col} | "
                    f"{predictor_label} lag+{lag}"
                )

                g = valid_pair(
                    master,
                    y_col,
                    x_col,
                )

                corr = correlations(
                    g,
                    y_col,
                    x_col,
                )

                slope = observed_equal_year_slope(
                    g,
                    y_col,
                    x_col,
                )

                # Independent deterministic RNG stream per test for the
                # year-cluster bootstrap.
                seed = (
                    RANDOM_SEED
                    + test_number
                    * 1009
                )

                rng_cluster = np.random.default_rng(
                    seed
                )

                cluster_slopes = year_cluster_bootstrap_slopes(
                    g,
                    y_col,
                    x_col,
                    rng_cluster,
                )

                block_slopes = block_bootstrap_slopes(
                    master,
                    y_col,
                    x_col,
                    block_samples,
                )

                (
                    cluster_low,
                    cluster_high,
                    cluster_n,
                ) = percentile_ci(
                    cluster_slopes
                )

                (
                    block_low,
                    block_high,
                    block_n,
                ) = percentile_ci(
                    block_slopes
                )

                cluster_sign = ci_sign(
                    cluster_low,
                    cluster_high,
                )

                block_sign = ci_sign(
                    block_low,
                    block_high,
                )

                robust = (
                    cluster_sign
                    != 0
                    and block_sign
                    != 0
                    and cluster_sign
                    == block_sign
                )

                rows.append({
                    "outcome_type":
                        outcome_type,
                    "outcome":
                        y_col,
                    "predictor":
                        predictor_label,
                    "lag_days_after_activity":
                        lag,
                    "lag_definition":
                        (
                            "response on D vs activity on D"
                            if lag == 0
                            else (
                                f"response on D vs activity on D-{lag}"
                            )
                        ),
                    **corr,
                    "equal_year_slope":
                        slope,
                    "year_cluster_ci_low":
                        cluster_low,
                    "year_cluster_ci_high":
                        cluster_high,
                    "year_cluster_boot_valid":
                        cluster_n,
                    "block7_ci_low":
                        block_low,
                    "block7_ci_high":
                        block_high,
                    "block7_boot_valid":
                        block_n,
                    "bootstrap_robust":
                        robust,
                    "bootstrap_sign":
                        (
                            cluster_sign
                            if robust
                            else 0
                        ),
                })

    out = pd.DataFrame(
        rows
    )

    out = add_bh_fdr(
        out,
        "pearson_p",
        "pearson_q_bh",
    )

    out = add_bh_fdr(
        out,
        "spearman_p",
        "spearman_q_bh",
    )

    return out


# =============================================================================
# BEST LAG / KEY SUMMARY
# =============================================================================

def best_lag_table(
    tests: pd.DataFrame,
):
    rows = []

    for (
        outcome,
        predictor,
    ), g in tests.groupby(
        [
            "outcome",
            "predictor",
        ]
    ):
        h = g.copy()

        h[
            "abs_spearman"
        ] = np.abs(
            pd.to_numeric(
                h[
                    "spearman_rho"
                ],
                errors="coerce",
            )
        )

        h[
            "abs_slope"
        ] = np.abs(
            pd.to_numeric(
                h[
                    "equal_year_slope"
                ],
                errors="coerce",
            )
        )

        h = h.sort_values(
            [
                "bootstrap_robust",
                "abs_spearman",
                "abs_slope",
            ],
            ascending=[
                False,
                False,
                False,
            ],
        )

        if h.empty:
            continue

        row = h.iloc[
            0
        ].drop(
            labels=[
                "abs_spearman",
                "abs_slope",
            ],
            errors="ignore",
        ).to_dict()

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    )


def key_validation_summary(
    tests: pd.DataFrame,
):
    patterns = [
        "pytecgg_veq_vs_igs_rmse_tecu",
        "pytecgg_vtec_vs_igs_rmse_tecu",
        "pyoasis_vtec_vs_igs_rmse_tecu",
        "pytecgg_veq_vs_madrigal_rmse_tecu",
        "pytecgg_vtec_vs_madrigal_rmse_tecu",
        "pyoasis_vtec_vs_madrigal_rmse_tecu",
        "igs_minus_madrigal_mean_tecu_background_residual",
    ]

    return tests[
        tests[
            "outcome"
        ].isin(
            patterns
        )
    ].sort_values(
        [
            "outcome",
            "predictor",
            "lag_days_after_activity",
        ]
    )


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
        320,
    )

    print(
        "=" * 126
    )

    print(
        "KOH2 LAGGED GEOMAGNETIC ROBUSTNESS ANALYSIS"
    )

    print(
        "=" * 126
    )

    print(
        "Daily master:",
        MASTER_FILE,
    )

    print(
        "Lags:",
        LAGS,
        "calendar days",
    )

    print(
        "Lag +1 means response today vs geomagnetic activity yesterday."
    )

    print(
        "Bootstrap replicates:",
        N_BOOTSTRAP,
    )

    print(
        "Moving block length:",
        BLOCK_LENGTH_DAYS,
        "calendar days",
    )

    master = load_master_with_lags()

    background, validation = select_outcomes(
        master
    )

    print()
    print(
        "Background-residual outcomes:",
        len(
            background
        ),
    )

    print(
        "Validation bias/RMSE outcomes:",
        len(
            validation
        ),
    )

    if not background:
        raise RuntimeError(
            "No *_background_residual columns were found in the daily master."
        )

    if not validation:
        print(
            "[WARNING] No validation bias/RMSE columns found; "
            "background residual analysis will still run."
        )

    outcomes = (
        [
            (
                "background_residual",
                col,
            )
            for col in background
        ]
        + [
            (
                "validation_error",
                col,
            )
            for col in validation
        ]
    )

    tests = run_all_tests(
        master,
        outcomes,
    )

    robust = tests[
        tests[
            "bootstrap_robust"
        ]
        == True
    ].copy()

    if not robust.empty:
        robust[
            "abs_spearman"
        ] = np.abs(
            pd.to_numeric(
                robust[
                    "spearman_rho"
                ],
                errors="coerce",
            )
        )

        robust = robust.sort_values(
            [
                "abs_spearman",
                "n",
            ],
            ascending=[
                False,
                False,
            ],
        ).drop(
            columns=[
                "abs_spearman",
            ]
        )

    best = best_lag_table(
        tests
    )

    key = key_validation_summary(
        tests
    )

    all_file = (
        OUTPUT_ROOT
        / "KOH2_lagged_geomagnetic_all_tests.csv"
    )

    robust_file = (
        OUTPUT_ROOT
        / "KOH2_lagged_geomagnetic_robust_associations.csv"
    )

    best_file = (
        OUTPUT_ROOT
        / "KOH2_lagged_geomagnetic_best_lag_per_outcome_predictor.csv"
    )

    key_file = (
        OUTPUT_ROOT
        / "KOH2_lagged_geomagnetic_key_validation_summary.csv"
    )

    report_file = (
        OUTPUT_ROOT
        / "KOH2_lagged_bootstrap_report.txt"
    )

    tests.to_csv(
        all_file,
        index=False,
    )

    robust.to_csv(
        robust_file,
        index=False,
    )

    best.to_csv(
        best_file,
        index=False,
    )

    key.to_csv(
        key_file,
        index=False,
    )

    with open(
        report_file,
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            "KOH2 LAGGED GEOMAGNETIC ROBUSTNESS ANALYSIS\n"
        )

        f.write(
            "=" * 110
            + "\n\n"
        )

        f.write(
            "Lag convention:\n"
            "  lag 0 = response D vs activity D\n"
            "  lag 1 = response D vs activity D-1\n"
            "  lag 2 = response D vs activity D-2\n\n"
        )

        f.write(
            f"Bootstrap replicates: {N_BOOTSTRAP}\n"
        )

        f.write(
            f"Moving block length: {BLOCK_LENGTH_DAYS} calendar days\n\n"
        )

        f.write(
            "A result is bootstrap_robust only when BOTH the year-cluster "
            "95% CI and the 7-day block-bootstrap 95% CI exclude zero with "
            "the same sign.\n\n"
        )

        f.write(
            "ROBUST ASSOCIATIONS\n"
        )

        f.write(
            "-" * 110
            + "\n"
        )

        if robust.empty:
            f.write(
                "None.\n"
            )
        else:
            f.write(
                robust.to_string(
                    index=False
                )
            )

        f.write(
            "\n\nKEY VALIDATION LAG SUMMARY\n"
        )

        f.write(
            "-" * 110
            + "\n"
        )

        f.write(
            key.to_string(
                index=False
            )
        )

    print()
    print(
        "=" * 126
    )

    print(
        "LAGGED BOOTSTRAP ANALYSIS COMPLETE"
    )

    print(
        "=" * 126
    )

    print(
        "All tests:",
        all_file,
    )

    print(
        "Robust associations:",
        robust_file,
    )

    print(
        "Best lag table:",
        best_file,
    )

    print(
        "Key validation summary:",
        key_file,
    )

    print(
        "Report:",
        report_file,
    )

    print()
    print(
        "BOOTSTRAP-ROBUST ASSOCIATIONS"
    )

    print(
        "-" * 126
    )

    if robust.empty:
        print(
            "None."
        )
    else:
        display_cols = [
            "outcome_type",
            "outcome",
            "predictor",
            "lag_days_after_activity",
            "n",
            "spearman_rho",
            "spearman_q_bh",
            "equal_year_slope",
            "year_cluster_ci_low",
            "year_cluster_ci_high",
            "block7_ci_low",
            "block7_ci_high",
        ]

        print(
            robust[
                display_cols
            ].head(
                30
            ).to_string(
                index=False
            )
        )

    print()
    print(
        "KEY VALIDATION LAG SUMMARY"
    )

    print(
        "-" * 126
    )

    if key.empty:
        print(
            "No key validation outcomes found."
        )
    else:
        display_cols = [
            "outcome",
            "predictor",
            "lag_days_after_activity",
            "n",
            "pearson_r",
            "spearman_rho",
            "spearman_q_bh",
            "equal_year_slope",
            "year_cluster_ci_low",
            "year_cluster_ci_high",
            "block7_ci_low",
            "block7_ci_high",
            "bootstrap_robust",
        ]

        print(
            key[
                display_cols
            ].to_string(
                index=False
            )
        )


if __name__ == "__main__":
    main()
