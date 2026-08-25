from __future__ import annotations

r"""
KOH2 SOLAR CYCLE 25 / GEOMAGNETIC ACTIVITY ANALYSIS, 2019-2026
==============================================================

PURPOSE
-------
Quantitatively test two different questions:

A. BACKGROUND SOLAR-CYCLE DEPENDENCE
   Does observed/reference VTEC increase with solar activity after controlling
   for season?

       VTEC ~ F10.7 + seasonal harmonics

   Daily Sunspot Number (SSN) is also reported as an independent solar-activity
   proxy, but F10.7 and SSN are NOT entered together in the same regression
   because they are strongly collinear.

B. GEOMAGNETIC-DISTURBANCE DEPENDENCE
   After removing the fitted F10.7 + seasonal background, do remaining VTEC
   anomalies or validation residuals depend on geomagnetic disturbance?

       background residual ~ Kp_max
       background residual ~ Ap_daily
       background residual ~ -Dst_min
       background residual ~ -SYM-H_min

This is deliberately different from simply correlating yearly means.

DATA ALREADY PRODUCED
---------------------
Main harmonized station-level series:
    
        KOH2_2019_2026_common_hour_values.csv

Additional method-validation daily statistics:
    
        KOH2_2019_2026_IGS_validation_daily_statistics.csv

    
        KOH2_2019_2026_Madrigal_daily_statistics.csv

SOLAR / GEOMAGNETIC SOURCES
---------------------------
1. NASA/SPDF OMNI2 hourly:
       Kp, ap, Dst, F10.7
   Yearly files:
       https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/
           omni2_YYYY.dat

   OMNI2 word numbers used:
       39 Kp, stored in tenths (e.g. 33 = 3+)
       41 Dst, nT
       50 ap index
       51 F10.7, sfu

   From the 3-hour ap sequence, a daily Ap-like value is calculated as the
   daily arithmetic mean of the valid 3-hour ap values.  Output label:
       Ap_daily_from_ap

2. NASA/SPDF High Resolution OMNI, 5-minute:
       SYM-H
   Yearly files:
       https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/
           omni_5minYYYY.asc

   HRO word number:
       42 SYM/H, nT

   5-minute data are used instead of 1-minute data to reduce downloads while
   retaining storm-time minima adequately for this study.

3. WDC-SILSO:
       Daily total Sunspot Number, Version 2.0
       SN_d_tot_V2.0.csv

PRIMARY DAILY TEC VARIABLES
---------------------------
Derived from the common-hour table:
    pytecgg_veq_mean_tecu
    igs_vtec_mean_tecu
    madrigal_vtec_mean_tecu
    pyiri_vtec_mean_tecu

Common-hour reference/method differences:
    pytecgg_minus_igs_tecu
    pytecgg_minus_madrigal_tecu
    igs_minus_madrigal_tecu
    pytecgg_minus_pyiri_tecu
    igs_minus_pyiri_tecu
    madrigal_minus_pyiri_tecu

Additional validation-performance variables (when available):
    PyTECGg vs IGS daily bias / RMSE
    pyOASIS vs IGS daily bias / RMSE
    PyTECGg vs Madrigal daily bias / RMSE
    pyOASIS vs Madrigal daily bias / RMSE

SCIENTIFIC CONTROLS
-------------------
Seasonality is represented by annual and semiannual Fourier terms:
    sin(2*pi*DOY/365.25)
    cos(2*pi*DOY/365.25)
    sin(4*pi*DOY/365.25)
    cos(4*pi*DOY/365.25)

Background regression:
    y = b0 + b1*F10.7 + seasonal terms

Quiet-day subset:
    Kp_max < 4
    Dst_min > -30 nT
    SYM-H_min > -30 nT

Storm-day flag:
    Kp_max >= 5 OR Dst_min <= -50 nT OR SYM-H_min <= -50 nT

Because the available KOH2 days are irregularly distributed among years,
all annual summaries are equal-day summaries.

OUTPUT ROOT
-----------
    <output-dir>

Main outputs:
    KOH2_2019_2026_solar_geomagnetic_daily_master.csv
    KOH2_2019_2026_yearly_equal_day_summary.csv
    KOH2_2019_2026_correlation_all_days.csv
    KOH2_2019_2026_correlation_quiet_days.csv
    KOH2_2019_2026_background_regression_F107_season.csv
    KOH2_2019_2026_background_residual_geomagnetic_correlations.csv
    KOH2_2019_2026_validation_error_geomagnetic_correlations.csv
    KOH2_2019_2026_storm_vs_quiet_summary.csv
    KOH2_2019_2026_analysis_report.txt

No empirical TEC bias correction is applied.
"""

from pathlib import Path
import argparse
from datetime import datetime, timezone
import math
import os
import re
import time

import numpy as np
import pandas as pd
import requests

try:
    from scipy import stats as scipy_stats
except ImportError as exc:
    raise RuntimeError(
        "\nSciPy is required for statistical significance calculations.\n"
        "Activate pytecgg_env and run:\n\n"
        "    python -m pip install scipy\n"
    ) from exc


# =============================================================================
# SETTINGS
# =============================================================================

YEARS = list(range(2019, 2027))
STATION = "KOH2"

# Runtime paths are configured from CLI arguments in configure_runtime().
COMMON_HOUR_FILE = Path(f"{STATION}_2019_2026_common_hour_values.csv")
IGS_DAILY_FILE: Path | None = None
MADRIGAL_DAILY_FILE: Path | None = None
OUTPUT_ROOT = Path("TEC_SOLAR_CYCLE25_ANALYSIS_2019_2026")
INDEX_CACHE = OUTPUT_ROOT / "_INDEX_CACHE"

OMNI2_BASE = (
    "https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni"
)

HRO_BASE = (
    "https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni"
)

SILSO_DAILY_URL = (
    "https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv"
)

REQUEST_TIMEOUT = (
    30,
    120,
)

# Delete the large ~33 MB annual HRO file after extracting just SYM-H.
# The tiny extracted annual CSV is kept, so reruns require no re-download.
DELETE_RAW_HRO_AFTER_EXTRACT = True

QUIET_KP_MAX = 4.0
QUIET_DST_MIN_NT = -30.0
QUIET_SYMH_MIN_NT = -30.0

STORM_KP_MAX = 5.0
STORM_DST_MIN_NT = -50.0
STORM_SYMH_MIN_NT = -50.0

MIN_CORRELATION_N = 4
MIN_REGRESSION_N = 8


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the KOH2 Solar Cycle 25 / geomagnetic activity analysis "
            "for the validated 2019-2026 common-hour products."
        )
    )
    parser.add_argument("--common-hour-file", required=True, type=Path,
                        help="Validated KOH2 common-hour values CSV.")
    parser.add_argument("--igs-daily-file", type=Path, default=None,
                        help="Optional IGS validation daily-statistics CSV.")
    parser.add_argument("--madrigal-daily-file", type=Path, default=None,
                        help="Optional Madrigal validation daily-statistics CSV.")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Directory for analysis products.")
    parser.add_argument(
        "--index-cache", type=Path, default=None,
        help=("OMNI/SYM-H/SILSO cache directory. Default: <output-dir>/_INDEX_CACHE. "
              "For equivalence testing, point this to the existing operational cache."),
    )
    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global COMMON_HOUR_FILE, IGS_DAILY_FILE, MADRIGAL_DAILY_FILE
    global OUTPUT_ROOT, INDEX_CACHE
    COMMON_HOUR_FILE = args.common_hour_file.expanduser().resolve()
    IGS_DAILY_FILE = (args.igs_daily_file.expanduser().resolve()
                      if args.igs_daily_file is not None else None)
    MADRIGAL_DAILY_FILE = (args.madrigal_daily_file.expanduser().resolve()
                           if args.madrigal_daily_file is not None else None)
    OUTPUT_ROOT = args.output_dir.expanduser().resolve()
    INDEX_CACHE = (args.index_cache.expanduser().resolve()
                   if args.index_cache is not None else OUTPUT_ROOT / "_INDEX_CACHE")


# =============================================================================
# DOWNLOAD HELPERS
# =============================================================================

def download_file(
    url: str,
    path: Path,
    min_bytes: int = 100,
):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if (
        path.is_file()
        and path.stat().st_size >= min_bytes
    ):
        return path

    tmp = path.with_suffix(
        path.suffix
        + ".part"
    )

    if tmp.exists():
        tmp.unlink()

    print(
        "Downloading:",
        url,
    )

    response = requests.get(
        url,
        stream=True,
        timeout=REQUEST_TIMEOUT,
    )

    response.raise_for_status()

    size = 0

    with open(
        tmp,
        "wb",
    ) as f:
        for chunk in response.iter_content(
            chunk_size=1024 * 1024
        ):
            if not chunk:
                continue

            f.write(
                chunk
            )

            size += len(
                chunk
            )

    if size < min_bytes:
        if tmp.exists():
            tmp.unlink()

        raise RuntimeError(
            f"Downloaded file too small ({size} bytes): {url}"
        )

    tmp.replace(
        path
    )

    print(
        f"  saved {size / (1024 ** 2):.1f} MB -> {path}"
    )

    return path


# =============================================================================
# OMNI2 HOURLY: Kp, ap, Dst, F10.7
# =============================================================================

OMNI2_COL_YEAR = 0
OMNI2_COL_DOY = 1
OMNI2_COL_HOUR = 2
OMNI2_COL_KP = 38
OMNI2_COL_DST = 40
OMNI2_COL_AP = 49
OMNI2_COL_F107 = 50


def parse_omni2_year(
    year: int,
):
    cache = (
        INDEX_CACHE
        / "OMNI2"
    )

    path = (
        cache
        / f"omni2_{year}.dat"
    )

    url = (
        f"{OMNI2_BASE}/omni2_{year}.dat"
    )

    download_file(
        url,
        path,
        min_bytes=100_000,
    )

    usecols = [
        OMNI2_COL_YEAR,
        OMNI2_COL_DOY,
        OMNI2_COL_HOUR,
        OMNI2_COL_KP,
        OMNI2_COL_DST,
        OMNI2_COL_AP,
        OMNI2_COL_F107,
    ]

    arr = np.loadtxt(
        path,
        usecols=usecols,
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
        "f107_sfu":
            arr[
                :,
                6
            ],
    })

    # OMNI2 fill values documented in omni2.text.
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

    df.loc[
        df[
            "f107_sfu"
        ]
        >= 999.9,
        "f107_sfu",
    ] = np.nan

    # OMNI Kp is stored in tenths:
    # 33 = 3+, 57 = 6-, 40 = 4, ...
    df[
        "kp"
    ] = (
        df[
            "kp_raw"
        ]
        / 10.0
    )

    # Build UTC timestamp.
    jan1 = pd.Timestamp(
        year=year,
        month=1,
        day=1,
        tz="UTC",
    )

    df[
        "epoch"
    ] = (
        jan1
        + pd.to_timedelta(
            df[
                "doy"
            ]
            - 1,
            unit="D",
        )
        + pd.to_timedelta(
            df[
                "hour"
            ],
            unit="h",
        )
    )

    return df


def daily_omni2_indices(
    years,
):
    rows = []

    for year in years:
        print()
        print(
            f"OMNI2 hourly indices: {year}"
        )

        df = parse_omni2_year(
            year
        )

        df[
            "date"
        ] = df[
            "epoch"
        ].dt.date.astype(
            str
        )

        for date, g in df.groupby(
            "date"
        ):
            kp = pd.to_numeric(
                g[
                    "kp"
                ],
                errors="coerce",
            )

            ap = pd.to_numeric(
                g[
                    "ap"
                ],
                errors="coerce",
            )

            dst = pd.to_numeric(
                g[
                    "dst_nt"
                ],
                errors="coerce",
            )

            f107 = pd.to_numeric(
                g[
                    "f107_sfu"
                ],
                errors="coerce",
            )

            rows.append({
                "date":
                    date,
                "omni_kp_max":
                    kp.max(
                        skipna=True
                    ),
                "omni_kp_mean":
                    kp.mean(
                        skipna=True
                    ),
                # Hourly records repeat each 3-hour ap value.
                # Mean over the day therefore reproduces the arithmetic
                # mean of the eight 3-hour ap values (daily Ap-like value).
                "omni_Ap_daily_from_ap":
                    ap.mean(
                        skipna=True
                    ),
                "omni_ap_max":
                    ap.max(
                        skipna=True
                    ),
                "omni_dst_min_nt":
                    dst.min(
                        skipna=True
                    ),
                "omni_dst_mean_nt":
                    dst.mean(
                        skipna=True
                    ),
                "omni_f107_sfu":
                    f107.median(
                        skipna=True
                    ),
                "n_omni_kp_hours":
                    int(
                        kp.notna().sum()
                    ),
                "n_omni_dst_hours":
                    int(
                        dst.notna().sum()
                    ),
            })

    return pd.DataFrame(
        rows
    )


# =============================================================================
# HIGH RES OMNI 5-MIN: SYM-H
# =============================================================================

# HRO word 42 => zero-based column 41.
HRO_COL_YEAR = 0
HRO_COL_DOY = 1
HRO_COL_HOUR = 2
HRO_COL_MINUTE = 3
HRO_COL_SYMH = 41


def extracted_symh_file(
    year: int,
):
    return (
        INDEX_CACHE
        / "OMNI_HRO_SYMH"
        / f"symh_5min_{year}.csv"
    )


def extract_symh_year(
    year: int,
):
    reduced = extracted_symh_file(
        year
    )

    if (
        reduced.is_file()
        and reduced.stat().st_size > 1000
    ):
        return pd.read_csv(
            reduced
        )

    raw_dir = (
        INDEX_CACHE
        / "OMNI_HRO_RAW"
    )

    raw = (
        raw_dir
        / f"omni_5min{year}.asc"
    )

    url = (
        f"{HRO_BASE}/omni_5min{year}.asc"
    )

    download_file(
        url,
        raw,
        min_bytes=1_000_000,
    )

    print(
        f"Extracting 5-minute SYM-H locally: {year}"
    )

    # Read only the five required fields from the large annual ASCII file.
    arr = np.loadtxt(
        raw,
        usecols=[
            HRO_COL_YEAR,
            HRO_COL_DOY,
            HRO_COL_HOUR,
            HRO_COL_MINUTE,
            HRO_COL_SYMH,
        ],
        dtype=float,
    )

    if arr.ndim == 1:
        arr = arr.reshape(
            1,
            -1,
        )

    symh = pd.DataFrame({
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
        "minute":
            arr[
                :,
                3
            ].astype(
                int
            ),
        "symh_nt":
            arr[
                :,
                4
            ],
    })

    # HRO geomagnetic index fill is I6; reject extreme fill/sentinel values.
    symh.loc[
        np.abs(
            symh[
                "symh_nt"
            ]
        )
        >= 9999,
        "symh_nt",
    ] = np.nan

    jan1 = pd.Timestamp(
        year=year,
        month=1,
        day=1,
        tz="UTC",
    )

    symh[
        "epoch"
    ] = (
        jan1
        + pd.to_timedelta(
            symh[
                "doy"
            ]
            - 1,
            unit="D",
        )
        + pd.to_timedelta(
            symh[
                "hour"
            ],
            unit="h",
        )
        + pd.to_timedelta(
            symh[
                "minute"
            ],
            unit="m",
        )
    )

    symh[
        "date"
    ] = symh[
        "epoch"
    ].dt.date.astype(
        str
    )

    reduced.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    symh[
        [
            "date",
            "epoch",
            "symh_nt",
        ]
    ].to_csv(
        reduced,
        index=False,
    )

    if DELETE_RAW_HRO_AFTER_EXTRACT:
        try:
            raw.unlink()

            print(
                "  deleted large raw HRO file after extraction"
            )
        except Exception as exc:
            print(
                "  [WARNING] could not delete raw HRO:",
                repr(
                    exc
                ),
            )

    return symh[
        [
            "date",
            "epoch",
            "symh_nt",
        ]
    ]


def daily_symh_indices(
    years,
):
    rows = []

    for year in years:
        print()
        print(
            f"High-resolution OMNI SYM-H: {year}"
        )

        df = extract_symh_year(
            year
        )

        df[
            "symh_nt"
        ] = pd.to_numeric(
            df[
                "symh_nt"
            ],
            errors="coerce",
        )

        for date, g in df.groupby(
            "date"
        ):
            s = g[
                "symh_nt"
            ]

            rows.append({
                "date":
                    date,
                "omni_symh_min_nt":
                    s.min(
                        skipna=True
                    ),
                "omni_symh_mean_nt":
                    s.mean(
                        skipna=True
                    ),
                "n_symh_5min":
                    int(
                        s.notna().sum()
                    ),
            })

    return pd.DataFrame(
        rows
    )


# =============================================================================
# SILSO DAILY SUNSPOT NUMBER
# =============================================================================

def load_silso_daily():
    path = (
        INDEX_CACHE
        / "SILSO"
        / "SN_d_tot_V2.0.csv"
    )

    download_file(
        SILSO_DAILY_URL,
        path,
        min_bytes=100_000,
    )

    df = pd.read_csv(
        path,
        sep=";",
        header=None,
        names=[
            "year",
            "month",
            "day",
            "decimal_year",
            "ssn",
            "ssn_std",
            "n_obs",
            "definitive",
        ],
    )

    for col in [
        "year",
        "month",
        "day",
        "ssn",
        "ssn_std",
        "n_obs",
        "definitive",
    ]:
        df[
            col
        ] = pd.to_numeric(
            df[
                col
            ],
            errors="coerce",
        )

    df.loc[
        df[
            "ssn"
        ]
        < 0,
        "ssn",
    ] = np.nan

    df[
        "date"
    ] = pd.to_datetime(
        {
            "year":
                df[
                    "year"
                ],
            "month":
                df[
                    "month"
                ],
            "day":
                df[
                    "day"
                ],
        },
        errors="coerce",
    ).dt.date.astype(
        str
    )

    return df[
        [
            "date",
            "ssn",
            "ssn_std",
            "n_obs",
            "definitive",
        ]
    ].rename(
        columns={
            "ssn":
                "silso_ssn",
            "ssn_std":
                "silso_ssn_std",
            "n_obs":
                "silso_n_obs",
            "definitive":
                "silso_definitive",
        }
    )


# =============================================================================
# DAILY TEC FROM COMMON-HOUR SERIES
# =============================================================================

def load_common_hour_daily():
    if not COMMON_HOUR_FILE.is_file():
        raise FileNotFoundError(
            f"Missing common-hour file:\n{COMMON_HOUR_FILE}"
        )

    df = pd.read_csv(
        COMMON_HOUR_FILE
    )

    required = {
        "epoch",
        "date",
        "year",
        "month",
        "doy",
        "f107_sfu",
        "iri_vtec_tecu",
        "pytecgg_veq_tecu",
        "igs_vtec_tecu",
        "madrigal_vtec_tecu",
    }

    missing = (
        required
        - set(
            df.columns
        )
    )

    if missing:
        raise RuntimeError(
            f"Common-hour file missing columns: {sorted(missing)}"
        )

    df[
        "epoch"
    ] = pd.to_datetime(
        df[
            "epoch"
        ],
        utc=True,
        errors="coerce",
    )

    numeric = [
        "f107_sfu",
        "iri_vtec_tecu",
        "pytecgg_veq_tecu",
        "igs_vtec_tecu",
        "madrigal_vtec_tecu",
    ]

    for col in numeric:
        df[
            col
        ] = pd.to_numeric(
            df[
                col
            ],
            errors="coerce",
        )

    # Common-hour differences.
    df[
        "pytecgg_minus_igs_tecu"
    ] = (
        df[
            "pytecgg_veq_tecu"
        ]
        - df[
            "igs_vtec_tecu"
        ]
    )

    df[
        "pytecgg_minus_madrigal_tecu"
    ] = (
        df[
            "pytecgg_veq_tecu"
        ]
        - df[
            "madrigal_vtec_tecu"
        ]
    )

    df[
        "igs_minus_madrigal_tecu"
    ] = (
        df[
            "igs_vtec_tecu"
        ]
        - df[
            "madrigal_vtec_tecu"
        ]
    )

    df[
        "pytecgg_minus_pyiri_tecu"
    ] = (
        df[
            "pytecgg_veq_tecu"
        ]
        - df[
            "iri_vtec_tecu"
        ]
    )

    df[
        "igs_minus_pyiri_tecu"
    ] = (
        df[
            "igs_vtec_tecu"
        ]
        - df[
            "iri_vtec_tecu"
        ]
    )

    df[
        "madrigal_minus_pyiri_tecu"
    ] = (
        df[
            "madrigal_vtec_tecu"
        ]
        - df[
            "iri_vtec_tecu"
        ]
    )

    value_cols = [
        "iri_vtec_tecu",
        "pytecgg_veq_tecu",
        "igs_vtec_tecu",
        "madrigal_vtec_tecu",
        "pytecgg_minus_igs_tecu",
        "pytecgg_minus_madrigal_tecu",
        "igs_minus_madrigal_tecu",
        "pytecgg_minus_pyiri_tecu",
        "igs_minus_pyiri_tecu",
        "madrigal_minus_pyiri_tecu",
    ]

    rows = []

    for (
        date,
        year,
        month,
        doy,
    ), g in df.groupby(
        [
            "date",
            "year",
            "month",
            "doy",
        ]
    ):
        row = {
            "date":
                str(
                    date
                ),
            "year":
                int(
                    year
                ),
            "month":
                int(
                    month
                ),
            "doy":
                int(
                    doy
                ),
            "common_hour_f107_sfu":
                g[
                    "f107_sfu"
                ].median(
                    skipna=True
                ),
        }

        for col in value_cols:
            s = g[
                col
            ]

            base = col.replace(
                "_tecu",
                ""
            )

            row[
                base
                + "_mean_tecu"
            ] = s.mean(
                skipna=True
            )

            row[
                base
                + "_median_tecu"
            ] = s.median(
                skipna=True
            )

            row[
                base
                + "_std_tecu"
            ] = s.std(
                skipna=True,
                ddof=0,
            )

            row[
                base
                + "_n"
            ] = int(
                s.notna().sum()
            )

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    ).sort_values(
        "date"
    )


# =============================================================================
# MERGE EXISTING VALIDATION ERROR METRICS
# =============================================================================

def merge_validation_stats(
    master: pd.DataFrame,
):
    out = master.copy()

    if IGS_DAILY_FILE is not None and IGS_DAILY_FILE.is_file():
        igs = pd.read_csv(
            IGS_DAILY_FILE
        )

        if {
            "date",
            "comparison",
        }.issubset(
            igs.columns
        ):
            keep_metrics = [
                col
                for col in [
                    "bias_tecu",
                    "mae_tecu",
                    "rmse_tecu",
                    "std_residual_tecu",
                    "pearson_r",
                ]
                if col in igs.columns
            ]

            wanted = {
                "PyTECGg_VEq_vs_IGS_station":
                    "pytecgg_veq_vs_igs",
                "PyTECGg_VTEC_vs_IGS_IPP":
                    "pytecgg_vtec_vs_igs",
                "pyOASIS_VTEC_vs_IGS_IPP":
                    "pyoasis_vtec_vs_igs",
            }

            for comparison, prefix in wanted.items():
                g = igs[
                    igs[
                        "comparison"
                    ]
                    == comparison
                ][
                    [
                        "date",
                    ]
                    + keep_metrics
                ].copy()

                if g.empty:
                    continue

                g = g.rename(
                    columns={
                        col:
                            f"{prefix}_{col}"
                        for col in keep_metrics
                    }
                )

                out = out.merge(
                    g,
                    on="date",
                    how="left",
                )

    if MADRIGAL_DAILY_FILE is not None and MADRIGAL_DAILY_FILE.is_file():
        mad = pd.read_csv(
            MADRIGAL_DAILY_FILE
        )

        if {
            "date",
            "comparison",
        }.issubset(
            mad.columns
        ):
            keep_metrics = [
                col
                for col in [
                    "bias_tecu",
                    "mae_tecu",
                    "rmse_tecu",
                    "pearson_r",
                    "match_fraction",
                ]
                if col in mad.columns
            ]

            wanted = {
                "PyTECGg_VEq_vs_Madrigal_station":
                    "pytecgg_veq_vs_madrigal",
                "PyTECGg_VTEC_vs_Madrigal_IPP":
                    "pytecgg_vtec_vs_madrigal",
                "pyOASIS_VTEC_vs_Madrigal_IPP":
                    "pyoasis_vtec_vs_madrigal",
            }

            for comparison, prefix in wanted.items():
                g = mad[
                    mad[
                        "comparison"
                    ]
                    == comparison
                ][
                    [
                        "date",
                    ]
                    + keep_metrics
                ].copy()

                if g.empty:
                    continue

                g = g.rename(
                    columns={
                        col:
                            f"{prefix}_{col}"
                        for col in keep_metrics
                    }
                )

                out = out.merge(
                    g,
                    on="date",
                    how="left",
                )

    return out


# =============================================================================
# SEASONAL FEATURES / FLAGS
# =============================================================================

def add_analysis_features(
    df: pd.DataFrame,
):
    out = df.copy()

    doy = pd.to_numeric(
        out[
            "doy"
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    theta = (
        2.0
        * np.pi
        * doy
        / 365.25
    )

    out[
        "season_sin1"
    ] = np.sin(
        theta
    )

    out[
        "season_cos1"
    ] = np.cos(
        theta
    )

    out[
        "season_sin2"
    ] = np.sin(
        2.0
        * theta
    )

    out[
        "season_cos2"
    ] = np.cos(
        2.0
        * theta
    )

    out[
        "dst_intensity_nt"
    ] = -pd.to_numeric(
        out[
            "omni_dst_min_nt"
        ],
        errors="coerce",
    )

    out[
        "symh_intensity_nt"
    ] = -pd.to_numeric(
        out[
            "omni_symh_min_nt"
        ],
        errors="coerce",
    )

    kp = pd.to_numeric(
        out[
            "omni_kp_max"
        ],
        errors="coerce",
    )

    dst = pd.to_numeric(
        out[
            "omni_dst_min_nt"
        ],
        errors="coerce",
    )

    symh = pd.to_numeric(
        out[
            "omni_symh_min_nt"
        ],
        errors="coerce",
    )

    out[
        "activity_data_complete"
    ] = (
        kp.notna()
        & dst.notna()
        & symh.notna()
    )

    out[
        "quiet_day"
    ] = (
        out[
            "activity_data_complete"
        ]
        & (
            kp
            < QUIET_KP_MAX
        )
        & (
            dst
            > QUIET_DST_MIN_NT
        )
        & (
            symh
            > QUIET_SYMH_MIN_NT
        )
    )

    out[
        "storm_day"
    ] = (
        out[
            "activity_data_complete"
        ]
        & (
            (
                kp
                >= STORM_KP_MAX
            )
            | (
                dst
                <= STORM_DST_MIN_NT
            )
            | (
                symh
                <= STORM_SYMH_MIN_NT
            )
        )
    )

    def activity_class(
        row,
    ):
        if not bool(
            row[
                "activity_data_complete"
            ]
        ):
            return "unknown"

        if bool(
            row[
                "storm_day"
            ]
        ):
            return "storm"

        if bool(
            row[
                "quiet_day"
            ]
        ):
            return "quiet"

        return "active_nonstorm"

    out[
        "activity_class"
    ] = out.apply(
        activity_class,
        axis=1,
    )

    return out


# =============================================================================
# CORRELATION HELPERS
# =============================================================================

def correlation_pair(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
):
    x = pd.to_numeric(
        df[
            x_col
        ],
        errors="coerce",
    )

    y = pd.to_numeric(
        df[
            y_col
        ],
        errors="coerce",
    )

    good = (
        x.notna()
        & y.notna()
    )

    x = x[
        good
    ].to_numpy(
        dtype=float
    )

    y = y[
        good
    ].to_numpy(
        dtype=float
    )

    n = len(
        x
    )

    if (
        n
        < MIN_CORRELATION_N
        or np.nanstd(
            x
        )
        == 0
        or np.nanstd(
            y
        )
        == 0
    ):
        return {
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

    pearson = scipy_stats.pearsonr(
        x,
        y,
    )

    spearman = scipy_stats.spearmanr(
        x,
        y,
    )

    return {
        "n":
            n,
        "pearson_r":
            float(
                pearson.statistic
            ),
        "pearson_p":
            float(
                pearson.pvalue
            ),
        "spearman_rho":
            float(
                spearman.statistic
            ),
        "spearman_p":
            float(
                spearman.pvalue
            ),
    }


def correlation_table(
    df: pd.DataFrame,
    x_cols,
    y_cols,
    subset_label: str,
):
    rows = []

    for y in y_cols:
        if y not in df.columns:
            continue

        for x in x_cols:
            if x not in df.columns:
                continue

            rows.append({
                "subset":
                    subset_label,
                "dependent":
                    y,
                "predictor":
                    x,
                **correlation_pair(
                    df,
                    x,
                    y,
                ),
            })

    return pd.DataFrame(
        rows
    )


def add_bh_fdr(
    table: pd.DataFrame,
    p_col: str,
    q_col: str,
):
    """
    Benjamini-Hochberg false-discovery-rate correction over the finite p-values
    in one correlation table.
    """
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

    # Enforce monotonicity from largest rank downward.
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


def add_correlation_fdr(
    table: pd.DataFrame,
):
    out = add_bh_fdr(
        table,
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
# OLS REGRESSION: F10.7 + SEASON
# =============================================================================

def fit_ols(
    df: pd.DataFrame,
    y_col: str,
    x_cols,
):
    cols = [
        y_col,
    ] + list(
        x_cols
    )

    g = df[
        cols
    ].apply(
        pd.to_numeric,
        errors="coerce",
    ).dropna()

    n = len(
        g
    )

    p = (
        len(
            x_cols
        )
        + 1
    )

    if n < max(
        MIN_REGRESSION_N,
        p
        + 2,
    ):
        return None

    y = g[
        y_col
    ].to_numpy(
        dtype=float
    )

    X_raw = g[
        list(
            x_cols
        )
    ].to_numpy(
        dtype=float
    )

    # Preserve raw scale for interpretable F10.7 slope.
    X = np.column_stack([
        np.ones(
            n
        ),
        X_raw,
    ])

    beta, _, _, _ = np.linalg.lstsq(
        X,
        y,
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

    sse = float(
        np.sum(
            residual
            ** 2
        )
    )

    sst = float(
        np.sum(
            (
                y
                - np.mean(
                    y
                )
            )
            ** 2
        )
    )

    r2 = (
        1.0
        - sse
        / sst
        if sst > 0
        else np.nan
    )

    dof = (
        n
        - p
    )

    sigma2 = (
        sse
        / dof
        if dof > 0
        else np.nan
    )

    try:
        cov = (
            sigma2
            * np.linalg.inv(
                X.T
                @ X
            )
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
        tstat = (
            beta
            / se
        )

    if dof > 0:
        pvals = (
            2.0
            * scipy_stats.t.sf(
                np.abs(
                    tstat
                ),
                df=dof,
            )
        )
    else:
        pvals = np.full(
            p,
            np.nan,
        )

    names = [
        "intercept",
    ] + list(
        x_cols
    )

    coefficients = []

    for i, name in enumerate(
        names
    ):
        coefficients.append({
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
                    tstat[
                        i
                    ]
                )
                if np.isfinite(
                    tstat[
                        i
                    ]
                )
                else np.nan,
            "p_value":
                float(
                    pvals[
                        i
                    ]
                )
                if np.isfinite(
                    pvals[
                        i
                    ]
                )
                else np.nan,
        })

    return {
        "n":
            n,
        "r2":
            r2,
        "rmse":
            float(
                np.sqrt(
                    np.mean(
                        residual
                        ** 2
                    )
                )
            ),
        "coefficients":
            coefficients,
        "index":
            g.index,
        "fitted":
            fitted,
        "residual":
            residual,
    }


def background_regressions(
    master: pd.DataFrame,
    dependent_cols,
):
    model_rows = []
    coefficient_rows = []
    out = master.copy()

    predictors = [
        "omni_f107_sfu",
        "season_sin1",
        "season_cos1",
        "season_sin2",
        "season_cos2",
    ]

    for y_col in dependent_cols:
        if y_col not in out.columns:
            continue

        fit = fit_ols(
            out,
            y_col,
            predictors,
        )

        if fit is None:
            continue

        model_rows.append({
            "dependent":
                y_col,
            "n":
                fit[
                    "n"
                ],
            "r2":
                fit[
                    "r2"
                ],
            "rmse_tecu":
                fit[
                    "rmse"
                ],
            "model":
                "F10.7 + annual/semiannual seasonal harmonics",
        })

        for row in fit[
            "coefficients"
        ]:
            coefficient_rows.append({
                "dependent":
                    y_col,
                **row,
            })

        residual_col = (
            y_col
            + "_background_residual"
        )

        fitted_col = (
            y_col
            + "_background_fitted"
        )

        out[
            residual_col
        ] = np.nan

        out[
            fitted_col
        ] = np.nan

        out.loc[
            fit[
                "index"
            ],
            residual_col,
        ] = fit[
            "residual"
        ]

        out.loc[
            fit[
                "index"
            ],
            fitted_col,
        ] = fit[
            "fitted"
        ]

    return (
        out,
        pd.DataFrame(
            model_rows
        ),
        pd.DataFrame(
            coefficient_rows
        ),
    )


# =============================================================================
# YEARLY EQUAL-DAY SUMMARY
# =============================================================================

def yearly_equal_day_summary(
    df: pd.DataFrame,
):
    desired = [
        "omni_f107_sfu",
        "silso_ssn",
        "omni_kp_max",
        "omni_Ap_daily_from_ap",
        "omni_dst_min_nt",
        "omni_symh_min_nt",
        "iri_vtec_mean_tecu",
        "pytecgg_veq_mean_tecu",
        "igs_vtec_mean_tecu",
        "madrigal_vtec_mean_tecu",
        "pytecgg_minus_igs_mean_tecu",
        "pytecgg_minus_madrigal_mean_tecu",
        "igs_minus_madrigal_mean_tecu",
    ]

    rows = []

    for year, g in df.groupby(
        "year"
    ):
        row = {
            "year":
                int(
                    year
                ),
            "n_days":
                int(
                    len(
                        g
                    )
                ),
            "n_quiet_days":
                int(
                    g[
                        "quiet_day"
                    ].sum()
                ),
            "n_storm_days":
                int(
                    g[
                        "storm_day"
                    ].sum()
                ),
        }

        for col in desired:
            if col not in g.columns:
                continue

            row[
                "mean_"
                + col
            ] = pd.to_numeric(
                g[
                    col
                ],
                errors="coerce",
            ).mean()

            row[
                "median_"
                + col
            ] = pd.to_numeric(
                g[
                    col
                ],
                errors="coerce",
            ).median()

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    ).sort_values(
        "year"
    )


# =============================================================================
# STORM VS QUIET
# =============================================================================

def storm_quiet_summary(
    df: pd.DataFrame,
    dependent_cols,
):
    rows = []

    for y_col in dependent_cols:
        if y_col not in df.columns:
            continue

        for activity in [
            "quiet",
            "active_nonstorm",
            "storm",
        ]:
            g = df[
                df[
                    "activity_class"
                ]
                == activity
            ]

            y = pd.to_numeric(
                g[
                    y_col
                ],
                errors="coerce",
            ).dropna()

            rows.append({
                "dependent":
                    y_col,
                "activity_class":
                    activity,
                "n":
                    len(
                        y
                    ),
                "mean":
                    y.mean(),
                "median":
                    y.median(),
                "std":
                    y.std(
                        ddof=0
                    ),
            })

    return pd.DataFrame(
        rows
    )


# =============================================================================
# REPORT HELPERS
# =============================================================================

def strongest_correlations(
    table: pd.DataFrame,
    n=20,
):
    if table.empty:
        return table

    g = table.copy()

    g[
        "abs_pearson_r"
    ] = np.abs(
        pd.to_numeric(
            g[
                "pearson_r"
            ],
            errors="coerce",
        )
    )

    return (
        g.sort_values(
            [
                "abs_pearson_r",
                "n",
            ],
            ascending=[
                False,
                False,
            ],
        )
        .head(
            n
        )
        .drop(
            columns=[
                "abs_pearson_r",
            ]
        )
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

    INDEX_CACHE.mkdir(
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

    print(
        "=" * 120
    )

    print(
        "KOH2 SOLAR CYCLE 25 / GEOMAGNETIC ACTIVITY ANALYSIS"
    )

    print(
        "=" * 120
    )

    print(
        "Main common-hour input:",
        COMMON_HOUR_FILE,
    )

    # -------------------------------------------------------------------------
    # 1. TEC DAILY TABLE
    # -------------------------------------------------------------------------

    print()
    print(
        "Building daily common-hour TEC metrics ..."
    )

    master = load_common_hour_daily()

    observation_dates = set(
        master[
            "date"
        ].astype(
            str
        )
    )

    print(
        "KOH2 observation days:",
        len(
            master
        ),
    )

    # -------------------------------------------------------------------------
    # 2. SOLAR / GEOMAGNETIC INDICES
    # -------------------------------------------------------------------------

    omni = daily_omni2_indices(
        YEARS
    )

    symh = daily_symh_indices(
        YEARS
    )

    silso = load_silso_daily()

    master = master.merge(
        omni,
        on="date",
        how="left",
    )

    master = master.merge(
        symh,
        on="date",
        how="left",
    )

    master = master.merge(
        silso,
        on="date",
        how="left",
    )

    master = merge_validation_stats(
        master
    )

    master = add_analysis_features(
        master
    )

    # -------------------------------------------------------------------------
    # 3. BACKGROUND SOLAR-CYCLE MODELS
    # -------------------------------------------------------------------------

    vtec_dependents = [
        "pytecgg_veq_mean_tecu",
        "igs_vtec_mean_tecu",
        "madrigal_vtec_mean_tecu",
        "iri_vtec_mean_tecu",
    ]

    reference_difference_dependents = [
        "pytecgg_minus_igs_mean_tecu",
        "pytecgg_minus_madrigal_mean_tecu",
        "igs_minus_madrigal_mean_tecu",
        "pytecgg_minus_pyiri_mean_tecu",
        "igs_minus_pyiri_mean_tecu",
        "madrigal_minus_pyiri_mean_tecu",
    ]

    background_dependents = (
        vtec_dependents
        + reference_difference_dependents
    )

    (
        master,
        background_models,
        background_coefficients,
    ) = background_regressions(
        master,
        background_dependents,
    )

    # -------------------------------------------------------------------------
    # 4. SIMPLE CORRELATIONS, ALL DAYS AND QUIET DAYS
    # -------------------------------------------------------------------------

    solar_predictors = [
        "omni_f107_sfu",
        "silso_ssn",
    ]

    geomag_predictors = [
        "omni_kp_max",
        "omni_Ap_daily_from_ap",
        "dst_intensity_nt",
        "symh_intensity_nt",
    ]

    all_predictors = (
        solar_predictors
        + geomag_predictors
    )

    all_y = (
        vtec_dependents
        + reference_difference_dependents
    )

    corr_all = correlation_table(
        master,
        all_predictors,
        all_y,
        "all_days",
    )

    corr_all = add_correlation_fdr(
        corr_all
    )

    quiet = master[
        master[
            "quiet_day"
        ]
    ].copy()

    corr_quiet = correlation_table(
        quiet,
        solar_predictors,
        all_y,
        "quiet_days",
    )

    corr_quiet = add_correlation_fdr(
        corr_quiet
    )

    # -------------------------------------------------------------------------
    # 5. GEOMAGNETIC CORRELATION OF BACKGROUND-REMOVED RESIDUALS
    # -------------------------------------------------------------------------

    residual_cols = [
        col
        + "_background_residual"
        for col in background_dependents
        if (
            col
            + "_background_residual"
        )
        in master.columns
    ]

    corr_background_residual = correlation_table(
        master,
        geomag_predictors,
        residual_cols,
        "F107_season_removed",
    )

    corr_background_residual = add_correlation_fdr(
        corr_background_residual
    )

    # -------------------------------------------------------------------------
    # 6. VALIDATION ERROR / RMSE VS GEOMAGNETIC ACTIVITY
    # -------------------------------------------------------------------------

    validation_cols = [
        col
        for col in master.columns
        if (
            (
                col.endswith(
                    "_bias_tecu"
                )
                or col.endswith(
                    "_rmse_tecu"
                )
            )
            and (
                "vs_igs"
                in col
                or "vs_madrigal"
                in col
            )
        )
    ]

    corr_validation = correlation_table(
        master,
        (
            solar_predictors
            + geomag_predictors
        ),
        validation_cols,
        "validation_error",
    )

    corr_validation = add_correlation_fdr(
        corr_validation
    )

    # -------------------------------------------------------------------------
    # 7. STORM VS QUIET SUMMARIES
    # -------------------------------------------------------------------------

    storm_quiet_cols = (
        all_y
        + validation_cols
    )

    storm_quiet = storm_quiet_summary(
        master,
        storm_quiet_cols,
    )

    # -------------------------------------------------------------------------
    # 8. YEARLY EQUAL-DAY SUMMARY
    # -------------------------------------------------------------------------

    yearly = yearly_equal_day_summary(
        master
    )

    # -------------------------------------------------------------------------
    # 9. SAVE
    # -------------------------------------------------------------------------

    master_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_solar_geomagnetic_daily_master.csv"
    )

    yearly_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_yearly_equal_day_summary.csv"
    )

    corr_all_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_correlation_all_days.csv"
    )

    corr_quiet_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_correlation_quiet_days.csv"
    )

    model_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_background_regression_F107_season.csv"
    )

    coeff_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_background_regression_coefficients.csv"
    )

    corr_residual_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_background_residual_geomagnetic_correlations.csv"
    )

    corr_validation_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_validation_error_geomagnetic_correlations.csv"
    )

    storm_quiet_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_storm_vs_quiet_summary.csv"
    )

    report_file = (
        OUTPUT_ROOT
        / f"{STATION}_2019_2026_analysis_report.txt"
    )

    master.to_csv(
        master_file,
        index=False,
    )

    yearly.to_csv(
        yearly_file,
        index=False,
    )

    corr_all.to_csv(
        corr_all_file,
        index=False,
    )

    corr_quiet.to_csv(
        corr_quiet_file,
        index=False,
    )

    background_models.to_csv(
        model_file,
        index=False,
    )

    background_coefficients.to_csv(
        coeff_file,
        index=False,
    )

    corr_background_residual.to_csv(
        corr_residual_file,
        index=False,
    )

    corr_validation.to_csv(
        corr_validation_file,
        index=False,
    )

    storm_quiet.to_csv(
        storm_quiet_file,
        index=False,
    )

    strongest_all = strongest_correlations(
        corr_all,
        25,
    )

    strongest_residual = strongest_correlations(
        corr_background_residual,
        25,
    )

    strongest_validation = strongest_correlations(
        corr_validation,
        25,
    )

    with open(
        report_file,
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            "KOH2 SOLAR CYCLE 25 / GEOMAGNETIC ACTIVITY ANALYSIS\n"
        )

        f.write(
            "=" * 100
            + "\n\n"
        )

        f.write(
            "Analysis period: 2019-2026 available KOH2 days\n"
        )

        f.write(
            f"Total daily records: {len(master)}\n"
        )

        f.write(
            f"Quiet days: {int(master['quiet_day'].sum())}\n"
        )

        f.write(
            f"Storm days: {int(master['storm_day'].sum())}\n"
        )

        f.write(
            "Days with incomplete activity indices: "
            f"{int((~master['activity_data_complete']).sum())}\n\n"
        )

        f.write(
            "Quiet criterion: "
            f"Kp_max < {QUIET_KP_MAX}, "
            f"Dst_min > {QUIET_DST_MIN_NT} nT, "
            f"SYM-H_min > {QUIET_SYMH_MIN_NT} nT\n"
        )

        f.write(
            "Storm criterion: "
            f"Kp_max >= {STORM_KP_MAX} OR "
            f"Dst_min <= {STORM_DST_MIN_NT} nT OR "
            f"SYM-H_min <= {STORM_SYMH_MIN_NT} nT\n\n"
        )

        f.write(
            "IMPORTANT: PyIRI uses F10.7 as model forcing. "
            "Therefore PyIRI-vs-F10.7 correlation is not an independent "
            "solar-cycle validation.\n\n"
        )

        f.write(
            "YEARLY EQUAL-DAY SUMMARY\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            yearly.to_string(
                index=False
            )
        )

        f.write(
            "\n\nBACKGROUND REGRESSION MODELS\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            background_models.to_string(
                index=False
            )
        )

        f.write(
            "\n\nSTRONGEST ALL-DAY CORRELATIONS\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            strongest_all.to_string(
                index=False
            )
        )

        f.write(
            "\n\nSTRONGEST GEOMAGNETIC CORRELATIONS AFTER "
            "F10.7+SEASON REMOVAL\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            strongest_residual.to_string(
                index=False
            )
        )

        f.write(
            "\n\nSTRONGEST VALIDATION-ERROR CORRELATIONS\n"
        )

        f.write(
            "-" * 100
            + "\n"
        )

        f.write(
            strongest_validation.to_string(
                index=False
            )
        )

        f.write(
            "\n"
        )

    # -------------------------------------------------------------------------
    # 10. CONSOLE SUMMARY
    # -------------------------------------------------------------------------

    print()
    print(
        "=" * 120
    )

    print(
        "SOLAR CYCLE 25 ANALYSIS COMPLETE"
    )

    print(
        "=" * 120
    )

    print(
        "Daily master:",
        master_file,
    )

    print(
        "Yearly equal-day:",
        yearly_file,
    )

    print(
        "All-day correlations:",
        corr_all_file,
    )

    print(
        "Quiet-day correlations:",
        corr_quiet_file,
    )

    print(
        "Background models:",
        model_file,
    )

    print(
        "Background coefficients:",
        coeff_file,
    )

    print(
        "Background residual vs geomagnetic:",
        corr_residual_file,
    )

    print(
        "Validation errors vs activity:",
        corr_validation_file,
    )

    print(
        "Storm vs quiet:",
        storm_quiet_file,
    )

    print(
        "Report:",
        report_file,
    )

    print()
    print(
        "YEARLY EQUAL-DAY SUMMARY"
    )

    print(
        "-" * 120
    )

    print(
        yearly.to_string(
            index=False
        )
    )

    print()
    print(
        "BACKGROUND MODEL FIT"
    )

    print(
        "-" * 120
    )

    print(
        background_models.to_string(
            index=False
        )
    )

    print()
    print(
        "STRONGEST GEOMAGNETIC ASSOCIATIONS AFTER F10.7+SEASON REMOVAL"
    )

    print(
        "-" * 120
    )

    print(
        strongest_residual.head(
            20
        ).to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
