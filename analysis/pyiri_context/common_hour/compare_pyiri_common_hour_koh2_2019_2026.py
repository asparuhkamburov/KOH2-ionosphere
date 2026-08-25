from __future__ import annotations

r"""
KOH2 COMMON-HOUR COMPARISON:
PyIRI climatology vs PyTECGg VEq vs IGS Final GIM vs Madrigal, 2019-2026
============================================================================

PURPOSE
-------
Place the four station-level VTEC series on the SAME hourly KOH2 epochs:

    1. PyIRI climatological VTEC (hourly)
    2. PyTECGg VEq at KOH2
    3. IGS Final GIM interpolated at KOH2
    4. Madrigal nearest available VTEC grid cell around KOH2

This is the cleanest bridge between the validation work and the later
Solar Cycle 25 analysis.

IMPORTANT
---------
PyIRI is not an independent solar-activity test because daily F10.7 is an
input to the model.  Here it is used as an empirical climatological reference.

Residual convention in pairwise statistics:
    external_or_method_value - PyIRI

Madrigal matching:
    nearest map epoch within 180 s
    nearest available spatial cell within 80 km

PyTECGg matching:
    nearest VEq epoch within 60 s of the exact hourly PyIRI epoch

OUTPUT
------
E:\KOH2data\TEC_COMMON_HOUR_PYIRI_2019_2026

    KOH2_2019_2026_common_hour_values.csv
    KOH2_2019_2026_common_hour_daily_statistics.csv
    KOH2_2019_2026_common_hour_yearly_equal_day_summary.csv
    KOH2_2019_2026_common_hour_monthly_summary.csv
    KOH2_2019_2026_common_hour_availability.csv

The values CSV is also intended to become the base table for the later
F10.7 / SSN / Kp / Ap / Dst / SYM-H analysis.
"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import importlib.util
import re

import numpy as np
import pandas as pd


# =============================================================================
# SETTINGS
# =============================================================================

BASE_ROOT = Path(".")
YEARS = list(range(2019, 2027))
STATION = "KOH2"
SELECTED_DOY = None

KOH2_LAT = -62.64008176
KOH2_LON = -60.36376872

PYIRI_ROOT = (
    BASE_ROOT
    / "TEC_REFERENCE_PYIRI_2019_2026"
)

PYIRI_HOURLY_ROOT = (
    PYIRI_ROOT
    / "HOURLY"
)

MADRIGAL_REGIONAL_ROOT = (
    BASE_ROOT
    / "TEC_VALIDATION_MADRIGAL_2019_2026"
    / "_MADRIGAL_CACHE"
    / "VTEC_REGIONAL_LOCAL"
)

IGS_CACHE_ROOT = (
    BASE_ROOT
    / "TEC_VALIDATION_IGS_2019_2026"
    / "_IGS_IONEX"
)

OUTPUT_ROOT = (
    BASE_ROOT
    / "TEC_COMMON_HOUR_PYIRI_2019_2026"
)

PYTECGG_MAX_TIME_S = 60.0
MADRIGAL_MAX_TIME_S = 180.0
MADRIGAL_MAX_DISTANCE_KM = 80.0


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Place KOH2 PyIRI, PyTECGg VEq, IGS Final GIM and Madrigal VTEC "
            "on common hourly epochs."
        )
    )

    parser.add_argument(
        "--data-root",
        required=True,
        type=Path,
        help="Root containing YEAR/MM/DD KOH2 production data.",
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: "
            "<data-root>/TEC_COMMON_HOUR_PYIRI_2019_2026"
        ),
    )

    parser.add_argument(
        "--pyiri-hourly-root",
        type=Path,
        default=None,
        help=(
            "Directory containing HOURLY/YEAR/KOH2_YEAR_DDD_"
            "PyIRI_URSI_hourly.csv. Default: "
            "<data-root>/TEC_REFERENCE_PYIRI_2019_2026/HOURLY"
        ),
    )

    parser.add_argument(
        "--madrigal-regional-root",
        type=Path,
        default=None,
        help=(
            "Directory containing VTEC_REGIONAL_LOCAL/YEAR/"
            "madrigal_regional_local_YEAR_DDD.parquet. Default: "
            "<data-root>/TEC_VALIDATION_MADRIGAL_2019_2026/"
            "_MADRIGAL_CACHE/VTEC_REGIONAL_LOCAL"
        ),
    )

    parser.add_argument(
        "--igs-cache-root",
        type=Path,
        default=None,
        help=(
            "Existing IGS Final IONEX cache root. Default: "
            "<data-root>/TEC_VALIDATION_IGS_2019_2026/_IGS_IONEX"
        ),
    )

    parser.add_argument(
        "--year",
        type=int,
        action="append",
        default=None,
        help=(
            "Year to process. May be supplied more than once. "
            "Default: 2019 through 2026."
        ),
    )

    parser.add_argument(
        "--date",
        type=str,
        default=None,
        help="Optional representative date in YYYY-MM-DD format.",
    )

    return parser.parse_args()


def configure_runtime(args: argparse.Namespace) -> None:
    global BASE_ROOT, YEARS, SELECTED_DOY
    global PYIRI_ROOT, PYIRI_HOURLY_ROOT
    global MADRIGAL_REGIONAL_ROOT, IGS_CACHE_ROOT, OUTPUT_ROOT

    BASE_ROOT = args.data_root.resolve()

    OUTPUT_ROOT = (
        args.output_root.resolve()
        if args.output_root is not None
        else BASE_ROOT / "TEC_COMMON_HOUR_PYIRI_2019_2026"
    )

    PYIRI_HOURLY_ROOT = (
        args.pyiri_hourly_root.resolve()
        if args.pyiri_hourly_root is not None
        else (
            BASE_ROOT
            / "TEC_REFERENCE_PYIRI_2019_2026"
            / "HOURLY"
        )
    )

    PYIRI_ROOT = PYIRI_HOURLY_ROOT.parent

    MADRIGAL_REGIONAL_ROOT = (
        args.madrigal_regional_root.resolve()
        if args.madrigal_regional_root is not None
        else (
            BASE_ROOT
            / "TEC_VALIDATION_MADRIGAL_2019_2026"
            / "_MADRIGAL_CACHE"
            / "VTEC_REGIONAL_LOCAL"
        )
    )

    IGS_CACHE_ROOT = (
        args.igs_cache_root.resolve()
        if args.igs_cache_root is not None
        else (
            BASE_ROOT
            / "TEC_VALIDATION_IGS_2019_2026"
            / "_IGS_IONEX"
        )
    )

    YEARS = (
        [int(y) for y in args.year]
        if args.year
        else list(range(2019, 2027))
    )

    SELECTED_DOY = None

    if args.date is not None:
        try:
            selected = datetime.strptime(
                args.date,
                "%Y-%m-%d",
            )
        except ValueError as exc:
            raise ValueError(
                "--date must use YYYY-MM-DD format"
            ) from exc

        if args.year and selected.year not in YEARS:
            raise ValueError(
                "--date year must be included in --year"
            )

        YEARS = [selected.year]
        SELECTED_DOY = int(selected.strftime("%j"))


# =============================================================================
# LOAD VALIDATED PUBLICATION IGS CORE
# =============================================================================

def load_igs_core():
    filename = "validate_tec_igs_2019_2026.py"

    candidates = [
        Path.cwd() / filename,
        Path(__file__).resolve().parent / filename,
    ]

    path = next(
        (
            p
            for p in candidates
            if p.is_file()
        ),
        None,
    )

    if path is None:
        raise FileNotFoundError(
            "\nCould not find:\n"
            f"    {filename}\n\n"
            "Place this script in the same directory as the validated publication "
            "IGS script."
        )

    spec = importlib.util.spec_from_file_location(
        "koh2_igs_v4_core",
        path,
    )

    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"Could not import IGS validation core: {path}"
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    print(
        "IGS validation core:",
        path,
    )

    return module


CORE = load_igs_core()

configure_year = CORE.configure_year
download_ionex = CORE.download_ionex
read_ionex = CORE.read_ionex
interpolate_ionex = CORE.interpolate_ionex
calculate_statistics = CORE.calculate_statistics
load_pytecgg = CORE.load_pytecgg


# =============================================================================
# HELPERS
# =============================================================================

def date_from_year_doy(
    year: int,
    doy: int,
):
    return datetime.strptime(
        f"{year}-{doy:03d}",
        "%Y-%j",
    ).replace(
        tzinfo=timezone.utc
    )


def great_circle_km(
    lat1,
    lon1,
    lat2,
    lon2,
):
    r = 6371.0

    p1 = np.radians(
        lat1
    )

    p2 = np.radians(
        lat2
    )

    dp = np.radians(
        lat2
        - lat1
    )

    dl = np.radians(
        (
            (
                lon2
                - lon1
                + 180.0
            )
            % 360.0
        )
        - 180.0
    )

    a = (
        np.sin(
            dp
            / 2.0
        )
        ** 2
        + np.cos(
            p1
        )
        * np.cos(
            p2
        )
        * np.sin(
            dl
            / 2.0
        )
        ** 2
    )

    return (
        2.0
        * r
        * np.arctan2(
            np.sqrt(
                a
            ),
            np.sqrt(
                1.0
                - a
            ),
        )
    )


def wrap_lon(
    lon,
):
    lon = np.asarray(
        lon,
        dtype=float,
    )

    return (
        (
            lon
            + 180.0
        )
        % 360.0
    ) - 180.0


def extract_doy_from_pytecgg(
    path: Path,
    year: int,
):
    m = re.search(
        r"_(\d{3})_(\d{4})_",
        path.name,
    )

    if not m:
        return None

    doy = int(
        m.group(
            1
        )
    )

    file_year = int(
        m.group(
            2
        )
    )

    if file_year != year:
        return None

    return doy


def discover_pytecgg(
    year: int,
):
    root = (
        BASE_ROOT
        / str(
            year
        )
    )

    result = {}

    if not root.is_dir():
        return result

    for path in root.glob(
        rf"*\*\PyTECGg_OUTPUT\{STATION}_*_{year}_PyTECGg_VEQ.parquet"
    ):
        doy = extract_doy_from_pytecgg(
            path,
            year,
        )

        if doy is not None:
            result[
                doy
            ] = path

    return result


def discover_pyiri(
    year: int,
):
    folder = (
        PYIRI_HOURLY_ROOT
        / str(
            year
        )
    )

    result = {}

    if not folder.is_dir():
        return result

    for path in folder.glob(
        f"{STATION}_{year}_*_PyIRI_URSI_hourly.csv"
    ):
        m = re.search(
            rf"{STATION}_{year}_(\d{{3}})_PyIRI_URSI_hourly\.csv$",
            path.name,
        )

        if not m:
            continue

        result[
            int(
                m.group(
                    1
                )
            )
        ] = path

    return dict(
        sorted(
            result.items()
        )
    )


def madrigal_file(
    year: int,
    doy: int,
):
    path = (
        MADRIGAL_REGIONAL_ROOT
        / str(
            year
        )
        / (
            f"madrigal_regional_local_{year}_{doy:03d}.parquet"
        )
    )

    return (
        path
        if path.is_file()
        else None
    )


# =============================================================================
# LOAD EACH SERIES
# =============================================================================

def load_pyiri_hourly(
    path: Path,
):
    df = pd.read_csv(
        path
    )

    required = {
        "epoch",
        "f107_sfu",
        "iri_vtec_tecu",
    }

    missing = (
        required
        - set(
            df.columns
        )
    )

    if missing:
        raise RuntimeError(
            f"{path} missing PyIRI columns: {sorted(missing)}"
        )

    df = df.copy()

    df[
        "epoch"
    ] = pd.to_datetime(
        df[
            "epoch"
        ],
        utc=True,
        errors="coerce",
    )

    df[
        "f107_sfu"
    ] = pd.to_numeric(
        df[
            "f107_sfu"
        ],
        errors="coerce",
    )

    df[
        "iri_vtec_tecu"
    ] = pd.to_numeric(
        df[
            "iri_vtec_tecu"
        ],
        errors="coerce",
    )

    return (
        df.dropna(
            subset=[
                "epoch",
                "iri_vtec_tecu",
            ]
        )
        .sort_values(
            "epoch"
        )
        .reset_index(
            drop=True
        )
    )


def match_pytecgg_veq(
    hourly: pd.DataFrame,
    path: Path | None,
):
    out = hourly.copy()

    out[
        "pytecgg_veq_tecu"
    ] = np.nan

    out[
        "pytecgg_time_offset_s"
    ] = np.nan

    if path is None:
        return out

    pt = load_pytecgg(
        path
    )

    station = (
        pt[
            [
                "epoch",
                "veq",
            ]
        ]
        .dropna(
            subset=[
                "epoch",
                "veq",
            ]
        )
        .drop_duplicates(
            subset=[
                "epoch",
            ]
        )
        .copy()
    )

    station[
        "epoch"
    ] = pd.to_datetime(
        station[
            "epoch"
        ],
        utc=True,
        errors="coerce",
    )

    station[
        "veq"
    ] = pd.to_numeric(
        station[
            "veq"
        ],
        errors="coerce",
    )

    station = (
        station.dropna()
        .sort_values(
            "epoch"
        )
        .reset_index(
            drop=True
        )
    )

    if station.empty:
        return out

    left = out[
        [
            "epoch",
        ]
    ].copy()

    right = station.rename(
        columns={
            "epoch":
                "pytecgg_epoch",
        }
    )

    matched = pd.merge_asof(
        left.sort_values(
            "epoch"
        ),
        right.sort_values(
            "pytecgg_epoch"
        ),
        left_on="epoch",
        right_on="pytecgg_epoch",
        direction="nearest",
        tolerance=pd.Timedelta(
            seconds=PYTECGG_MAX_TIME_S
        ),
    )

    out[
        "pytecgg_veq_tecu"
    ] = matched[
        "veq"
    ].to_numpy()

    dt = (
        matched[
            "pytecgg_epoch"
        ]
        - matched[
            "epoch"
        ]
    ).dt.total_seconds().abs()

    out[
        "pytecgg_time_offset_s"
    ] = dt.to_numpy()

    return out


def add_igs_station(
    hourly: pd.DataFrame,
    year: int,
    doy: int,
):
    out = hourly.copy()

    ionex_path = download_ionex(
        year,
        doy,
    )

    grid = read_ionex(
        ionex_path
    )

    out[
        "igs_vtec_tecu"
    ] = interpolate_ionex(
        grid,
        out[
            "epoch"
        ],
        np.full(
            len(
                out
            ),
            KOH2_LAT,
        ),
        np.full(
            len(
                out
            ),
            KOH2_LON,
        ),
    )

    return out


def load_madrigal(
    path: Path | None,
):
    if path is None:
        return pd.DataFrame()

    df = pd.read_parquet(
        path
    )

    if df.empty:
        return df

    df = df.copy()

    df[
        "epoch"
    ] = pd.to_datetime(
        df[
            "epoch"
        ],
        utc=True,
        errors="coerce",
    )

    for col in [
        "gdlat",
        "glon",
        "tec",
    ]:
        df[
            col
        ] = pd.to_numeric(
            df[
                col
            ],
            errors="coerce",
        )

    if "dtec" in df.columns:
        df[
            "dtec"
        ] = pd.to_numeric(
            df[
                "dtec"
            ],
            errors="coerce",
        )
    else:
        df[
            "dtec"
        ] = np.nan

    df[
        "glon"
    ] = wrap_lon(
        df[
            "glon"
        ].to_numpy(
            dtype=float
        )
    )

    return (
        df.dropna(
            subset=[
                "epoch",
                "gdlat",
                "glon",
                "tec",
            ]
        )
        .sort_values(
            "epoch"
        )
        .reset_index(
            drop=True
        )
    )


def match_madrigal_station(
    hourly: pd.DataFrame,
    madrigal: pd.DataFrame,
):
    out = hourly.copy()

    out[
        "madrigal_vtec_tecu"
    ] = np.nan

    out[
        "madrigal_dtec_tecu"
    ] = np.nan

    out[
        "madrigal_time_offset_s"
    ] = np.nan

    out[
        "madrigal_distance_km"
    ] = np.nan

    if madrigal.empty:
        return out

    map_times = pd.DatetimeIndex(
        madrigal[
            "epoch"
        ].dropna().unique()
    ).sort_values()

    if len(
        map_times
    ) == 0:
        return out

    map_ns = map_times.as_unit(
        "ns"
    ).asi8

    obs_times = pd.DatetimeIndex(
        out[
            "epoch"
        ]
    )

    obs_ns = obs_times.as_unit(
        "ns"
    ).asi8

    pos = np.searchsorted(
        map_ns,
        obs_ns,
        side="left",
    )

    lo = np.clip(
        pos - 1,
        0,
        len(
            map_ns
        ) - 1,
    )

    hi = np.clip(
        pos,
        0,
        len(
            map_ns
        ) - 1,
    )

    dlo = np.abs(
        obs_ns
        - map_ns[
            lo
        ]
    )

    dhi = np.abs(
        obs_ns
        - map_ns[
            hi
        ]
    )

    choose_hi = (
        dhi
        < dlo
    )

    chosen = np.where(
        choose_hi,
        hi,
        lo,
    )

    dt_s = (
        np.abs(
            obs_ns
            - map_ns[
                chosen
            ]
        )
        / 1e9
    )

    for i in range(
        len(
            out
        )
    ):
        if (
            not np.isfinite(
                dt_s[
                    i
                ]
            )
            or dt_s[
                i
            ]
            > MADRIGAL_MAX_TIME_S
        ):
            continue

        mt = map_times[
            chosen[
                i
            ]
        ]

        grid = madrigal[
            madrigal[
                "epoch"
            ]
            == mt
        ]

        if grid.empty:
            continue

        distances = great_circle_km(
            KOH2_LAT,
            KOH2_LON,
            grid[
                "gdlat"
            ].to_numpy(
                dtype=float
            ),
            grid[
                "glon"
            ].to_numpy(
                dtype=float
            ),
        )

        if len(
            distances
        ) == 0:
            continue

        j = int(
            np.nanargmin(
                distances
            )
        )

        distance = float(
            distances[
                j
            ]
        )

        if distance > MADRIGAL_MAX_DISTANCE_KM:
            continue

        row = grid.iloc[
            j
        ]

        out.loc[
            out.index[
                i
            ],
            "madrigal_vtec_tecu",
        ] = float(
            row[
                "tec"
            ]
        )

        if pd.notna(
            row[
                "dtec"
            ]
        ):
            out.loc[
                out.index[
                    i
                ],
                "madrigal_dtec_tecu",
            ] = float(
                row[
                    "dtec"
                ]
            )

        out.loc[
            out.index[
                i
            ],
            "madrigal_time_offset_s",
        ] = float(
            dt_s[
                i
            ]
        )

        out.loc[
            out.index[
                i
            ],
            "madrigal_distance_km",
        ] = distance

    return out


# =============================================================================
# STATISTICS
# =============================================================================

PAIRWISE_COLUMNS = {
    "PyTECGg_VEq_minus_PyIRI": (
        "pytecgg_veq_tecu",
        "iri_vtec_tecu",
    ),
    "IGS_minus_PyIRI": (
        "igs_vtec_tecu",
        "iri_vtec_tecu",
    ),
    "Madrigal_minus_PyIRI": (
        "madrigal_vtec_tecu",
        "iri_vtec_tecu",
    ),
    "PyTECGg_VEq_minus_IGS": (
        "pytecgg_veq_tecu",
        "igs_vtec_tecu",
    ),
    "PyTECGg_VEq_minus_Madrigal": (
        "pytecgg_veq_tecu",
        "madrigal_vtec_tecu",
    ),
    "IGS_minus_Madrigal": (
        "igs_vtec_tecu",
        "madrigal_vtec_tecu",
    ),
}


def daily_pairwise_rows(
    day: pd.DataFrame,
    year: int,
    doy: int,
):
    rows = []

    date = date_from_year_doy(
        year,
        doy,
    )

    for name, (
        a_col,
        b_col,
    ) in PAIRWISE_COLUMNS.items():
        valid = day[
            [
                a_col,
                b_col,
            ]
        ].dropna()

        stats = calculate_statistics(
            valid[
                a_col
            ],
            valid[
                b_col
            ],
        )

        rows.append({
            "date":
                date.date().isoformat(),
            "year":
                year,
            "month":
                date.month,
            "doy":
                doy,
            "comparison":
                name,
            **stats,
            "median_a_tecu":
                (
                    float(
                        valid[
                            a_col
                        ].median()
                    )
                    if len(
                        valid
                    )
                    else np.nan
                ),
            "median_b_tecu":
                (
                    float(
                        valid[
                            b_col
                        ].median()
                    )
                    if len(
                        valid
                    )
                    else np.nan
                ),
            "f107_sfu":
                (
                    float(
                        day[
                            "f107_sfu"
                        ].dropna().median()
                    )
                    if day[
                        "f107_sfu"
                    ].notna().any()
                    else np.nan
                ),
            "n_all4_common_hours":
                int(
                    day[
                        [
                            "iri_vtec_tecu",
                            "pytecgg_veq_tecu",
                            "igs_vtec_tecu",
                            "madrigal_vtec_tecu",
                        ]
                    ]
                    .dropna()
                    .shape[
                        0
                    ]
                ),
        })

    return rows


def yearly_equal_day_summary(
    daily: pd.DataFrame,
):
    if daily.empty:
        return pd.DataFrame()

    rows = []

    for (
        year,
        comparison,
    ), g in daily.groupby(
        [
            "year",
            "comparison",
        ]
    ):
        good = g[
            pd.to_numeric(
                g[
                    "n"
                ],
                errors="coerce",
            )
            > 0
        ].copy()

        rows.append({
            "year":
                int(
                    year
                ),
            "comparison":
                comparison,
            "n_days":
                int(
                    len(
                        good
                    )
                ),
            "n_points_total":
                int(
                    pd.to_numeric(
                        good[
                            "n"
                        ],
                        errors="coerce",
                    ).fillna(
                        0
                    ).sum()
                ),
            "mean_daily_bias_tecu":
                good[
                    "bias_tecu"
                ].mean(),
            "median_daily_bias_tecu":
                good[
                    "bias_tecu"
                ].median(),
            "mean_daily_rmse_tecu":
                good[
                    "rmse_tecu"
                ].mean(),
            "median_daily_rmse_tecu":
                good[
                    "rmse_tecu"
                ].median(),
            "median_daily_pearson_r":
                good[
                    "pearson_r"
                ].median(),
            "median_daily_f107_sfu":
                good[
                    "f107_sfu"
                ].median(),
            "median_all4_common_hours_per_day":
                good[
                    "n_all4_common_hours"
                ].median(),
        })

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "comparison",
        ]
    )


def monthly_summary(
    daily: pd.DataFrame,
):
    if daily.empty:
        return pd.DataFrame()

    rows = []

    for (
        year,
        month,
        comparison,
    ), g in daily.groupby(
        [
            "year",
            "month",
            "comparison",
        ]
    ):
        good = g[
            pd.to_numeric(
                g[
                    "n"
                ],
                errors="coerce",
            )
            > 0
        ].copy()

        rows.append({
            "year":
                int(
                    year
                ),
            "month":
                int(
                    month
                ),
            "comparison":
                comparison,
            "n_days":
                int(
                    len(
                        good
                    )
                ),
            "mean_daily_bias_tecu":
                good[
                    "bias_tecu"
                ].mean(),
            "median_daily_bias_tecu":
                good[
                    "bias_tecu"
                ].median(),
            "mean_daily_rmse_tecu":
                good[
                    "rmse_tecu"
                ].mean(),
            "median_daily_rmse_tecu":
                good[
                    "rmse_tecu"
                ].median(),
            "median_daily_pearson_r":
                good[
                    "pearson_r"
                ].median(),
            "median_f107_sfu":
                good[
                    "f107_sfu"
                ].median(),
        })

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "month",
            "comparison",
        ]
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()
    configure_runtime(args)

    # The imported publication IGS core owns configure_year()/download_ionex().
    # Redirect its cache globals to the selected validated IONEX cache.
    CORE.BASE_ROOT = BASE_ROOT
    CORE.MASTER_OUTPUT = OUTPUT_ROOT / "_UNUSED_IGS_MASTER"
    CORE.CACHE_ROOT = IGS_CACHE_ROOT
    CORE.BY_YEAR_OUTPUT = OUTPUT_ROOT / "_UNUSED_IGS_BY_YEAR"

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
        260,
    )

    print(
        "=" * 118
    )

    print(
        "KOH2 COMMON-HOUR PyIRI / PyTECGg / IGS / MADRIGAL COMPARISON"
    )

    print(
        "=" * 118
    )

    print(
        f"KOH2: lat={KOH2_LAT:.8f}, lon={KOH2_LON:.8f}"
    )

    print(
        f"PyTECGg tolerance: {PYTECGG_MAX_TIME_S:.0f} s"
    )

    print(
        f"Madrigal tolerance: {MADRIGAL_MAX_TIME_S:.0f} s, "
        f"{MADRIGAL_MAX_DISTANCE_KM:.0f} km"
    )

    print(
        "Data root:",
        BASE_ROOT,
    )

    print(
        "PyIRI hourly root:",
        PYIRI_HOURLY_ROOT,
    )

    print(
        "Madrigal regional root:",
        MADRIGAL_REGIONAL_ROOT,
    )

    print(
        "IGS IONEX cache:",
        IGS_CACHE_ROOT,
    )

    print(
        "Output:",
        OUTPUT_ROOT,
    )

    print(
        "Selected DOY:",
        (
            f"{SELECTED_DOY:03d}"
            if SELECTED_DOY is not None
            else "ALL AVAILABLE"
        ),
    )

    print(
        "=" * 118
    )

    all_values = []
    daily_rows = []
    availability_rows = []

    for year in YEARS:
        pyiri_files = discover_pyiri(
            year
        )

        pt_files = discover_pytecgg(
            year
        )

        if SELECTED_DOY is not None:
            pyiri_files = {
                doy: path
                for doy, path in pyiri_files.items()
                if doy == SELECTED_DOY
            }

            pt_files = {
                doy: path
                for doy, path in pt_files.items()
                if doy == SELECTED_DOY
            }

        print()
        print(
            "=" * 118
        )

        print(
            f"{year}: PyIRI days={len(pyiri_files)}, "
            f"PyTECGg days={len(pt_files)}"
        )

        print(
            "=" * 118
        )

        configure_year(
            year
        )

        for i, (
            doy,
            iri_path,
        ) in enumerate(
            pyiri_files.items(),
            1,
        ):
            print(
                f"{i:3d}/{len(pyiri_files):3d}  "
                f"{year} DOY {doy:03d}"
            )

            status = {
                "year":
                    year,
                "doy":
                    doy,
                "date":
                    date_from_year_doy(
                        year,
                        doy,
                    ).date().isoformat(),
                "pyiri":
                    True,
                "pytecgg":
                    doy in pt_files,
                "madrigal":
                    madrigal_file(
                        year,
                        doy,
                    ) is not None,
                "status":
                    "FAILED",
            }

            try:
                day = load_pyiri_hourly(
                    iri_path
                )

                day = match_pytecgg_veq(
                    day,
                    pt_files.get(
                        doy
                    ),
                )

                day = add_igs_station(
                    day,
                    year,
                    doy,
                )

                mad = load_madrigal(
                    madrigal_file(
                        year,
                        doy,
                    )
                )

                day = match_madrigal_station(
                    day,
                    mad,
                )

                date = date_from_year_doy(
                    year,
                    doy,
                )

                day[
                    "date"
                ] = date.date().isoformat()

                day[
                    "year"
                ] = year

                day[
                    "month"
                ] = date.month

                day[
                    "doy"
                ] = doy

                day[
                    "all4_available"
                ] = (
                    day[
                        [
                            "iri_vtec_tecu",
                            "pytecgg_veq_tecu",
                            "igs_vtec_tecu",
                            "madrigal_vtec_tecu",
                        ]
                    ].notna().all(
                        axis=1
                    )
                )

                all_values.append(
                    day
                )

                rows = daily_pairwise_rows(
                    day,
                    year,
                    doy,
                )

                daily_rows.extend(
                    rows
                )

                status.update({
                    "status":
                        "OK",
                    "n_hourly":
                        len(
                            day
                        ),
                    "n_pytecgg":
                        int(
                            day[
                                "pytecgg_veq_tecu"
                            ].notna().sum()
                        ),
                    "n_igs":
                        int(
                            day[
                                "igs_vtec_tecu"
                            ].notna().sum()
                        ),
                    "n_madrigal":
                        int(
                            day[
                                "madrigal_vtec_tecu"
                            ].notna().sum()
                        ),
                    "n_all4":
                        int(
                            day[
                                "all4_available"
                            ].sum()
                        ),
                })

                print(
                    f"    hours: PyIRI={len(day)}, "
                    f"PyTECGg={status['n_pytecgg']}, "
                    f"IGS={status['n_igs']}, "
                    f"Madrigal={status['n_madrigal']}, "
                    f"ALL4={status['n_all4']}"
                )

            except KeyboardInterrupt:
                raise

            except Exception as exc:
                status[
                    "error"
                ] = repr(
                    exc
                )

                print(
                    "    [ERROR]",
                    repr(
                        exc
                    ),
                )

            availability_rows.append(
                status
            )

    values = (
        pd.concat(
            all_values,
            ignore_index=True,
        )
        if all_values
        else pd.DataFrame()
    )

    daily = pd.DataFrame(
        daily_rows
    )

    availability = pd.DataFrame(
        availability_rows
    )

    yearly = yearly_equal_day_summary(
        daily
    )

    monthly = monthly_summary(
        daily
    )

    values_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "common_hour_values.csv"
        )
    )

    daily_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "common_hour_daily_statistics.csv"
        )
    )

    yearly_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "common_hour_yearly_equal_day_summary.csv"
        )
    )

    monthly_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "common_hour_monthly_summary.csv"
        )
    )

    availability_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "common_hour_availability.csv"
        )
    )

    values.to_csv(
        values_file,
        index=False,
    )

    daily.to_csv(
        daily_file,
        index=False,
    )

    yearly.to_csv(
        yearly_file,
        index=False,
    )

    monthly.to_csv(
        monthly_file,
        index=False,
    )

    availability.to_csv(
        availability_file,
        index=False,
    )

    print()
    print(
        "=" * 118
    )

    print(
        "COMMON-HOUR COMPARISON COMPLETE"
    )

    print(
        "=" * 118
    )

    print(
        "Values:",
        values_file,
    )

    print(
        "Daily:",
        daily_file,
    )

    print(
        "Yearly equal-day:",
        yearly_file,
    )

    print(
        "Monthly:",
        monthly_file,
    )

    print(
        "Availability:",
        availability_file,
    )

    if not yearly.empty:
        print()
        print(
            "YEARLY EQUAL-DAY SUMMARY"
        )

        print(
            "-" * 118
        )

        print(
            yearly.to_string(
                index=False
            )
        )

        print(
            "-" * 118
        )


if __name__ == "__main__":
    main()
