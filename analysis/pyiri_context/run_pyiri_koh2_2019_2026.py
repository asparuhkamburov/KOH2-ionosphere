from __future__ import annotations

r"""
KOH2 PyIRI / IRI CLIMATOLOGICAL VTEC, 2019-2026
================================================

PURPOSE
-------
Generate an empirical/climatological IRI VTEC reference at KOH2 for the
available GNSS observation days.

Implementation:
    PyIRI (pure-Python implementation of International Reference Ionosphere)

IMPORTANT SCIENTIFIC LABEL
--------------------------
This is NOT the native CCMC/Fortran IRI-2020 executable.

For the monograph, describe it as:

    "PyIRI implementation of IRI climatology (URSI option)"

or, more fully:

    "IRI climatological VTEC generated with the PyIRI implementation,
     using URSI foF2 coefficients and daily F10.7 forcing."

The model is evaluated at the fixed KOH2 geographic coordinates.

VTEC:
    electron density integrated vertically from 90 km to 2000 km

Time:
    1-hour UT sampling: 00:00 ... 23:00

Solar forcing:
    NASA OMNI daily observed F10.7.

    PyIRI's IRI_density_1day() API accepts one user-provided F10.7 value.
    Therefore the production workflow uses the NASA OMNI daily F10.7 series
    directly and records that source in every output row.

Coefficient option:
    URSI (ccir_or_ursi = 1)

Why URSI:
    KOH2 is a high-latitude maritime/Antarctic site and IRI documentation
    traditionally recommends URSI when one global/oceanic option is desired.

PRODUCTION RUN
--------------
By default the script processes all available KOH2 observation days from
2019 through 2026. Runtime paths and representative-date selection are
provided through command-line arguments.

OUTPUT
------
E:\KOH2data\TEC_REFERENCE_PYIRI_2019_2026

    HOURLY\
        KOH2_YYYY_DDD_PyIRI_URSI_hourly.csv

    KOH2_2019_2026_PyIRI_daily_statistics.csv
    KOH2_2019_2026_PyIRI_availability.csv
    KOH2_2019_2026_PyIRI_yearly_equal_day_summary.csv
    KOH2_2019_2026_PyIRI_monthly_summary.csv

INSTALL
-------
In pytecgg_env:

    python -m pip install PyIRI requests

"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import importlib.metadata
import re
import sys
import time

import numpy as np
import pandas as pd
import requests


# =============================================================================
# SETTINGS
# =============================================================================

BASE_ROOT = Path(".")
YEARS = list(range(2019, 2027))
STATION = "KOH2"
SELECTED_DOY = None

# KOH2 coordinates established in the GNSS workflow.
KOH2_LAT = -62.64008176
KOH2_LON = -60.36376872
KOH2_H_M = 52.996

# Hourly UT climatological profile.
UT_HOURS = np.arange(
    0.0,
    24.0,
    1.0,
    dtype=float,
)

# Vertical integration grid.
ALT_MIN_KM = 90.0
ALT_MAX_KM = 2000.0
ALT_STEP_KM = 5.0

ALT_KM = np.arange(
    ALT_MIN_KM,
    ALT_MAX_KM + 0.5 * ALT_STEP_KM,
    ALT_STEP_KM,
    dtype=float,
)

# 0 = CCIR, 1 = URSI
CCIR_OR_URSI = 1
COEFFICIENT_LABEL = "URSI"

OUTPUT_ROOT = (
    BASE_ROOT
    / "TEC_REFERENCE_PYIRI_2019_2026"
)

HOURLY_ROOT = (
    OUTPUT_ROOT
    / "HOURLY"
)

INDICES_ROOT = (
    OUTPUT_ROOT
    / "_INDICES"
)

# Current daily IRI indices source linked from the IRI project.
IRI_APF107_URL = (
    "https://chain-new.chain-project.net/"
    "echaim_downloads/apf107.dat"
)

# NASA OMNI daily averages fallback.
NASA_OMNI_DAILY_URL = (
    "https://spdf.gsfc.nasa.gov/pub/data/omni/"
    "low_res_omni/omni_01_av.dat"
)

REQUEST_TIMEOUT = (
    30,
    90,
)


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate PyIRI/IRI climatological VTEC for available KOH2 "
            "observation days."
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
            "<data-root>/TEC_REFERENCE_PYIRI_2019_2026"
        ),
    )
    parser.add_argument(
        "--indices-root",
        type=Path,
        default=None,
        help=(
            "Solar-index cache directory. Default: <output-root>/_INDICES. "
            "For wrapper validation, point this to the existing operational "
            "PyIRI _INDICES directory."
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
    global OUTPUT_ROOT, HOURLY_ROOT, INDICES_ROOT

    BASE_ROOT = args.data_root.resolve()

    OUTPUT_ROOT = (
        args.output_root.resolve()
        if args.output_root is not None
        else BASE_ROOT / "TEC_REFERENCE_PYIRI_2019_2026"
    )

    HOURLY_ROOT = OUTPUT_ROOT / "HOURLY"

    INDICES_ROOT = (
        args.indices_root.resolve()
        if args.indices_root is not None
        else OUTPUT_ROOT / "_INDICES"
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
# LOAD PYIRI
# =============================================================================

def load_pyiri():
    try:
        import PyIRI
        import PyIRI.main_library as iri
    except ImportError as exc:
        raise RuntimeError(
            "\nPyIRI is not installed in this Python environment.\n\n"
            "Activate pytecgg_env and run:\n\n"
            "    python -m pip install PyIRI requests\n"
        ) from exc

    try:
        version = importlib.metadata.version(
            "PyIRI"
        )
    except Exception:
        version = "unknown"

    print(
        "PyIRI version:",
        version,
    )

    return (
        PyIRI,
        iri,
        version,
    )


PYIRI, IRI, PYIRI_VERSION = load_pyiri()


# =============================================================================
# DATE / OBSERVATION-DAY DISCOVERY
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


def extract_doy(
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


def discover_observation_days(
    year: int,
):
    root = (
        BASE_ROOT
        / str(
            year
        )
    )

    days = set()

    if not root.is_dir():
        return []

    for path in root.glob(
        rf"*\*\PyTECGg_OUTPUT\{STATION}_*_{year}_PyTECGg_VEQ.parquet"
    ):
        doy = extract_doy(
            path,
            year,
        )

        if doy is not None:
            days.add(
                doy
            )

    for path in root.glob(
        rf"*\*\pyOASIS_OUTPUT\INDICES\TEC\{STATION}_*_{year}_L1L2.TEC"
    ):
        doy = extract_doy(
            path,
            year,
        )

        if doy is not None:
            days.add(
                doy
            )

    return sorted(
        days
    )


# =============================================================================
# SOLAR INDEX DOWNLOAD
# =============================================================================

def download_file(
    url: str,
    path: Path,
):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    response = requests.get(
        url,
        timeout=REQUEST_TIMEOUT,
    )

    response.raise_for_status()

    content = response.content

    if len(
        content
    ) < 1000:
        raise RuntimeError(
            f"Downloaded file is unexpectedly small: {url}"
        )

    path.write_bytes(
        content
    )

    return path


def normalize_apf107_year(
    raw_year: int,
):
    """
    IRI apf107.dat uses an I3 year field in historical distributions.

    Common forms are:
        19   -> 2019-like 2-digit representation
        119  -> 2019 (years since 1900)
        2019 -> full year in newer/generated variants

    Handle all forms defensively.
    """
    if raw_year >= 1900:
        return raw_year

    if raw_year >= 100:
        return (
            1900
            + raw_year
        )

    # We only need 2019-2026 here.
    if raw_year <= 50:
        return (
            2000
            + raw_year
        )

    return (
        1900
        + raw_year
    )


def parse_apf107(
    path: Path,
):
    rows = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as f:
        for line in f:
            parts = line.split()

            # Expected fields:
            # y m d + 8 x ap3h + ap_daily + -11 + F107d + F107_81 + F107_365
            if len(
                parts
            ) < 16:
                continue

            try:
                raw_year = int(
                    parts[
                        0
                    ]
                )

                year = normalize_apf107_year(
                    raw_year
                )

                month = int(
                    parts[
                        1
                    ]
                )

                day = int(
                    parts[
                        2
                    ]
                )

                f107_daily = float(
                    parts[
                        -3
                    ]
                )

                f107_81 = float(
                    parts[
                        -2
                    ]
                )

                f107_365 = float(
                    parts[
                        -1
                    ]
                )

                ap_daily = float(
                    parts[
                        11
                    ]
                )

                date = pd.Timestamp(
                    year=year,
                    month=month,
                    day=day,
                    tz="UTC",
                )

            except Exception:
                continue

            if not (
                40.0
                <= f107_daily
                <= 500.0
            ):
                continue

            rows.append({
                "date":
                    date,
                "f107_daily":
                    f107_daily,
                "f107_81":
                    f107_81,
                "f107_365":
                    f107_365,
                "ap_daily":
                    ap_daily,
                "solar_index_source":
                    "IRI_apf107.dat",
            })

    if not rows:
        raise RuntimeError(
            f"No valid records parsed from {path}"
        )

    df = pd.DataFrame(
        rows
    )

    df = (
        df.drop_duplicates(
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

    return df


def parse_omni_daily(
    path: Path,
):
    """
    NASA OMNI2 daily average format.

    Relevant words:
        1  Year
        2  DOY
        3  Hour
        51 F10.7 index

    Daily average records use the same field structure as OMNI2 hourly files.
    """
    rows = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as f:
        for line in f:
            parts = line.split()

            if len(
                parts
            ) < 51:
                continue

            try:
                year = int(
                    parts[
                        0
                    ]
                )

                doy = int(
                    parts[
                        1
                    ]
                )

                f107 = float(
                    parts[
                        50
                    ]
                )

                date = pd.Timestamp(
                    date_from_year_doy(
                        year,
                        doy,
                    )
                )

            except Exception:
                continue

            if not (
                40.0
                <= f107
                <= 500.0
            ):
                continue

            rows.append({
                "date":
                    date,
                "f107_daily":
                    f107,
                "f107_81":
                    np.nan,
                "f107_365":
                    np.nan,
                "ap_daily":
                    np.nan,
                "solar_index_source":
                    "NASA_OMNI_daily_F10.7",
            })

    if not rows:
        raise RuntimeError(
            f"No valid OMNI F10.7 records parsed from {path}"
        )

    df = pd.DataFrame(
        rows
    )

    return (
        df.drop_duplicates(
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


def load_solar_indices():
    """
    Production choice for PyIRI:
    use NASA OMNI daily F10.7 directly.

    PyIRI IRI_density_1day() accepts one user-provided F10.7 value, so the
    additional apf107.dat fields are not required for this calculation.
    """
    omni_path = (
        INDICES_ROOT
        / "omni_01_av.dat"
    )

    if (
        not omni_path.is_file()
        or omni_path.stat().st_size < 1000
    ):
        print(
            "Downloading NASA OMNI daily F10.7 ..."
        )

        download_file(
            NASA_OMNI_DAILY_URL,
            omni_path,
        )
    else:
        print(
            "Using cached NASA OMNI daily F10.7:",
            omni_path,
        )

    df = parse_omni_daily(
        omni_path
    )

    print(
        "Solar index source: NASA OMNI daily F10.7"
    )

    print(
        "Index date range:",
        df[
            "date"
        ].min(),
        "to",
        df[
            "date"
        ].max(),
    )

    return df


def solar_row_for_day(
    solar: pd.DataFrame,
    year: int,
    doy: int,
):
    target = pd.Timestamp(
        date_from_year_doy(
            year,
            doy,
        )
    )

    found = solar[
        solar[
            "date"
        ]
        == target
    ]

    if found.empty:
        raise RuntimeError(
            f"No F10.7 value found for {year} DOY {doy:03d}"
        )

    return found.iloc[
        0
    ]


# =============================================================================
# IRI CALCULATION
# =============================================================================

def run_pyiri_day(
    year: int,
    doy: int,
    f107: float,
):
    date = date_from_year_doy(
        year,
        doy,
    )

    alon = np.array(
        [
            KOH2_LON,
        ],
        dtype=float,
    )

    alat = np.array(
        [
            KOH2_LAT,
        ],
        dtype=float,
    )

    t0 = time.perf_counter()

    (
        f2,
        f1,
        epeak,
        es_peak,
        sun,
        mag,
        edp,
    ) = IRI.IRI_density_1day(
        date.year,
        date.month,
        date.day,
        UT_HOURS,
        alon,
        alat,
        ALT_KM,
        float(
            f107
        ),
        PYIRI.coeff_dir,
        CCIR_OR_URSI,
    )

    # PyIRI returns EDP with dimensions:
    # [N_time, N_altitude, N_horizontal]
    vtec = IRI.edp_to_vtec(
        edp,
        ALT_KM,
        min_alt=ALT_MIN_KM,
        max_alt=ALT_MAX_KM,
    )

    vtec = np.asarray(
        vtec,
        dtype=float,
    ).reshape(
        -1
    )

    elapsed = (
        time.perf_counter()
        - t0
    )

    epochs = [
        pd.Timestamp(
            year=date.year,
            month=date.month,
            day=date.day,
            hour=int(
                ut
            ),
            tz="UTC",
        )
        for ut in UT_HOURS
    ]

    out = pd.DataFrame({
        "epoch":
            epochs,
        "ut_hour":
            UT_HOURS,
        "lat_deg":
            KOH2_LAT,
        "lon_deg":
            KOH2_LON,
        "f107_sfu":
            float(
                f107
            ),
        "iri_vtec_tecu":
            vtec,
    })

    return (
        out,
        elapsed,
    )


# =============================================================================
# DAILY / MONTHLY / YEARLY SUMMARIES
# =============================================================================

def daily_statistics(
    year: int,
    doy: int,
    hourly: pd.DataFrame,
    solar_row,
    elapsed_s: float,
):
    x = pd.to_numeric(
        hourly[
            "iri_vtec_tecu"
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    x = x[
        np.isfinite(
            x
        )
    ]

    date = date_from_year_doy(
        year,
        doy,
    )

    return {
        "date":
            date.date().isoformat(),
        "year":
            year,
        "month":
            date.month,
        "doy":
            doy,
        "station":
            STATION,
        "lat_deg":
            KOH2_LAT,
        "lon_deg":
            KOH2_LON,
        "implementation":
            "PyIRI",
        "pyiri_version":
            PYIRI_VERSION,
        "coefficient_option":
            COEFFICIENT_LABEL,
        "time_step_hours":
            1.0,
        "tec_alt_min_km":
            ALT_MIN_KM,
        "tec_alt_max_km":
            ALT_MAX_KM,
        "tec_alt_step_km":
            ALT_STEP_KM,
        "solar_index_source":
            solar_row[
                "solar_index_source"
            ],
        "f107_daily_sfu":
            float(
                solar_row[
                    "f107_daily"
                ]
            ),
        "f107_81_sfu":
            float(
                solar_row[
                    "f107_81"
                ]
            )
            if pd.notna(
                solar_row[
                    "f107_81"
                ]
            )
            else np.nan,
        "f107_365_sfu":
            float(
                solar_row[
                    "f107_365"
                ]
            )
            if pd.notna(
                solar_row[
                    "f107_365"
                ]
            )
            else np.nan,
        "ap_daily":
            float(
                solar_row[
                    "ap_daily"
                ]
            )
            if pd.notna(
                solar_row[
                    "ap_daily"
                ]
            )
            else np.nan,
        "n_hourly":
            int(
                len(
                    x
                )
            ),
        "iri_vtec_mean_tecu":
            float(
                np.mean(
                    x
                )
            )
            if len(
                x
            )
            else np.nan,
        "iri_vtec_median_tecu":
            float(
                np.median(
                    x
                )
            )
            if len(
                x
            )
            else np.nan,
        "iri_vtec_min_tecu":
            float(
                np.min(
                    x
                )
            )
            if len(
                x
            )
            else np.nan,
        "iri_vtec_max_tecu":
            float(
                np.max(
                    x
                )
            )
            if len(
                x
            )
            else np.nan,
        "iri_vtec_std_tecu":
            float(
                np.std(
                    x,
                    ddof=0,
                )
            )
            if len(
                x
            )
            else np.nan,
        "runtime_seconds":
            elapsed_s,
    }


def equal_day_yearly_summary(
    daily: pd.DataFrame,
):
    if daily.empty:
        return pd.DataFrame()

    rows = []

    for year, g in daily.groupby(
        "year"
    ):
        rows.append({
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
            "mean_daily_f107_sfu":
                g[
                    "f107_daily_sfu"
                ].mean(),
            "median_daily_f107_sfu":
                g[
                    "f107_daily_sfu"
                ].median(),
            "mean_daily_iri_vtec_tecu":
                g[
                    "iri_vtec_mean_tecu"
                ].mean(),
            "median_of_daily_iri_median_tecu":
                g[
                    "iri_vtec_median_tecu"
                ].median(),
            "mean_daily_iri_max_tecu":
                g[
                    "iri_vtec_max_tecu"
                ].mean(),
            "median_daily_iri_std_tecu":
                g[
                    "iri_vtec_std_tecu"
                ].median(),
        })

    return pd.DataFrame(
        rows
    ).sort_values(
        "year"
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
    ), g in daily.groupby(
        [
            "year",
            "month",
        ]
    ):
        rows.append({
            "year":
                int(
                    year
                ),
            "month":
                int(
                    month
                ),
            "n_days":
                int(
                    len(
                        g
                    )
                ),
            "mean_f107_sfu":
                g[
                    "f107_daily_sfu"
                ].mean(),
            "median_f107_sfu":
                g[
                    "f107_daily_sfu"
                ].median(),
            "mean_iri_vtec_tecu":
                g[
                    "iri_vtec_mean_tecu"
                ].mean(),
            "median_iri_vtec_tecu":
                g[
                    "iri_vtec_median_tecu"
                ].median(),
            "mean_daily_max_iri_vtec_tecu":
                g[
                    "iri_vtec_max_tecu"
                ].mean(),
        })

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "month",
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

    HOURLY_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    INDICES_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    pd.set_option(
        "display.max_columns",
        None,
    )

    pd.set_option(
        "display.width",
        240,
    )

    print(
        "=" * 112
    )

    print(
        "KOH2 PyIRI / IRI CLIMATOLOGICAL VTEC"
    )

    print(
        "=" * 112
    )

    print(
        "Coordinates:",
        f"{KOH2_LAT:.8f}, {KOH2_LON:.8f}",
    )

    print(
        "Coefficient option:",
        COEFFICIENT_LABEL,
    )

    print(
        "UT samples:",
        len(
            UT_HOURS
        ),
        "(hourly)",
    )

    print(
        "TEC integration:",
        f"{ALT_MIN_KM:.0f}-{ALT_MAX_KM:.0f} km, "
        f"{ALT_STEP_KM:.1f}-km vertical step",
    )

    print(
        "Output:",
        OUTPUT_ROOT,
    )

    print(
        "Solar-index cache:",
        INDICES_ROOT,
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
        "=" * 112
    )

    solar = load_solar_indices()

    daily_rows = []
    availability_rows = []

    years = YEARS

    for year in years:
        days = discover_observation_days(
            year
        )

        if SELECTED_DOY is not None:
            days = [
                doy
                for doy in days
                if doy == SELECTED_DOY
            ]

        print()
        print(
            "=" * 112
        )

        print(
            f"{year}: days to model = {len(days)}"
        )

        if not days:
            print(
                "    No KOH2 production days found for this year."
            )

        print(
            "=" * 112
        )

        for i, doy in enumerate(
            days,
            1,
        ):
            print(
                f"{i:3d}/{len(days):3d}  "
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
                "status":
                    "FAILED",
            }

            try:
                srow = solar_row_for_day(
                    solar,
                    year,
                    doy,
                )

                f107 = float(
                    srow[
                        "f107_daily"
                    ]
                )

                print(
                    f"    F10.7 = {f107:.1f} sfu "
                    f"({srow['solar_index_source']})"
                )

                hourly, elapsed = run_pyiri_day(
                    year,
                    doy,
                    f107,
                )

                hourly_file = (
                    HOURLY_ROOT
                    / str(
                        year
                    )
                    / (
                        f"{STATION}_{year}_{doy:03d}_"
                        f"PyIRI_{COEFFICIENT_LABEL}_hourly.csv"
                    )
                )

                hourly_file.parent.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                hourly.to_csv(
                    hourly_file,
                    index=False,
                )

                row = daily_statistics(
                    year,
                    doy,
                    hourly,
                    srow,
                    elapsed,
                )

                daily_rows.append(
                    row
                )

                status.update({
                    "status":
                        "OK",
                    "f107_sfu":
                        f107,
                    "n_hourly":
                        row[
                            "n_hourly"
                        ],
                    "iri_vtec_mean_tecu":
                        row[
                            "iri_vtec_mean_tecu"
                        ],
                    "iri_vtec_median_tecu":
                        row[
                            "iri_vtec_median_tecu"
                        ],
                    "iri_vtec_min_tecu":
                        row[
                            "iri_vtec_min_tecu"
                        ],
                    "iri_vtec_max_tecu":
                        row[
                            "iri_vtec_max_tecu"
                        ],
                    "runtime_seconds":
                        elapsed,
                    "hourly_file":
                        str(
                            hourly_file
                        ),
                })

                print(
                    f"    mean   = {row['iri_vtec_mean_tecu']:.3f} TECU"
                )

                print(
                    f"    median = {row['iri_vtec_median_tecu']:.3f} TECU"
                )

                print(
                    f"    range  = {row['iri_vtec_min_tecu']:.3f} .. "
                    f"{row['iri_vtec_max_tecu']:.3f} TECU"
                )

                print(
                    f"    runtime= {elapsed:.2f} s"
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

            # Restart-friendly snapshots.
            pd.DataFrame(
                daily_rows
            ).to_csv(
                OUTPUT_ROOT
                / (
                    f"{STATION}_2019_2026_"
                    "PyIRI_daily_statistics.csv"
                ),
                index=False,
            )

            pd.DataFrame(
                availability_rows
            ).to_csv(
                OUTPUT_ROOT
                / (
                    f"{STATION}_2019_2026_"
                    "PyIRI_availability.csv"
                ),
                index=False,
            )

    daily = pd.DataFrame(
        daily_rows
    )

    availability = pd.DataFrame(
        availability_rows
    )

    yearly = equal_day_yearly_summary(
        daily
    )

    monthly = monthly_summary(
        daily
    )

    daily_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "PyIRI_daily_statistics.csv"
        )
    )

    availability_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "PyIRI_availability.csv"
        )
    )

    yearly_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "PyIRI_yearly_equal_day_summary.csv"
        )
    )

    monthly_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "PyIRI_monthly_summary.csv"
        )
    )

    daily.to_csv(
        daily_file,
        index=False,
    )

    availability.to_csv(
        availability_file,
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

    print()
    print(
        "=" * 112
    )

    print(
        "PYIRI RUN COMPLETE"
    )

    print(
        "=" * 112
    )

    print(
        "Daily statistics:",
        daily_file,
    )

    print(
        "Availability:",
        availability_file,
    )

    print(
        "Yearly equal-day:",
        yearly_file,
    )

    print(
        "Monthly:",
        monthly_file,
    )

    if not yearly.empty:
        print()
        print(
            "YEARLY EQUAL-DAY SUMMARY"
        )

        print(
            "-" * 112
        )

        print(
            yearly.to_string(
                index=False
            )
        )

        print(
            "-" * 112
        )


if __name__ == "__main__":
    main()
