from __future__ import annotations

r"""
KOH2 2019-2026 comparison with MIT Haystack / CEDAR Madrigal GNSS VTEC
=====================================================================

PURPOSE
-------
Compare the already-produced KOH2 TEC results from:

    1. PyTECGg VEq at the KOH2 station
    2. PyTECGg VTEC at PyTECGg IPPs
    3. pyOASIS VTEC at pyOASIS IPPs

against the MIT Haystack MAPGPS / Madrigal gridded GNSS VTEC product:

    Instrument: 8000  World-wide GNSS Receiver Network
    Kindat:     3500  TEC binned 1 degree by 1 degree by 5 min

The script also checks kindat 3506 (daily receiver-site list), when available,
to determine whether KOH2 itself contributed observations to the Madrigal map.

SCIENTIFIC INTERPRETATION
-------------------------
Madrigal/MAPGPS is an observation-derived GNSS product with much finer nominal
sampling than the IGS combined GIM.  However, it is still derived from a global
GNSS receiver network.  If KOH2 is found in the daily site list, that day's
comparison must NOT be described as fully independent validation.

MATCHING METHOD
---------------
Madrigal is not treated as a gap-free analytic GIM.

For each KOH2-derived VTEC/VEq point:
  * choose the nearest Madrigal 5-minute map epoch;
  * require |dt| <= 180 s;
  * choose the nearest AVAILABLE Madrigal grid cell;
  * require horizontal separation <= 80 km.

This strict nearest-bin method avoids interpolating across missing Madrigal
coverage in Antarctica.

Residual convention:
    KOH2 method TEC - Madrigal VTEC

IMPORTANT
---------
This production version uses TEST_ONLY=False and processes all available
KOH2 days from 2019 through 2026. It downloads each native Madrigal daily HDF5
once and performs the Antarctic geographic filtering locally.

It also records the accepted Madrigal match fraction for every daily
comparison so that Antarctic spatial-coverage limitations can be quantified.

The script reuses the already validated readers/statistics from:
    validate_tec_igs_2019_2026_E_drive_V4.py

Place both scripts in the same directory.

OUTPUT ROOT
-----------
    E:\KOH2data\TEC_VALIDATION_MADRIGAL_2019_2026

MADRIGAL DATA ACCESS
--------------------
The Madrigal service requires a user's name, email, and affiliation for data
access logging.  This script asks for them interactively at startup and does
not save them in the script.

Required packages in pytecgg_env:
    python -m pip install madrigalweb h5py
"""

from pathlib import Path
from datetime import datetime, timedelta, timezone
import argparse
import importlib.util
import math
import sys
import time

import h5py
import numpy as np
import pandas as pd


# =============================================================================
# SETTINGS
# =============================================================================

BASE_ROOT = Path(".")
YEARS = list(range(2019, 2027))
STATION = "KOH2"
SELECTED_DOY = None

MADRIGAL_URL = "https://cedar.openmadrigal.org"

MADRIGAL_INSTRUMENT = 8000
MADRIGAL_VTEC_KINDAT = 3500
MADRIGAL_SITE_KINDAT = 3506

MIN_ELEVATION_DEG = 30.0

# Madrigal maps are nominally 5 min × 1 degree × 1 degree.
MAX_TIME_DIFFERENCE_S = 180.0
MAX_SPATIAL_DISTANCE_KM = 80.0

# Add a margin around the actual KOH2 IPP footprint before requesting Madrigal.
# This controls transfer size, not the final nearest-cell acceptance radius.
REGION_MARGIN_DEG = 2.0

# Safe geographic fallback around KOH2 if a method has no usable IPPs.
FALLBACK_LAT_MIN = -75.0
FALLBACK_LAT_MAX = -50.0
FALLBACK_LON_MIN = -85.0
FALLBACK_LON_MAX = -35.0

MASTER_OUTPUT = (
    BASE_ROOT
    / "TEC_VALIDATION_MADRIGAL_2019_2026"
)

CACHE_ROOT = (
    MASTER_OUTPUT
    / "_MADRIGAL_CACHE"
)

BY_YEAR_OUTPUT = (
    MASTER_OUTPUT
    / "BY_YEAR"
)

SAVE_MATCHED_POINTS = False

# FAST MODE:
# Download the native daily Madrigal HDF5 once, then filter it locally.
# This avoids the slow server-side isprint geographic scan.
USE_NATIVE_HDF5_LOCAL_FILTER = True
LOCAL_HDF5_CHUNK_ROWS = 1_000_000


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare KOH2 PyTECGg and pyOASIS TEC products with "
            "MIT Haystack / CEDAR Madrigal MAPGPS VTEC."
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
            "<data-root>/TEC_VALIDATION_MADRIGAL_2019_2026"
        ),
    )
    parser.add_argument(
        "--madrigal-cache-root",
        type=Path,
        default=None,
        help=(
            "Optional existing Madrigal cache root. Default: "
            "<output-root>/_MADRIGAL_CACHE"
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
    global MASTER_OUTPUT, CACHE_ROOT, BY_YEAR_OUTPUT

    BASE_ROOT = args.data_root.resolve()

    MASTER_OUTPUT = (
        args.output_root.resolve()
        if args.output_root is not None
        else BASE_ROOT / "TEC_VALIDATION_MADRIGAL_2019_2026"
    )

    CACHE_ROOT = (
        args.madrigal_cache_root.resolve()
        if args.madrigal_cache_root is not None
        else MASTER_OUTPUT / "_MADRIGAL_CACHE"
    )

    BY_YEAR_OUTPUT = MASTER_OUTPUT / "BY_YEAR"

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
# IMPORT THE VALIDATED PUBLICATION IGS CORE
# =============================================================================

def load_validation_core():
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
            "Place this Madrigal script in the same folder as the validated "
            "publication IGS script."
        )

    spec = importlib.util.spec_from_file_location(
        "koh2_validation_v4_core",
        path,
    )

    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"Could not import validation core: {path}"
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    print(
        "Validation core:",
        path,
    )

    return module


CORE = load_validation_core()

KOH2_LAT = CORE.KOH2_LAT
KOH2_LON = CORE.KOH2_LON
KOH2_H = CORE.KOH2_H

load_pytecgg = CORE.load_pytecgg
load_pyoasis = CORE.load_pyoasis
calculate_statistics = CORE.calculate_statistics
GlobalAccumulator = CORE.GlobalAccumulator


# =============================================================================
# MADRIGAL API
# =============================================================================

def load_madrigal_api():
    try:
        import madrigalWeb.madrigalWeb as madrigal_api
    except ImportError as exc:
        raise RuntimeError(
            "\nThe Madrigal Python API is not installed in this environment.\n"
            "Activate pytecgg_env and run:\n\n"
            "    python -m pip install madrigalweb h5py\n"
        ) from exc

    return madrigal_api


MADRIGAL_API = load_madrigal_api()


def get_user_identity():
    print()
    print(
        "Madrigal requires user identification for data-access logging."
    )

    fullname = input(
        "Your full name       : "
    ).strip()

    email = input(
        "Your email           : "
    ).strip()

    affiliation = input(
        "Your affiliation     : "
    ).strip()

    if not fullname or not email or not affiliation:
        raise RuntimeError(
            "Name, email, and affiliation are required by Madrigal."
        )

    return (
        fullname,
        email,
        affiliation,
    )


# =============================================================================
# KOH2 FILE DISCOVERY
# =============================================================================

def extract_doy(
    path: Path,
    year: int,
):
    import re

    m = re.search(
        r"_(\d{3})_(\d{4})_",
        path.name,
    )

    if not m:
        return None

    doy = int(
        m.group(1)
    )

    found_year = int(
        m.group(2)
    )

    if found_year != year:
        return None

    return doy


def discover_year_files(
    year: int,
):
    root = (
        BASE_ROOT
        / str(year)
    )

    pytecgg = {}
    pyoasis = {}

    if not root.is_dir():
        return (
            pytecgg,
            pyoasis,
            [],
        )

    for path in root.glob(
        rf"*\*\PyTECGg_OUTPUT\{STATION}_*_{year}_PyTECGg_VEQ.parquet"
    ):
        doy = extract_doy(
            path,
            year,
        )

        if doy is not None:
            pytecgg[
                doy
            ] = path

    for path in root.glob(
        rf"*\*\pyOASIS_OUTPUT\INDICES\TEC\{STATION}_*_{year}_L1L2.TEC"
    ):
        doy = extract_doy(
            path,
            year,
        )

        if doy is not None:
            pyoasis[
                doy
            ] = path

    days = sorted(
        set(
            pytecgg
        )
        | set(
            pyoasis
        )
    )

    return (
        pytecgg,
        pyoasis,
        days,
    )


def datetime_from_year_doy(
    year: int,
    doy: int,
):
    return datetime.strptime(
        f"{year}-{doy:03d}",
        "%Y-%j",
    ).replace(
        tzinfo=timezone.utc
    )


# =============================================================================
# DETERMINE THE REGIONAL REQUEST BOX
# =============================================================================

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


def get_daily_region(
    pt: pd.DataFrame | None,
    po: pd.DataFrame | None,
):
    lats = [
        KOH2_LAT,
    ]

    lons = [
        KOH2_LON,
    ]

    if pt is not None:
        g = pt[
            pd.to_numeric(
                pt[
                    "ele"
                ],
                errors="coerce",
            )
            >= MIN_ELEVATION_DEG
        ]

        lat = pd.to_numeric(
            g[
                "lat_ipp"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        lon = pd.to_numeric(
            g[
                "lon_ipp"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        good = (
            np.isfinite(
                lat
            )
            & np.isfinite(
                lon
            )
        )

        if np.any(
            good
        ):
            lats.extend(
                lat[
                    good
                ].tolist()
            )

            lons.extend(
                wrap_lon(
                    lon[
                        good
                    ]
                ).tolist()
            )

    if po is not None:
        g = po[
            pd.to_numeric(
                po[
                    "elevation"
                ],
                errors="coerce",
            )
            >= MIN_ELEVATION_DEG
        ]

        lat = pd.to_numeric(
            g[
                "lat_ipp"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        lon = pd.to_numeric(
            g[
                "lon_ipp"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        good = (
            np.isfinite(
                lat
            )
            & np.isfinite(
                lon
            )
        )

        if np.any(
            good
        ):
            lats.extend(
                lat[
                    good
                ].tolist()
            )

            lons.extend(
                wrap_lon(
                    lon[
                        good
                    ]
                ).tolist()
            )

    if len(
        lats
    ) <= 1:
        return (
            FALLBACK_LAT_MIN,
            FALLBACK_LAT_MAX,
            FALLBACK_LON_MIN,
            FALLBACK_LON_MAX,
        )

    lat_min = max(
        -90.0,
        float(
            np.nanmin(
                lats
            )
        )
        - REGION_MARGIN_DEG,
    )

    lat_max = min(
        90.0,
        float(
            np.nanmax(
                lats
            )
        )
        + REGION_MARGIN_DEG,
    )

    lon_arr = wrap_lon(
        lons
    )

    lon_min = max(
        -180.0,
        float(
            np.nanmin(
                lon_arr
            )
        )
        - REGION_MARGIN_DEG,
    )

    lon_max = min(
        180.0,
        float(
            np.nanmax(
                lon_arr
            )
        )
        + REGION_MARGIN_DEG,
    )

    return (
        lat_min,
        lat_max,
        lon_min,
        lon_max,
    )


# =============================================================================
# MADRIGAL EXPERIMENT / FILE DISCOVERY
# =============================================================================

def experiment_start_dt(
    exp,
):
    return datetime(
        int(exp.startyear),
        int(exp.startmonth),
        int(exp.startday),
        int(exp.starthour),
        int(exp.startmin),
        int(exp.startsec),
        tzinfo=timezone.utc,
    )


def experiment_end_dt(
    exp,
):
    return datetime(
        int(exp.endyear),
        int(exp.endmonth),
        int(exp.endday),
        int(exp.endhour),
        int(exp.endmin),
        int(exp.endsec),
        tzinfo=timezone.utc,
    )


def choose_experiment_for_target_day(
    exps,
    target_start,
    target_end,
):
    """
    Madrigal getExperiments() returns experiments whose time spans overlap the
    requested interval.  Around midnight/year boundaries this can include the
    previous day's experiment.  Prefer an experiment that STARTS on the target
    UTC date, then maximize temporal overlap with the requested day.
    """
    if not exps:
        return None

    tec_exps = [
        exp
        for exp in exps
        if "tec" in str(
            getattr(
                exp,
                "name",
                "",
            )
        ).lower()
    ]

    candidates = (
        tec_exps
        if tec_exps
        else list(
            exps
        )
    )

    scored = []

    for exp in candidates:
        try:
            s = experiment_start_dt(
                exp
            )

            e = experiment_end_dt(
                exp
            )
        except Exception:
            continue

        overlap_start = max(
            s,
            target_start,
        )

        overlap_end = min(
            e,
            target_end,
        )

        overlap_s = max(
            0.0,
            (
                overlap_end
                - overlap_start
            ).total_seconds(),
        )

        exact_start_date = int(
            s.date()
            == target_start.date()
        )

        contains_midday = int(
            s
            <= (
                target_start
                + timedelta(
                    hours=12
                )
            )
            <= e
        )

        # exact target-date start is the strongest discriminator;
        # overlap and midday containment break ties.
        score = (
            exact_start_date,
            contains_midday,
            overlap_s,
            s,
        )

        scored.append(
            (
                score,
                exp,
            )
        )

    if not scored:
        return None

    scored.sort(
        key=lambda x:
            x[
                0
            ],
        reverse=True,
    )

    return scored[
        0
    ][
        1
    ]


def experiment_for_day(
    client,
    year: int,
    doy: int,
):
    start = datetime_from_year_doy(
        year,
        doy,
    )

    end = (
        start
        + timedelta(
            days=1
        )
    )

    exps = client.getExperiments(
        MADRIGAL_INSTRUMENT,
        start.year,
        start.month,
        start.day,
        0,
        0,
        0,
        end.year,
        end.month,
        end.day,
        0,
        0,
        0,
        local=0,
    )

    if not exps:
        raise RuntimeError(
            f"No Madrigal experiment found for {year} DOY {doy:03d}"
        )

    exp = choose_experiment_for_target_day(
        exps,
        start,
        end,
    )

    if exp is None:
        raise RuntimeError(
            f"No suitable Madrigal TEC experiment found for "
            f"{year} DOY {doy:03d}"
        )

    exp_id = int(
        getattr(
            exp,
            "id",
            -1,
        )
    )

    if exp_id != -1:
        return (
            client,
            exp,
        )

    # OpenMadrigal can return a remote pointer with id=-1.
    remote_url = getattr(
        exp,
        "madrigalUrl",
        None,
    )

    if not remote_url:
        raise RuntimeError(
            "Madrigal returned a remote experiment without madrigalUrl."
        )

    remote_client = MADRIGAL_API.MadrigalData(
        remote_url
    )

    local_exps = remote_client.getExperiments(
        MADRIGAL_INSTRUMENT,
        start.year,
        start.month,
        start.day,
        0,
        0,
        0,
        end.year,
        end.month,
        end.day,
        0,
        0,
        0,
        local=1,
    )

    if not local_exps:
        raise RuntimeError(
            f"Remote Madrigal site had no local TEC experiment for "
            f"{year} DOY {doy:03d}"
        )

    local_exp = choose_experiment_for_target_day(
        local_exps,
        start,
        end,
    )

    if local_exp is None:
        raise RuntimeError(
            f"Remote Madrigal site had no suitable exact-day TEC experiment "
            f"for {year} DOY {doy:03d}"
        )

    return (
        remote_client,
        local_exp,
    )


def select_kindat_file(
    client,
    experiment,
    kindat: int,
):
    files = client.getExperimentFiles(
        int(
            experiment.id
        )
    )

    candidates = [
        f
        for f in files
        if int(
            getattr(
                f,
                "kindat",
                -9999,
            )
        )
        == kindat
    ]

    if not candidates:
        return None

    # Prefer category 1 (default/final) when available.
    defaults = [
        f
        for f in candidates
        if int(
            getattr(
                f,
                "category",
                -1,
            )
        )
        == 1
    ]

    if defaults:
        candidates = defaults

    # Prefer explicit final status if multiple remain.
    finals = [
        f
        for f in candidates
        if "final" in str(
            getattr(
                f,
                "status",
                "",
            )
        ).lower()
    ]

    if finals:
        candidates = finals

    return candidates[
        0
    ]


# =============================================================================
# FAST NATIVE HDF5 DOWNLOAD + LOCAL REGIONAL FILTER
# =============================================================================

def native_vtec_cache_file(
    year: int,
    doy: int,
):
    return (
        CACHE_ROOT
        / "VTEC_NATIVE_HDF5"
        / str(
            year
        )
        / (
            f"madrigal_native_vtec_{year}_{doy:03d}.hdf5"
        )
    )


def local_regional_cache_file(
    year: int,
    doy: int,
):
    return (
        CACHE_ROOT
        / "VTEC_REGIONAL_LOCAL"
        / str(
            year
        )
        / (
            f"madrigal_regional_local_{year}_{doy:03d}.parquet"
        )
    )


def download_native_vtec_hdf5(
    client,
    madrigal_file,
    year: int,
    doy: int,
    fullname: str,
    email: str,
    affiliation: str,
):
    """
    Download the stored daily HDF5 file without an isprint geographic scan.
    The local file is cached permanently for restart/rerun speed.
    """
    path = native_vtec_cache_file(
        year,
        doy,
    )

    if (
        path.is_file()
        and path.stat().st_size > 10_000
    ):
        print(
            "    Native Madrigal HDF5 cache:",
            path,
        )

        print(
            f"      size = {path.stat().st_size / (1024 ** 2):.1f} MB"
        )

        return path

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = path.with_suffix(
        ".hdf5.part"
    )

    if tmp.exists():
        tmp.unlink()

    print(
        "    Downloading native daily Madrigal HDF5 "
        "(no server-side spatial scan)..."
    )

    t0 = time.perf_counter()

    try:
        client.downloadFile(
            madrigal_file.name,
            str(
                tmp
            ),
            fullname,
            email,
            affiliation,
            format="hdf5",
        )
    except TypeError:
        # Compatibility with older madrigalWeb releases that expect
        # the final argument positionally.
        client.downloadFile(
            madrigal_file.name,
            str(
                tmp
            ),
            fullname,
            email,
            affiliation,
            "hdf5",
        )

    if (
        not tmp.is_file()
        or tmp.stat().st_size < 10_000
    ):
        raise RuntimeError(
            f"Native Madrigal HDF5 download failed or is too small: {tmp}"
        )

    tmp.replace(
        path
    )

    elapsed = time.perf_counter() - t0

    print(
        f"      downloaded {path.stat().st_size / (1024 ** 2):.1f} MB "
        f"in {elapsed:.1f} s"
    )

    return path


def _field_name(
    available,
    *candidates,
):
    lookup = {
        str(
            name
        ).lower():
            name
        for name in available
    }

    for candidate in candidates:
        found = lookup.get(
            candidate.lower()
        )

        if found is not None:
            return found

    return None


def load_native_regional_vtec(
    path: Path,
    year: int,
    doy: int,
    region,
):
    """
    Read the Madrigal HDF5 Table Layout in chunks and retain only the
    Antarctic region needed by KOH2.  This work runs on the user's computer,
    not on the Madrigal server.
    """
    regional_cache = local_regional_cache_file(
        year,
        doy,
    )

    if (
        regional_cache.is_file()
        and regional_cache.stat().st_size > 1000
    ):
        print(
            "    Local regional cache:",
            regional_cache,
        )

        return pd.read_parquet(
            regional_cache
        )

    regional_cache.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    lat_min, lat_max, lon_min, lon_max = region

    target_start = datetime_from_year_doy(
        year,
        doy,
    )

    target_end = (
        target_start
        + timedelta(
            days=1
        )
    )

    start_unix = target_start.timestamp()
    end_unix = target_end.timestamp()

    print(
        "    Filtering native HDF5 locally:"
    )

    print(
        f"      lat {lat_min:.2f} .. {lat_max:.2f}"
    )

    print(
        f"      lon {lon_min:.2f} .. {lon_max:.2f}"
    )

    t0 = time.perf_counter()

    chunks = []

    with h5py.File(
        path,
        "r",
    ) as h5:
        if (
            "Data" not in h5
            or "Table Layout" not in h5[
                "Data"
            ]
        ):
            raise RuntimeError(
                "Madrigal HDF5 has no Data/Table Layout dataset."
            )

        table = h5[
            "Data"
        ][
            "Table Layout"
        ]

        available = list(
            table.dtype.names
            or []
        )

        lat_field = _field_name(
            available,
            "gdlat",
        )

        lon_field = _field_name(
            available,
            "glon",
        )

        tec_field = _field_name(
            available,
            "tec",
        )

        dtec_field = _field_name(
            available,
            "dtec",
        )

        unix_field = _field_name(
            available,
            "ut1_unix",
        )

        year_field = _field_name(
            available,
            "year",
        )

        month_field = _field_name(
            available,
            "month",
        )

        day_field = _field_name(
            available,
            "day",
        )

        hour_field = _field_name(
            available,
            "hour",
        )

        min_field = _field_name(
            available,
            "min",
            "minute",
        )

        sec_field = _field_name(
            available,
            "sec",
            "second",
        )

        if (
            lat_field is None
            or lon_field is None
            or tec_field is None
        ):
            raise RuntimeError(
                "Required Madrigal HDF5 fields were not found. "
                f"Available fields include: {available[:80]}"
            )

        if unix_field is None:
            required_time = [
                year_field,
                month_field,
                day_field,
                hour_field,
                min_field,
                sec_field,
            ]

            if any(
                field is None
                for field in required_time
            ):
                raise RuntimeError(
                    "Neither ut1_unix nor complete calendar-time fields "
                    "were found in Madrigal HDF5."
                )

        wanted = [
            lat_field,
            lon_field,
            tec_field,
        ]

        if dtec_field is not None:
            wanted.append(
                dtec_field
            )

        if unix_field is not None:
            wanted.append(
                unix_field
            )
        else:
            wanted.extend([
                year_field,
                month_field,
                day_field,
                hour_field,
                min_field,
                sec_field,
            ])

        # Preserve order while removing duplicates.
        wanted = list(
            dict.fromkeys(
                wanted
            )
        )

        projected = table.fields(
            wanted
        )

        nrows = len(
            table
        )

        print(
            f"      native Table Layout rows = {nrows:,}"
        )

        for start in range(
            0,
            nrows,
            LOCAL_HDF5_CHUNK_ROWS,
        ):
            stop = min(
                nrows,
                start
                + LOCAL_HDF5_CHUNK_ROWS,
            )

            block = projected[
                start:stop
            ]

            lat = np.asarray(
                block[
                    lat_field
                ],
                dtype=float,
            )

            lon = wrap_lon(
                np.asarray(
                    block[
                        lon_field
                    ],
                    dtype=float,
                )
            )

            tec = np.asarray(
                block[
                    tec_field
                ],
                dtype=float,
            )

            mask = (
                np.isfinite(
                    lat
                )
                & np.isfinite(
                    lon
                )
                & np.isfinite(
                    tec
                )
                & (
                    lat
                    >= lat_min
                )
                & (
                    lat
                    <= lat_max
                )
                & (
                    lon
                    >= lon_min
                )
                & (
                    lon
                    <= lon_max
                )
            )

            if unix_field is not None:
                unix = np.asarray(
                    block[
                        unix_field
                    ],
                    dtype=float,
                )

                mask &= (
                    np.isfinite(
                        unix
                    )
                    & (
                        unix
                        >= start_unix
                    )
                    & (
                        unix
                        < end_unix
                    )
                )
            else:
                y = np.asarray(
                    block[
                        year_field
                    ],
                    dtype=int,
                )

                mo = np.asarray(
                    block[
                        month_field
                    ],
                    dtype=int,
                )

                d = np.asarray(
                    block[
                        day_field
                    ],
                    dtype=int,
                )

                mask &= (
                    y
                    == target_start.year
                ) & (
                    mo
                    == target_start.month
                ) & (
                    d
                    == target_start.day
                )

            idx = np.flatnonzero(
                mask
            )

            if len(
                idx
            ) == 0:
                continue

            if unix_field is not None:
                epoch = pd.to_datetime(
                    unix[
                        idx
                    ],
                    unit="s",
                    utc=True,
                    errors="coerce",
                )
            else:
                epoch = pd.to_datetime(
                    {
                        "year":
                            np.asarray(
                                block[
                                    year_field
                                ],
                                dtype=int,
                            )[
                                idx
                            ],
                        "month":
                            np.asarray(
                                block[
                                    month_field
                                ],
                                dtype=int,
                            )[
                                idx
                            ],
                        "day":
                            np.asarray(
                                block[
                                    day_field
                                ],
                                dtype=int,
                            )[
                                idx
                            ],
                        "hour":
                            np.asarray(
                                block[
                                    hour_field
                                ],
                                dtype=int,
                            )[
                                idx
                            ],
                        "minute":
                            np.asarray(
                                block[
                                    min_field
                                ],
                                dtype=int,
                            )[
                                idx
                            ],
                        "second":
                            np.asarray(
                                block[
                                    sec_field
                                ],
                                dtype=float,
                            )[
                                idx
                            ],
                    },
                    utc=True,
                    errors="coerce",
                )

            if dtec_field is not None:
                dtec = np.asarray(
                    block[
                        dtec_field
                    ],
                    dtype=float,
                )[
                    idx
                ]
            else:
                dtec = np.full(
                    len(
                        idx
                    ),
                    np.nan,
                    dtype=float,
                )

            chunks.append(
                pd.DataFrame({
                    "epoch":
                        epoch,
                    "gdlat":
                        lat[
                            idx
                        ],
                    "glon":
                        lon[
                            idx
                        ],
                    "tec":
                        tec[
                            idx
                        ],
                    "dtec":
                        dtec,
                })
            )

        if not chunks:
            return pd.DataFrame(
                columns=[
                    "epoch",
                    "gdlat",
                    "glon",
                    "tec",
                    "dtec",
                ]
            )

    df = pd.concat(
        chunks,
        ignore_index=True,
    )

    df = df.dropna(
        subset=[
            "epoch",
            "gdlat",
            "glon",
            "tec",
        ]
    )

    # Collapse any duplicate grid bins conservatively.
    df = (
        df.groupby(
            [
                "epoch",
                "gdlat",
                "glon",
            ],
            as_index=False,
        )
        .agg({
            "tec":
                "median",
            "dtec":
                "median",
        })
        .sort_values(
            [
                "epoch",
                "gdlat",
                "glon",
            ]
        )
        .reset_index(
            drop=True
        )
    )

    df.to_parquet(
        regional_cache,
        index=False,
    )

    elapsed = time.perf_counter() - t0

    print(
        f"      local regional rows = {len(df):,}"
    )

    print(
        f"      local filtering time = {elapsed:.1f} s"
    )

    print(
        "      saved regional cache:",
        regional_cache,
    )

    return df


# =============================================================================
# MADRIGAL REGIONAL VTEC DOWNLOAD THROUGH ISPRINT
# =============================================================================

VTEC_COLUMNS = [
    "year",
    "month",
    "day",
    "hour",
    "minute",
    "second",
    "gdlat",
    "glon",
    "tec",
    "dtec",
]


def regional_cache_file(
    year: int,
    doy: int,
):
    return (
        CACHE_ROOT
        / "VTEC_REGIONAL"
        / str(
            year
        )
        / (
            f"madrigal_vtec_V2_{year}_{doy:03d}.txt"
        )
    )


def site_cache_file(
    year: int,
    doy: int,
):
    return (
        CACHE_ROOT
        / "SITE_LIST"
        / str(
            year
        )
        / (
            f"madrigal_sites_V2_{year}_{doy:03d}.txt"
        )
    )


def download_regional_vtec(
    client,
    madrigal_file,
    year: int,
    doy: int,
    region,
    fullname: str,
    email: str,
    affiliation: str,
):
    path = regional_cache_file(
        year,
        doy,
    )

    if (
        path.is_file()
        and path.stat().st_size > 100
    ):
        print(
            "    Madrigal regional cache:",
            path,
        )

        return path

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    lat_min, lat_max, lon_min, lon_max = region

    parms = (
        "year,month,day,hour,min,sec,"
        "gdlat,glon,tec,dtec"
    )

    filters = (
        f"filter=gdlat,{lat_min:.4f},{lat_max:.4f} "
        f"filter=glon,{lon_min:.4f},{lon_max:.4f}"
    )

    print(
        "    Requesting Madrigal regional VTEC:"
    )

    print(
        f"      lat {lat_min:.2f} .. {lat_max:.2f}"
    )

    print(
        f"      lon {lon_min:.2f} .. {lon_max:.2f}"
    )

    print(
        "      This server-side scan may take several minutes."
    )

    client.isprint(
        madrigal_file.name,
        parms,
        filters,
        fullname,
        email,
        affiliation,
        outputFile=str(
            path
        ),
    )

    if (
        not path.is_file()
        or path.stat().st_size < 20
    ):
        raise RuntimeError(
            f"Madrigal regional extraction was empty or invalid: {path}"
        )

    return path


def load_regional_vtec(
    path: Path,
    year: int | None = None,
    doy: int | None = None,
):
    try:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=VTEC_COLUMNS,
            engine="python",
            comment="#",
        )
    except Exception as exc:
        raise RuntimeError(
            f"Could not parse Madrigal regional file: {path}"
        ) from exc

    if len(
        df
    ) == 0:
        return df

    for col in VTEC_COLUMNS:
        df[
            col
        ] = pd.to_numeric(
            df[
                col
            ],
            errors="coerce",
        )

    # Madrigal's reported YEAR/MONTH/DAY/HOUR/MIN/SEC represent the
    # measurement time. Support fractional seconds defensively.
    base = pd.to_datetime(
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
            "hour":
                df[
                    "hour"
                ],
            "minute":
                df[
                    "minute"
                ],
        },
        utc=True,
        errors="coerce",
    )

    df[
        "epoch"
    ] = (
        base
        + pd.to_timedelta(
            df[
                "second"
            ].fillna(
                0.0
            ),
            unit="s",
        )
    )

    df[
        "glon"
    ] = wrap_lon(
        df[
            "glon"
        ].to_numpy(
            dtype=float
        )
    )

    df = df.dropna(
        subset=[
            "epoch",
            "gdlat",
            "glon",
            "tec",
        ]
    ).copy()

    # Duplicate bins are collapsed conservatively using the median.
    agg = {
        "tec":
            "median",
    }

    if "dtec" in df.columns:
        agg[
            "dtec"
        ] = "median"

    df = (
        df.groupby(
            [
                "epoch",
                "gdlat",
                "glon",
            ],
            as_index=False,
        )
        .agg(
            agg
        )
        .sort_values(
            [
                "epoch",
                "gdlat",
                "glon",
            ]
        )
        .reset_index(
            drop=True
        )
    )

    if (
        year is not None
        and doy is not None
        and not df.empty
    ):
        target_start = datetime_from_year_doy(
            year,
            doy,
        )

        target_end = (
            target_start
            + timedelta(
                days=1
            )
        )

        epoch_utc = pd.to_datetime(
            df[
                "epoch"
            ],
            utc=True,
            errors="coerce",
        )

        exact_day = (
            (
                epoch_utc
                >= target_start
            )
            & (
                epoch_utc
                < target_end
            )
        )

        df = df[
            exact_day
        ].copy()

    return df


# =============================================================================
# MADRIGAL DAILY SITE-LIST CHECK
# =============================================================================

def download_site_list(
    client,
    site_file,
    year: int,
    doy: int,
    fullname: str,
    email: str,
    affiliation: str,
):
    if site_file is None:
        return None

    path = site_cache_file(
        year,
        doy,
    )

    if (
        path.is_file()
        and path.stat().st_size > 20
    ):
        return path

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    client.isprint(
        site_file.name,
        "gps_site,gdlatr,gdlonr",
        "",
        fullname,
        email,
        affiliation,
        outputFile=str(
            path
        ),
    )

    if (
        not path.is_file()
        or path.stat().st_size < 5
    ):
        return None

    return path


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


def inspect_site_list(
    path: Path | None,
):
    result = {
        "site_list_available":
            False,
        "n_madrigal_input_sites":
            np.nan,
        "koh2_exact_input":
            np.nan,
        "nearest_input_site":
            "",
        "nearest_input_distance_km":
            np.nan,
    }

    if path is None:
        return result

    rows = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as f:
        for line in f:
            items = line.strip().split()

            if len(
                items
            ) < 3:
                continue

            site = items[
                0
            ].strip(
                "'\"b"
            )

            try:
                lat = float(
                    items[
                        1
                    ]
                )

                lon = float(
                    items[
                        2
                    ]
                )
            except Exception:
                continue

            rows.append(
                (
                    site,
                    lat,
                    lon,
                )
            )

    if not rows:
        return result

    result[
        "site_list_available"
    ] = True

    result[
        "n_madrigal_input_sites"
    ] = len(
        rows
    )

    exact = any(
        site.upper()
        == STATION.upper()
        for site, _, _ in rows
    )

    result[
        "koh2_exact_input"
    ] = bool(
        exact
    )

    distances = np.asarray(
        [
            great_circle_km(
                KOH2_LAT,
                KOH2_LON,
                lat,
                lon,
            )
            for _, lat, lon in rows
        ],
        dtype=float,
    )

    i = int(
        np.nanargmin(
            distances
        )
    )

    result[
        "nearest_input_site"
    ] = rows[
        i
    ][
        0
    ]

    result[
        "nearest_input_distance_km"
    ] = float(
        distances[
            i
        ]
    )

    return result


# =============================================================================
# STRICT NEAREST MADRIGAL BIN MATCHING
# =============================================================================

def datetime_index_ns_int(
    values,
):
    t = pd.DatetimeIndex(
        pd.to_datetime(
            values,
            utc=True,
            errors="coerce",
        )
    )

    # pandas may preserve us or ns. Normalize explicitly.
    return (
        t.as_unit(
            "ns"
        ).asi8
    )


def nearest_time_indices(
    obs_epoch,
    map_epoch,
):
    obs_ns = datetime_index_ns_int(
        obs_epoch
    )

    map_index = pd.DatetimeIndex(
        pd.to_datetime(
            map_epoch,
            utc=True,
            errors="coerce",
        )
    ).sort_values()

    map_ns = (
        map_index.as_unit(
            "ns"
        ).asi8
    )

    if len(
        map_ns
    ) == 0:
        return (
            np.full(
                len(
                    obs_ns
                ),
                -1,
                dtype=int,
            ),
            np.full(
                len(
                    obs_ns
                ),
                np.nan,
                dtype=float,
            ),
            map_index,
        )

    pos = np.searchsorted(
        map_ns,
        obs_ns,
        side="left",
    )

    lo = np.clip(
        pos
        - 1,
        0,
        len(
            map_ns
        )
        - 1,
    )

    hi = np.clip(
        pos,
        0,
        len(
            map_ns
        )
        - 1,
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

    invalid_obs = (
        obs_ns
        == np.iinfo(
            np.int64
        ).min
    )

    chosen[
        invalid_obs
    ] = -1

    dt_s[
        invalid_obs
    ] = np.nan

    return (
        chosen,
        dt_s,
        map_index,
    )


def nearest_spatial_indices(
    target_lat,
    target_lon,
    grid_lat,
    grid_lon,
):
    """
    Return nearest available Madrigal cell for all target points in one map
    epoch.  Local equirectangular distance is entirely adequate over the
    small Antarctic regional box used here.
    """
    target_lat = np.asarray(
        target_lat,
        dtype=float,
    )

    target_lon = wrap_lon(
        target_lon
    )

    grid_lat = np.asarray(
        grid_lat,
        dtype=float,
    )

    grid_lon = wrap_lon(
        grid_lon
    )

    n = len(
        target_lat
    )

    idx_out = np.full(
        n,
        -1,
        dtype=int,
    )

    dist_out = np.full(
        n,
        np.nan,
        dtype=float,
    )

    if (
        n == 0
        or len(
            grid_lat
        )
        == 0
    ):
        return (
            idx_out,
            dist_out,
        )

    # Chunk target points to cap temporary matrix memory.
    chunk = 500

    for start in range(
        0,
        n,
        chunk,
    ):
        end = min(
            n,
            start
            + chunk,
        )

        lat = target_lat[
            start:end
        ][
            :,
            None,
        ]

        lon = target_lon[
            start:end
        ][
            :,
            None,
        ]

        glat = grid_lat[
            None,
            :,
        ]

        glon = grid_lon[
            None,
            :,
        ]

        dlat_km = (
            glat
            - lat
        ) * 111.32

        dlon_deg = (
            (
                glon
                - lon
                + 180.0
            )
            % 360.0
        ) - 180.0

        mean_lat = np.radians(
            (
                glat
                + lat
            )
            / 2.0
        )

        dlon_km = (
            dlon_deg
            * 111.32
            * np.cos(
                mean_lat
            )
        )

        d2 = (
            dlat_km
            * dlat_km
            + dlon_km
            * dlon_km
        )

        local_idx = np.argmin(
            d2,
            axis=1,
        )

        local_dist = np.sqrt(
            d2[
                np.arange(
                    end
                    - start
                ),
                local_idx,
            ]
        )

        idx_out[
            start:end
        ] = local_idx

        dist_out[
            start:end
        ] = local_dist

    return (
        idx_out,
        dist_out,
    )


def match_to_madrigal(
    observations: pd.DataFrame,
    madrigal: pd.DataFrame,
    lat_col: str,
    lon_col: str,
    method_col: str,
):
    out = observations.copy()

    out[
        "madrigal_vtec"
    ] = np.nan

    out[
        "madrigal_dtec"
    ] = np.nan

    out[
        "madrigal_time_offset_s"
    ] = np.nan

    out[
        "madrigal_distance_km"
    ] = np.nan

    if (
        out.empty
        or madrigal.empty
    ):
        return out

    unique_map_times = pd.DatetimeIndex(
        madrigal[
            "epoch"
        ].dropna().unique()
    ).sort_values()

    (
        time_idx,
        time_offset_s,
        _,
    ) = nearest_time_indices(
        out[
            "epoch"
        ],
        unique_map_times,
    )

    out[
        "madrigal_time_offset_s"
    ] = time_offset_s

    time_valid = (
        time_idx
        >= 0
    ) & (
        time_offset_s
        <= MAX_TIME_DIFFERENCE_S
    )

    # Work one 5-minute map epoch at a time.
    for map_i in np.unique(
        time_idx[
            time_valid
        ]
    ):
        obs_mask = (
            time_valid
            & (
                time_idx
                == map_i
            )
        )

        obs_positions = np.flatnonzero(
            obs_mask
        )

        if len(
            obs_positions
        ) == 0:
            continue

        map_time = unique_map_times[
            map_i
        ]

        grid = madrigal[
            madrigal[
                "epoch"
            ]
            == map_time
        ]

        if grid.empty:
            continue

        nearest_idx, distances = nearest_spatial_indices(
            pd.to_numeric(
                out.iloc[
                    obs_positions
                ][
                    lat_col
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
            pd.to_numeric(
                out.iloc[
                    obs_positions
                ][
                    lon_col
                ],
                errors="coerce",
            ).to_numpy(
                dtype=float
            ),
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

        accepted = (
            nearest_idx
            >= 0
        ) & (
            distances
            <= MAX_SPATIAL_DISTANCE_KM
        )

        if not np.any(
            accepted
        ):
            continue

        accepted_obs_pos = obs_positions[
            accepted
        ]

        accepted_grid_idx = nearest_idx[
            accepted
        ]

        grid_tec = grid[
            "tec"
        ].to_numpy(
            dtype=float
        )

        out.iloc[
            accepted_obs_pos,
            out.columns.get_loc(
                "madrigal_vtec"
            ),
        ] = grid_tec[
            accepted_grid_idx
        ]

        if "dtec" in grid.columns:
            grid_dtec = grid[
                "dtec"
            ].to_numpy(
                dtype=float
            )

            out.iloc[
                accepted_obs_pos,
                out.columns.get_loc(
                    "madrigal_dtec"
                ),
            ] = grid_dtec[
                accepted_grid_idx
            ]

        out.iloc[
            accepted_obs_pos,
            out.columns.get_loc(
                "madrigal_distance_km"
            ),
        ] = distances[
            accepted
        ]

    out[
        "residual_tecu"
    ] = (
        pd.to_numeric(
            out[
                method_col
            ],
            errors="coerce",
        )
        - pd.to_numeric(
            out[
                "madrigal_vtec"
            ],
            errors="coerce",
        )
    )

    return out


# =============================================================================
# DAILY COMPARISONS
# =============================================================================

def daily_stat_row(
    year,
    doy,
    comparison,
    matched,
    method_col,
    n_candidates,
):
    valid = matched.dropna(
        subset=[
            method_col,
            "madrigal_vtec",
        ]
    )

    stats = calculate_statistics(
        valid[
            method_col
        ],
        valid[
            "madrigal_vtec"
        ],
    )

    date = datetime_from_year_doy(
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
        "comparison":
            comparison,
        "n_candidates":
            int(
                n_candidates
            ),
        "match_fraction":
            (
                float(
                    stats[
                        "n"
                    ]
                )
                / float(
                    n_candidates
                )
                if n_candidates
                else np.nan
            ),
        **stats,
        "method_median_tecu":
            (
                float(
                    np.nanmedian(
                        valid[
                            method_col
                        ]
                    )
                )
                if len(
                    valid
                )
                else np.nan
            ),
        "madrigal_median_tecu":
            (
                float(
                    np.nanmedian(
                        valid[
                            "madrigal_vtec"
                        ]
                    )
                )
                if len(
                    valid
                )
                else np.nan
            ),
        "median_madrigal_dtec_tecu":
            (
                float(
                    np.nanmedian(
                        valid[
                            "madrigal_dtec"
                        ]
                    )
                )
                if (
                    len(
                        valid
                    )
                    and np.any(
                        np.isfinite(
                            valid[
                                "madrigal_dtec"
                            ]
                        )
                    )
                )
                else np.nan
            ),
        "median_time_offset_s":
            (
                float(
                    np.nanmedian(
                        valid[
                            "madrigal_time_offset_s"
                        ]
                    )
                )
                if len(
                    valid
                )
                else np.nan
            ),
        "median_spatial_distance_km":
            (
                float(
                    np.nanmedian(
                        valid[
                            "madrigal_distance_km"
                        ]
                    )
                )
                if len(
                    valid
                )
                else np.nan
            ),
    }


def compare_day(
    year,
    doy,
    pt,
    po,
    madrigal,
):
    rows = []
    matched_products = {}

    if pt is not None:
        # PyTECGg VEq: one value per epoch, evaluated at KOH2 itself.
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
            .sort_values(
                "epoch"
            )
            .copy()
        )

        station[
            "station_lat"
        ] = KOH2_LAT

        station[
            "station_lon"
        ] = KOH2_LON

        m_station = match_to_madrigal(
            station,
            madrigal,
            "station_lat",
            "station_lon",
            "veq",
        )

        rows.append(
            daily_stat_row(
                year,
                doy,
                "PyTECGg_VEq_vs_Madrigal_station",
                m_station,
                "veq",
                len(
                    station
                ),
            )
        )

        matched_products[
            "PyTECGg_VEq"
        ] = m_station

        # PyTECGg VTEC at IPP.
        ipp = pt[
            [
                "epoch",
                "sv",
                "lat_ipp",
                "lon_ipp",
                "ele",
                "vtec",
            ]
        ].copy()

        ipp = ipp[
            pd.to_numeric(
                ipp[
                    "ele"
                ],
                errors="coerce",
            )
            >= MIN_ELEVATION_DEG
        ]

        ipp = ipp.dropna(
            subset=[
                "epoch",
                "lat_ipp",
                "lon_ipp",
                "vtec",
            ]
        )

        m_pt = match_to_madrigal(
            ipp,
            madrigal,
            "lat_ipp",
            "lon_ipp",
            "vtec",
        )

        rows.append(
            daily_stat_row(
                year,
                doy,
                "PyTECGg_VTEC_vs_Madrigal_IPP",
                m_pt,
                "vtec",
                len(
                    ipp
                ),
            )
        )

        matched_products[
            "PyTECGg_VTEC"
        ] = m_pt

    if po is not None:
        ipp = po[
            [
                "epoch",
                "sat",
                "lat_ipp",
                "lon_ipp",
                "elevation",
                "vtec",
            ]
        ].copy()

        ipp = ipp[
            pd.to_numeric(
                ipp[
                    "elevation"
                ],
                errors="coerce",
            )
            >= MIN_ELEVATION_DEG
        ]

        ipp = ipp.dropna(
            subset=[
                "epoch",
                "lat_ipp",
                "lon_ipp",
                "vtec",
            ]
        )

        m_po = match_to_madrigal(
            ipp,
            madrigal,
            "lat_ipp",
            "lon_ipp",
            "vtec",
        )

        rows.append(
            daily_stat_row(
                year,
                doy,
                "pyOASIS_VTEC_vs_Madrigal_IPP",
                m_po,
                "vtec",
                len(
                    ipp
                ),
            )
        )

        matched_products[
            "pyOASIS_VTEC"
        ] = m_po

    return (
        rows,
        matched_products,
    )


# =============================================================================
# YEARLY SUMMARIES
# =============================================================================

def equal_day_summary(
    daily,
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
        rows.append({
            "year":
                int(
                    year
                ),
            "comparison":
                comparison,
            "n_days":
                int(
                    (
                        pd.to_numeric(
                            g[
                                "n"
                            ],
                            errors="coerce",
                        )
                        > 0
                    ).sum()
                ),
            "n_points_total":
                int(
                    pd.to_numeric(
                        g[
                            "n"
                        ],
                        errors="coerce",
                    ).fillna(
                        0
                    ).sum()
                ),
            "n_candidates_total":
                int(
                    pd.to_numeric(
                        g[
                            "n_candidates"
                        ],
                        errors="coerce",
                    ).fillna(
                        0
                    ).sum()
                ),
            "overall_match_fraction":
                (
                    float(
                        pd.to_numeric(
                            g[
                                "n"
                            ],
                            errors="coerce",
                        ).fillna(
                            0
                        ).sum()
                    )
                    / float(
                        pd.to_numeric(
                            g[
                                "n_candidates"
                            ],
                            errors="coerce",
                        ).fillna(
                            0
                        ).sum()
                    )
                    if pd.to_numeric(
                        g[
                            "n_candidates"
                        ],
                        errors="coerce",
                    ).fillna(
                        0
                    ).sum()
                    > 0
                    else np.nan
                ),
            "median_daily_match_fraction":
                g[
                    "match_fraction"
                ].median(),
            "mean_daily_bias_tecu":
                g[
                    "bias_tecu"
                ].mean(),
            "median_daily_bias_tecu":
                g[
                    "bias_tecu"
                ].median(),
            "mean_daily_rmse_tecu":
                g[
                    "rmse_tecu"
                ].mean(),
            "median_daily_rmse_tecu":
                g[
                    "rmse_tecu"
                ].median(),
            "median_daily_pearson_r":
                g[
                    "pearson_r"
                ].median(),
            "median_daily_madrigal_dtec_tecu":
                g[
                    "median_madrigal_dtec_tecu"
                ].median(),
            "median_time_offset_s":
                g[
                    "median_time_offset_s"
                ].median(),
            "median_spatial_distance_km":
                g[
                    "median_spatial_distance_km"
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


def point_weighted_summary(
    accumulators,
):
    rows = []

    for (
        year,
        comparison,
    ), acc in accumulators.items():
        rows.append({
            "year":
                year,
            "comparison":
                comparison,
            **acc.statistics(),
        })

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "comparison",
        ]
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()
    configure_runtime(args)

    MASTER_OUTPUT.mkdir(
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

    pd.set_option(
        "display.float_format",
        lambda x:
            f"{x:.6f}",
    )

    print(
        "=" * 118
    )

    print(
        "KOH2 vs MIT HAYSTACK / MADRIGAL GNSS VTEC"
    )

    print(
        "=" * 118
    )

    print(
        "Root:",
        BASE_ROOT,
    )

    print(
        "Madrigal:",
        MADRIGAL_URL,
    )

    print(
        "Instrument:",
        MADRIGAL_INSTRUMENT,
        "(World-wide GNSS Receiver Network)",
    )

    print(
        "VTEC kindat:",
        MADRIGAL_VTEC_KINDAT,
        "(1 deg x 1 deg x 5 min)",
    )

    print(
        "Match tolerances:"
    )

    print(
        f"  time    <= {MAX_TIME_DIFFERENCE_S:.0f} s"
    )

    print(
        f"  distance<= {MAX_SPATIAL_DISTANCE_KM:.0f} km"
    )

    print(
        "Output:",
        MASTER_OUTPUT,
    )

    print(
        "Madrigal cache:",
        CACHE_ROOT,
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
        "Native HDF5 local filtering:",
        USE_NATIVE_HDF5_LOCAL_FILTER,
    )

    print(
        "=" * 118
    )

    fullname, email, affiliation = get_user_identity()

    client = MADRIGAL_API.MadrigalData(
        MADRIGAL_URL
    )

    daily_rows = []
    availability_rows = []
    accumulators = {}

    years = YEARS

    for year in years:
        (
            pt_files,
            po_files,
            days,
        ) = discover_year_files(
            year
        )

        if SELECTED_DOY is not None:
            days = [
                doy
                for doy in days
                if doy == SELECTED_DOY
            ]

        if not days:
            print(
                f"No KOH2 data days found for {year}"
            )
            continue

        print()
        print(
            "=" * 118
        )

        print(
            f"{year}: PyTECGg={len(pt_files)} days, "
            f"pyOASIS={len(po_files)} days, "
            f"union={len(set(pt_files) | set(po_files))} days"
        )

        print(
            "=" * 118
        )

        for i, doy in enumerate(
            days,
            1,
        ):
            print()
            print(
                "-" * 118
            )

            print(
                f"{i}/{len(days)}  {year} DOY {doy:03d}"
            )

            print(
                "-" * 118
            )

            pt = None
            po = None

            if doy in pt_files:
                print(
                    "  PyTECGg:",
                    pt_files[
                        doy
                    ],
                )

                pt = load_pytecgg(
                    pt_files[
                        doy
                    ]
                )

            if doy in po_files:
                print(
                    "  pyOASIS:",
                    po_files[
                        doy
                    ],
                )

                po = load_pyoasis(
                    po_files[
                        doy
                    ]
                )

            if (
                pt is None
                and po is None
            ):
                print(
                    "  No method data for this day."
                )
                continue

            region = get_daily_region(
                pt,
                po,
            )

            status = {
                "year":
                    year,
                "doy":
                    doy,
                "date":
                    datetime_from_year_doy(
                        year,
                        doy,
                    ).date().isoformat(),
                "pytecgg_present":
                    pt is not None,
                "pyoasis_present":
                    po is not None,
                "madrigal_status":
                    "FAILED",
            }

            try:
                day_client, experiment = experiment_for_day(
                    client,
                    year,
                    doy,
                )

                print(
                    "  Madrigal experiment:",
                    getattr(
                        experiment,
                        "name",
                        "",
                    ),
                )

                print(
                    "  Experiment start:",
                    experiment_start_dt(
                        experiment
                    ).isoformat(),
                )

                print(
                    "  Experiment end  :",
                    experiment_end_dt(
                        experiment
                    ).isoformat(),
                )

                vtec_file = select_kindat_file(
                    day_client,
                    experiment,
                    MADRIGAL_VTEC_KINDAT,
                )

                if vtec_file is None:
                    raise RuntimeError(
                        "No Madrigal kindat 3500 gridded VTEC file found."
                    )

                print(
                    "  Madrigal VTEC file:",
                    vtec_file.name,
                )

                site_file = select_kindat_file(
                    day_client,
                    experiment,
                    MADRIGAL_SITE_KINDAT,
                )

                if USE_NATIVE_HDF5_LOCAL_FILTER:
                    native_path = download_native_vtec_hdf5(
                        day_client,
                        vtec_file,
                        year,
                        doy,
                        fullname,
                        email,
                        affiliation,
                    )

                    madrigal = load_native_regional_vtec(
                        native_path,
                        year,
                        doy,
                        region,
                    )

                else:
                    regional_path = download_regional_vtec(
                        day_client,
                        vtec_file,
                        year,
                        doy,
                        region,
                        fullname,
                        email,
                        affiliation,
                    )

                    madrigal = load_regional_vtec(
                        regional_path,
                        year=year,
                        doy=doy,
                    )

                if madrigal.empty:
                    raise RuntimeError(
                        "Madrigal returned no regional TEC bins for the "
                        "requested UTC day."
                    )

                print(
                    "  Madrigal regional rows:",
                    f"{len(madrigal):,}",
                )

                print(
                    "  Madrigal epochs:",
                    madrigal[
                        "epoch"
                    ].nunique(),
                )

                print(
                    "  Madrigal epoch first:",
                    madrigal[
                        "epoch"
                    ].min(),
                )

                print(
                    "  Madrigal epoch last :",
                    madrigal[
                        "epoch"
                    ].max(),
                )

                print(
                    "  Madrigal latitude range:",
                    f"{madrigal['gdlat'].min():.2f} .. "
                    f"{madrigal['gdlat'].max():.2f}",
                )

                print(
                    "  Madrigal longitude range:",
                    f"{madrigal['glon'].min():.2f} .. "
                    f"{madrigal['glon'].max():.2f}",
                )

                # Receiver-list independence check.
                site_path = download_site_list(
                    day_client,
                    site_file,
                    year,
                    doy,
                    fullname,
                    email,
                    affiliation,
                )

                site_info = inspect_site_list(
                    site_path
                )

                status.update(
                    site_info
                )

                if site_info[
                    "site_list_available"
                ]:
                    print(
                        "  Madrigal input sites:",
                        site_info[
                            "n_madrigal_input_sites"
                        ],
                    )

                    print(
                        "  KOH2 exact input:",
                        site_info[
                            "koh2_exact_input"
                        ],
                    )

                    print(
                        "  Nearest Madrigal input receiver:",
                        site_info[
                            "nearest_input_site"
                        ],
                        f"({site_info['nearest_input_distance_km']:.1f} km)",
                    )
                else:
                    print(
                        "  Madrigal daily site list unavailable."
                    )

                rows, products = compare_day(
                    year,
                    doy,
                    pt,
                    po,
                    madrigal,
                )

                for row in rows:
                    daily_rows.append(
                        row
                    )

                    print()
                    print(
                        " ",
                        row[
                            "comparison"
                        ],
                    )

                    print(
                        f"    n       = {row['n']:,} / "
                        f"{row['n_candidates']:,} "
                        f"({100.0 * row['match_fraction']:.1f}% matched)"
                    )

                    print(
                        f"    bias    = {row['bias_tecu']:.6f} TECU"
                    )

                    print(
                        f"    MAE     = {row['mae_tecu']:.6f} TECU"
                    )

                    print(
                        f"    RMSE    = {row['rmse_tecu']:.6f} TECU"
                    )

                    print(
                        f"    r       = {row['pearson_r']:.6f}"
                    )

                    print(
                        f"    median dt = "
                        f"{row['median_time_offset_s']:.1f} s"
                    )

                    print(
                        f"    median distance = "
                        f"{row['median_spatial_distance_km']:.1f} km"
                    )

                    print(
                        f"    median Madrigal DTEC = "
                        f"{row['median_madrigal_dtec_tecu']:.3f} TECU"
                    )

                    key = (
                        year,
                        row[
                            "comparison"
                        ],
                    )

                    if key not in accumulators:
                        accumulators[
                            key
                        ] = GlobalAccumulator()

                    product_key = None
                    method_col = None

                    if row[
                        "comparison"
                    ] == "PyTECGg_VEq_vs_Madrigal_station":
                        product_key = "PyTECGg_VEq"
                        method_col = "veq"

                    elif row[
                        "comparison"
                    ] == "PyTECGg_VTEC_vs_Madrigal_IPP":
                        product_key = "PyTECGg_VTEC"
                        method_col = "vtec"

                    elif row[
                        "comparison"
                    ] == "pyOASIS_VTEC_vs_Madrigal_IPP":
                        product_key = "pyOASIS_VTEC"
                        method_col = "vtec"

                    if (
                        product_key is not None
                        and product_key in products
                    ):
                        g = products[
                            product_key
                        ].dropna(
                            subset=[
                                method_col,
                                "madrigal_vtec",
                            ]
                        )

                        accumulators[
                            key
                        ].add(
                            g[
                                method_col
                            ],
                            g[
                                "madrigal_vtec"
                            ],
                        )

                        if SAVE_MATCHED_POINTS:
                            out_dir = (
                                BY_YEAR_OUTPUT
                                / str(
                                    year
                                )
                                / "DAILY_MATCHED"
                            )

                            out_dir.mkdir(
                                parents=True,
                                exist_ok=True,
                            )

                            g.to_parquet(
                                out_dir
                                / (
                                    f"{STATION}_{year}_{doy:03d}_"
                                    f"{product_key}_vs_Madrigal.parquet"
                                ),
                                index=False,
                            )

                status[
                    "madrigal_status"
                ] = "OK"

                status[
                    "regional_rows"
                ] = len(
                    madrigal
                )

                for row in rows:
                    status[
                        row[
                            "comparison"
                        ]
                        + "_n"
                    ] = row[
                        "n"
                    ]

            except KeyboardInterrupt:
                raise

            except Exception as exc:
                print(
                    "  [MADRIGAL ERROR]",
                    repr(
                        exc
                    ),
                )

                status[
                    "error"
                ] = repr(
                    exc
                )

            availability_rows.append(
                status
            )

            # Restart-friendly snapshots.
            pd.DataFrame(
                daily_rows
            ).to_csv(
                MASTER_OUTPUT
                / (
                    f"{STATION}_2019_2026_"
                    "Madrigal_daily_statistics.csv"
                ),
                index=False,
            )

            pd.DataFrame(
                availability_rows
            ).to_csv(
                MASTER_OUTPUT
                / (
                    f"{STATION}_2019_2026_"
                    "Madrigal_availability.csv"
                ),
                index=False,
            )

    daily_df = pd.DataFrame(
        daily_rows
    )

    availability_df = pd.DataFrame(
        availability_rows
    )

    point_df = point_weighted_summary(
        accumulators
    )

    equal_df = equal_day_summary(
        daily_df
    )

    daily_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "Madrigal_daily_statistics.csv"
        )
    )

    availability_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "Madrigal_availability.csv"
        )
    )

    point_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "Madrigal_yearly_point_weighted_summary.csv"
        )
    )

    equal_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "Madrigal_equal_day_yearly_summary.csv"
        )
    )

    daily_df.to_csv(
        daily_file,
        index=False,
    )

    availability_df.to_csv(
        availability_file,
        index=False,
    )

    point_df.to_csv(
        point_file,
        index=False,
    )

    equal_df.to_csv(
        equal_file,
        index=False,
    )

    print()
    print(
        "=" * 118
    )

    print(
        "MADRIGAL COMPARISON COMPLETE"
    )

    print(
        "=" * 118
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
        "Point-weighted:",
        point_file,
    )

    print(
        "Equal-day:",
        equal_file,
    )

    if not point_df.empty:
        print()
        print(
            "YEARLY POINT-WEIGHTED SUMMARY"
        )

        print(
            "-" * 118
        )

        print(
            point_df.to_string(
                index=False
            )
        )

        print(
            "-" * 118
        )


if __name__ == "__main__":
    main()
