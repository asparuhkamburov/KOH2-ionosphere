from __future__ import annotations

r"""
KOH2 2019-2026 TEC validation against IGS Final GIM (IONEX)
======================================================

Purpose
-------
Validate all available KOH2 observation days from 2019-2026 processed with:

  1) PyTECGg
     ...\MM\DD\PyTECGg_OUTPUT\KOH2_DDD_2025_PyTECGg_VEQ.parquet

  2) pyOASIS
     ...\MM\DD\pyOASIS_OUTPUT\INDICES\TEC\KOH2_DDD_2025_L1L2.TEC

against:

  3) IGS Final Global Ionosphere Maps (GIM), IONEX
     IGS0OPSFIN_YYYYDDD0000_01D_02H_GIM.INX.gz

The script automatically discovers available days in each year. It does NOT assume
that they are consecutive.

Comparisons
-----------
A. PyTECGg VEq vs IGS GIM at the KOH2 station coordinates.
B. PyTECGg per-satellite VTEC vs IGS GIM interpolated at each IPP.
C. pyOASIS per-satellite VTEC vs IGS GIM interpolated at each IPP.

For IPP comparisons, only observations with elevation >= 30 degrees are used.

Interpolation
-------------
- IGS GIM is bilinearly interpolated in latitude/longitude.
- It is linearly interpolated in time between the 2-hour IONEX maps.
- The script also makes a 2-hour-scale PyTECGg VEq vs IGS comparison.

Statistics
----------
Bias = method - IGS
MAE
RMSE
standard deviation of residual
median residual
5th / 95th percentile residual
Pearson correlation r
N matched values

Authoring note
--------------
This is an independent validation script. It does not modify PyTECGg or pyOASIS.
"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import gzip
import math
import netrc
import re
import shutil
import sys
import time
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import polars as pl
import requests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


# ============================================================================
# USER SETTINGS
# ============================================================================

# Default scientific study interval. Runtime selection can be narrowed with CLI.
YEARS = list(range(2019, 2027))
STATION = "KOH2"

# Runtime paths are configured by configure_runtime().
BASE_ROOT = Path(".")
SELECTED_DOY = None

# KOH2 coordinates from the RINEX/header metadata already used in processing.
KOH2_X = 1453335.2992
KOH2_Y = -2554570.1548
KOH2_Z = -5641700.7402

MIN_ELEVATION_DEG = 30.0

# IGS Final combined GIM:
CDDIS_BASE = "https://cddis.nasa.gov/archive/gnss/products/ionex"

MASTER_OUTPUT = BASE_ROOT / "TEC_VALIDATION_IGS_2019_2026"
CACHE_ROOT = MASTER_OUTPUT / "_IGS_IONEX"
BY_YEAR_OUTPUT = MASTER_OUTPUT / "BY_YEAR"

# These globals are configured for each year by configure_year().
YEAR = None
ROOT = None
IONEX_CACHE = None
OUTPUT_ROOT = None
DAILY_OUTPUT = None

# Save detailed matched data day-by-day.
SAVE_DAILY_PARQUET = True

# Number of point pairs retained per day for diagnostic scatter plots.
# This affects only the plots, not the statistics.
SCATTER_SAMPLE_PER_DAY = 1500

# Model-scale station comparison.
MODEL_SCALE = "2h"

# Plotting.
DPI = 180

# CDDIS download robustness for a long 2019-2026 batch.
# Each candidate filename is retried before the script moves on to the next
# naming convention. A failed day is still skipped cleanly by the yearly loop.
DOWNLOAD_ATTEMPTS_PER_CANDIDATE = 3
DOWNLOAD_CONNECT_TIMEOUT_S = 30
DOWNLOAD_READ_TIMEOUT_S = 60
DOWNLOAD_RETRY_WAIT_S = 5
DOWNLOAD_CHUNK_BYTES = 256 * 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate KOH2 PyTECGg and pyOASIS VTEC against "
            "IGS Final GIM (IONEX)."
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
            "Validation output directory. Default: "
            "<data-root>/TEC_VALIDATION_IGS_2019_2026"
        ),
    )
    parser.add_argument(
        "--ionex-cache-root",
        type=Path,
        default=None,
        help=(
            "Optional existing IGS IONEX cache root containing year folders. "
            "Default: <output-root>/_IGS_IONEX"
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
    global BASE_ROOT, MASTER_OUTPUT, CACHE_ROOT, BY_YEAR_OUTPUT
    global YEARS, SELECTED_DOY

    BASE_ROOT = args.data_root.resolve()

    MASTER_OUTPUT = (
        args.output_root.resolve()
        if args.output_root is not None
        else BASE_ROOT / "TEC_VALIDATION_IGS_2019_2026"
    )

    CACHE_ROOT = (
        args.ionex_cache_root.resolve()
        if args.ionex_cache_root is not None
        else MASTER_OUTPUT / "_IGS_IONEX"
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
        SELECTED_DOY = int(
            selected.strftime("%j")
        )


def configure_year(year: int):
    global YEAR, ROOT, IONEX_CACHE, OUTPUT_ROOT, DAILY_OUTPUT

    YEAR = int(year)
    ROOT = BASE_ROOT / str(YEAR)

    IONEX_CACHE = (
        CACHE_ROOT
        / str(YEAR)
    )

    OUTPUT_ROOT = (
        BY_YEAR_OUTPUT
        / str(YEAR)
    )

    DAILY_OUTPUT = (
        OUTPUT_ROOT
        / "DAILY_MATCHED"
    )


# ============================================================================
# EARTHDATA / CDDIS AUTHENTICATION
# ============================================================================

EARTHDATA_HOST = "urs.earthdata.nasa.gov"


def find_netrc_credentials():
    r"""
    Find Earthdata credentials in:
      %USERPROFILE%\_netrc   (normal Windows CDDIS setup)
      %USERPROFILE%\.netrc
    """
    candidates = [
        Path.home() / "_netrc",
        Path.home() / ".netrc",
    ]

    for path in candidates:
        if not path.is_file():
            continue

        try:
            auth = netrc.netrc(str(path)).authenticators(EARTHDATA_HOST)
        except Exception as exc:
            print(f"[WARNING] Could not read {path}: {exc}")
            continue

        if auth:
            username, _, password = auth
            return username, password, path

    return None, None, None


class SessionWithHeaderRedirection(requests.Session):
    """
    Earthdata Login compatible requests session.

    Keeps Authorization when redirecting to/from urs.earthdata.nasa.gov,
    while avoiding forwarding credentials to unrelated hosts.
    """

    AUTH_HOST = EARTHDATA_HOST

    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)

    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers

        if "Authorization" not in headers:
            return

        original_parsed = urlparse(response.request.url)
        redirect_parsed = urlparse(prepared_request.url)

        if (
            original_parsed.hostname != redirect_parsed.hostname
            and redirect_parsed.hostname != self.AUTH_HOST
            and original_parsed.hostname != self.AUTH_HOST
        ):
            del headers["Authorization"]


_CDDIS_SESSION = None


def get_cddis_session():
    global _CDDIS_SESSION

    if _CDDIS_SESSION is not None:
        return _CDDIS_SESSION

    username, password, netrc_path = find_netrc_credentials()

    if not username or not password:
        raise RuntimeError(
            "Earthdata credentials were not found.\n"
            "Expected a Windows _netrc file such as:\n"
            r"  <USER_HOME>\_netrc" "\n"
            "containing:\n"
            "  machine urs.earthdata.nasa.gov login YOUR_LOGIN password YOUR_PASSWORD"
        )

    print(f"Earthdata credentials: {netrc_path}")
    _CDDIS_SESSION = SessionWithHeaderRedirection(username, password)

    # Retry connection-level and selected transient HTTP failures.
    # Manual retries in download_ionex() additionally cover failures that
    # happen while streaming the response body.
    try:
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry

        retry = Retry(
            total=3,
            connect=3,
            read=3,
            status=3,
            backoff_factor=1.0,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET", "HEAD"]),
            raise_on_status=False,
        )

        adapter = HTTPAdapter(
            max_retries=retry,
            pool_connections=4,
            pool_maxsize=4,
        )

        _CDDIS_SESSION.mount(
            "https://",
            adapter,
        )
    except Exception as exc:
        print(
            "[WARNING] Could not configure HTTP retry adapter:",
            repr(exc),
        )

    return _CDDIS_SESSION


# ============================================================================
# STATION COORDINATES
# ============================================================================

def ecef_to_geodetic(x, y, z):
    """WGS84 ECEF -> geodetic latitude, longitude, ellipsoidal height."""
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)

    lon = math.atan2(y, x)
    p = math.hypot(x, y)

    lat = math.atan2(z, p * (1.0 - e2))

    for _ in range(20):
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        h = p / math.cos(lat) - N

        lat_new = math.atan2(
            z,
            p * (1.0 - e2 * N / (N + h))
        )

        if abs(lat_new - lat) < 1e-14:
            lat = lat_new
            break

        lat = lat_new

    sin_lat = math.sin(lat)
    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    h = p / math.cos(lat) - N

    return math.degrees(lat), math.degrees(lon), h


KOH2_LAT, KOH2_LON, KOH2_H = ecef_to_geodetic(
    KOH2_X, KOH2_Y, KOH2_Z
)


# ============================================================================
# DATE / FILE DISCOVERY
# ============================================================================

def date_from_year_doy(year: int, doy: int) -> datetime:
    dt = datetime.strptime(f"{year}-{doy:03d}", "%Y-%j")
    return dt.replace(tzinfo=timezone.utc)


def extract_doy(path: Path):
    m = re.search(r"_(\d{3})_(\d{4})_", path.name)
    if not m:
        return None

    doy = int(m.group(1))
    year = int(m.group(2))

    if year != YEAR:
        return None

    return doy


def discover_data_files():
    pytecgg = {}
    pyoasis = {}

    for path in ROOT.glob(
        rf"*\*\PyTECGg_OUTPUT\{STATION}_*_{YEAR}_PyTECGg_VEQ.parquet"
    ):
        doy = extract_doy(path)
        if doy is not None:
            pytecgg[doy] = path

    for path in ROOT.glob(
        rf"*\*\pyOASIS_OUTPUT\INDICES\TEC\{STATION}_*_{YEAR}_L1L2.TEC"
    ):
        doy = extract_doy(path)
        if doy is not None:
            pyoasis[doy] = path

    all_days = sorted(set(pytecgg) | set(pyoasis))

    return pytecgg, pyoasis, all_days


# ============================================================================
# IGS IONEX DOWNLOAD
# ============================================================================

def ionex_candidate_names(year: int, doy: int):
    """
    Return IGS Final combined GIM candidate names covering both the legacy
    short-name era and the long-name transition.

    IGS switched operational products to long filenames from GPS week 2238
    (late November 2022). During the transition, both the older ION/IOX
    content naming and the newer GIM/INX naming may occur.
    """
    yy = year % 100

    return [
        # Current long filename/content convention.
        f"IGS0OPSFIN_{year}{doy:03d}0000_01D_02H_GIM.INX.gz",

        # Alternate 4-character AC spelling occasionally encountered.
        f"IGS00PSFIN_{year}{doy:03d}0000_01D_02H_GIM.INX.gz",

        # Transitional long content/extension convention.
        f"IGS0OPSFIN_{year}{doy:03d}0000_01D_02H_ION.IOX.gz",
        f"IGS00PSFIN_{year}{doy:03d}0000_01D_02H_ION.IOX.gz",

        # Legacy final combined IONEX.
        f"igsg{doy:03d}0.{yy:02d}i.Z",
    ]


def ionex_url(year: int, doy: int, name: str):
    return f"{CDDIS_BASE}/{year}/{doy:03d}/{name}"


def decompressed_path_for(path: Path) -> Path:
    low = path.name.lower()

    if low.endswith(".gz"):
        return path.with_suffix("")

    if path.name.endswith(".Z"):
        return path.with_suffix("")

    return path


def decompress_ionex(compressed: Path, output: Path):
    low = compressed.name.lower()

    if low.endswith(".gz"):
        with gzip.open(compressed, "rb") as src, open(output, "wb") as dst:
            shutil.copyfileobj(src, dst)

    elif compressed.name.endswith(".Z"):
        try:
            from unlzw3 import unlzw
        except ImportError as exc:
            raise RuntimeError(
                "Historical IGS IONEX files before the long-filename transition "
                "are commonly UNIX-compress .Z files.\n"
                "Install the small decompressor once in pytecgg_env:\n"
                "  python -m pip install unlzw3"
            ) from exc

        raw = compressed.read_bytes()
        output.write_bytes(unlzw(raw))

    else:
        shutil.copyfile(compressed, output)

    if not output.is_file() or output.stat().st_size < 10000:
        raise RuntimeError(
            f"Decompressed IONEX appears invalid or too small: {output}"
        )

    return output


def try_existing_ionex(year: int, doy: int, names):
    r"""
    Reuse either the new multi-year cache or an older year-local cache such as
    ...\2025\_IGS_IONEX created by the original validator.
    """
    search_dirs = [
        IONEX_CACHE,
        ROOT / "_IGS_IONEX",
    ]

    for directory in search_dirs:
        if not directory.is_dir():
            continue

        for name in names:
            compressed = directory / name
            decompressed = decompressed_path_for(compressed)

            if decompressed.is_file() and decompressed.stat().st_size > 10000:
                print(f"  IONEX cached: {decompressed}")
                return decompressed

            if compressed.is_file() and compressed.stat().st_size > 10000:
                IONEX_CACHE.mkdir(parents=True, exist_ok=True)

                local_compressed = IONEX_CACHE / compressed.name

                if compressed.resolve() != local_compressed.resolve():
                    shutil.copy2(compressed, local_compressed)
                else:
                    local_compressed = compressed

                local_output = decompressed_path_for(local_compressed)

                print(f"  Decompressing cached: {local_compressed.name}")
                return decompress_ionex(local_compressed, local_output)

    return None


def download_ionex(year: int, doy: int) -> Path:
    r"""
    Obtain one IGS Final IONEX file.

    The E: production tree is allowed to contain only PyTECGg_OUTPUT and
    pyOASIS_OUTPUT. IONEX products are therefore downloaded into the validator's
    own cache under:

        <IONEX_CACHE_ROOT>\<year>

    Robustness:
    - reuse already decompressed/downloaded files;
    - try all relevant historical/current IGS filename conventions;
    - retry each candidate several times if HTTPS streaming is interrupted;
    - delete only incomplete .part files;
    - leave successful compressed/decompressed cache files for later reruns.
    """
    IONEX_CACHE.mkdir(parents=True, exist_ok=True)

    names = ionex_candidate_names(year, doy)

    existing = try_existing_ionex(
        year,
        doy,
        names,
    )

    if existing is not None:
        return existing

    session = get_cddis_session()
    errors = []

    for name in names:
        url = ionex_url(
            year,
            doy,
            name,
        )

        compressed = IONEX_CACHE / name
        output = decompressed_path_for(compressed)

        for attempt in range(
            1,
            DOWNLOAD_ATTEMPTS_PER_CANDIDATE + 1,
        ):
            print("  Trying IGS Final IONEX:")
            print("   ", url)

            if DOWNLOAD_ATTEMPTS_PER_CANDIDATE > 1:
                print(
                    f"    download attempt "
                    f"{attempt}/{DOWNLOAD_ATTEMPTS_PER_CANDIDATE}"
                )

            tmp = compressed.with_suffix(
                compressed.suffix + ".part"
            )

            # A stale partial download must never be mistaken for a valid file.
            if tmp.exists():
                try:
                    tmp.unlink()
                except Exception:
                    pass

            response = None

            try:
                response = session.get(
                    url,
                    stream=True,
                    timeout=(
                        DOWNLOAD_CONNECT_TIMEOUT_S,
                        DOWNLOAD_READ_TIMEOUT_S,
                    ),
                )

                if response.status_code == 404:
                    errors.append(
                        f"404 {name}"
                    )
                    # A 404 is not transient for this filename: go directly
                    # to the next naming convention.
                    break

                response.raise_for_status()

                content_type = response.headers.get(
                    "Content-Type",
                    "",
                ).lower()

                bytes_written = 0
                first_chunk = True

                with open(tmp, "wb") as f:
                    for chunk in response.iter_content(
                        chunk_size=DOWNLOAD_CHUNK_BYTES
                    ):
                        if not chunk:
                            continue

                        if first_chunk:
                            head = chunk[:500].lower()

                            if (
                                b"<html" in head
                                or b"earthdata login" in head
                                or "text/html" in content_type
                            ):
                                raise RuntimeError(
                                    "CDDIS returned an HTML page instead of "
                                    "IONEX. Check Earthdata authentication."
                                )

                            first_chunk = False

                        f.write(chunk)
                        bytes_written += len(chunk)

                if bytes_written < 1000:
                    raise RuntimeError(
                        f"Downloaded file is too small "
                        f"({bytes_written} bytes): {name}"
                    )

                tmp.replace(compressed)

                if compressed.stat().st_size < 1000:
                    raise RuntimeError(
                        f"Downloaded file is too small: {compressed}"
                    )

                print(
                    "  Downloaded:",
                    compressed.name,
                    f"({compressed.stat().st_size:,} bytes)",
                )

                return decompress_ionex(
                    compressed,
                    output,
                )

            except KeyboardInterrupt:
                # Preserve normal Ctrl+C behavior rather than hiding it.
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except Exception:
                        pass
                raise

            except Exception as exc:
                message = (
                    f"{name}, attempt {attempt}/"
                    f"{DOWNLOAD_ATTEMPTS_PER_CANDIDATE}: {exc!r}"
                )

                errors.append(
                    message
                )

                print(
                    "    [DOWNLOAD ERROR]",
                    repr(exc),
                )

                if tmp.exists():
                    try:
                        tmp.unlink()
                    except Exception:
                        pass

                if attempt < DOWNLOAD_ATTEMPTS_PER_CANDIDATE:
                    print(
                        f"    Retrying in "
                        f"{DOWNLOAD_RETRY_WAIT_S} s..."
                    )
                    time.sleep(
                        DOWNLOAD_RETRY_WAIT_S
                    )

            finally:
                if response is not None:
                    try:
                        response.close()
                    except Exception:
                        pass

    raise RuntimeError(
        f"No IGS Final IONEX candidate could be obtained for "
        f"{year} DOY {doy:03d}.\n"
        + "\n".join(
            f"    {x}"
            for x in errors
        )
    )


# ============================================================================
# IONEX PARSER
# ============================================================================

class IonexGrid:
    def __init__(self, epochs, latitudes, longitudes, tec):
        self.epochs = epochs
        self.latitudes = np.asarray(latitudes, dtype=float)
        self.longitudes = np.asarray(longitudes, dtype=float)
        self.tec = np.asarray(tec, dtype=float)

        self.time_seconds = np.asarray(
            [e.timestamp() for e in self.epochs],
            dtype=float,
        )


def _line_label(line: str) -> str:
    if len(line) < 61:
        return ""
    return line[60:].strip()


IONEX_NUMBER_RE = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
)


def ionex_numbers(line: str, limit: int = 60):
    """
    Robust numeric extraction for fixed-width IONEX records.

    IONEX can legally contain adjacent signed fields, e.g.:
        87.5-180.0
    which cannot be parsed safely with whitespace splitting.
    """
    return [
        float(x)
        for x in IONEX_NUMBER_RE.findall(
            line[:limit]
        )
    ]


def read_ionex(path: Path) -> IonexGrid:
    """
    Read 2-D TEC maps from an IONEX file.

    Values marked 9999 are treated as missing.
    The IONEX EXPONENT is applied to all TEC values.
    """
    with open(path, "r", encoding="ascii", errors="replace") as f:
        lines = f.readlines()

    exponent = 0
    i = 0

    # Header
    while i < len(lines):
        line = lines[i]
        label = _line_label(line)

        if label == "EXPONENT":
            vals = ionex_numbers(line)

            if vals:
                exponent = int(
                    round(vals[0])
                )

        if label == "END OF HEADER":
            i += 1
            break

        i += 1

    factor = 10.0 ** exponent

    map_epochs = []
    map_matrices = []
    master_lats = None
    master_lons = None

    while i < len(lines):
        line = lines[i]
        label = _line_label(line)

        if label != "START OF TEC MAP":
            i += 1
            continue

        i += 1
        epoch = None
        rows = []

        while i < len(lines):
            line = lines[i]
            label = _line_label(line)

            if label == "EPOCH OF CURRENT MAP":
                vals = ionex_numbers(
                    line,
                    limit=36,
                )

                if len(vals) >= 6:
                    year, month, day, hour, minute, second = [
                        int(round(v))
                        for v in vals[:6]
                    ]

                    epoch = datetime(
                        year,
                        month,
                        day,
                        hour,
                        minute,
                        second,
                        tzinfo=timezone.utc,
                    )

                i += 1
                continue

            if label == "LAT/LON1/LON2/DLON/H":
                vals = ionex_numbers(line)

                if len(vals) < 5:
                    raise RuntimeError(
                        f"Invalid LAT/LON row in {path.name}: "
                        f"{line.rstrip()}"
                    )

                lat, lon1, lon2, dlon, _height = vals[:5]

                nlon = (
                    int(
                        round(
                            (lon2 - lon1)
                            / dlon
                        )
                    )
                    + 1
                )

                raw_values = []
                i += 1

                while (
                    i < len(lines)
                    and len(raw_values) < nlon
                ):
                    data_line = lines[i]
                    data_label = _line_label(data_line)

                    if data_label in {
                        "LAT/LON1/LON2/DLON/H",
                        "END OF TEC MAP",
                        "EPOCH OF CURRENT MAP",
                    }:
                        break

                    for tok in data_line.split():
                        try:
                            raw_values.append(
                                float(tok)
                            )
                        except ValueError:
                            pass

                    i += 1

                raw = np.asarray(
                    raw_values[:nlon],
                    dtype=float,
                )

                if raw.size != nlon:
                    raise RuntimeError(
                        f"Incomplete IONEX longitude row in {path.name}: "
                        f"expected {nlon}, got {raw.size}"
                    )

                # IONEX uses the exact raw integer value 9999 as the
                # missing-TEC sentinel.  Do not use np.isclose() here:
                # with some recent NumPy 2.x builds that call can fail in
                # this long-running environment even though raw is numeric.
                # A direct comparison is also scientifically more appropriate
                # because 9999 is a discrete format flag, not an approximate
                # measurement.
                missing = (
                    raw
                    == 9999.0
                )

                values = raw * factor
                values[missing] = np.nan

                lons = (
                    lon1
                    + np.arange(
                        nlon,
                        dtype=float,
                    )
                    * dlon
                )

                rows.append(
                    (
                        lat,
                        lons,
                        values,
                    )
                )

                continue

            if label == "END OF TEC MAP":
                i += 1
                break

            i += 1

        if epoch is None or not rows:
            continue

        # Sort latitude ascending.
        rows.sort(
            key=lambda r:
                r[0]
        )

        lats = np.asarray(
            [
                r[0]
                for r in rows
            ],
            dtype=float,
        )

        lons = np.asarray(
            rows[0][1],
            dtype=float,
        )

        matrix = np.vstack(
            [
                r[2]
                for r in rows
            ]
        )

        # Ensure longitude ascending.
        if lons[0] > lons[-1]:
            lons = lons[::-1]
            matrix = matrix[:, ::-1]

        if master_lats is None:
            master_lats = lats
            master_lons = lons

        else:
            if (
                len(lats) != len(master_lats)
                or len(lons) != len(master_lons)
                or not np.allclose(
                    lats,
                    master_lats,
                )
                or not np.allclose(
                    lons,
                    master_lons,
                )
            ):
                raise RuntimeError(
                    f"IONEX grid changed inside {path.name}; "
                    "this validator expects one common 2-D grid."
                )

        map_epochs.append(epoch)
        map_matrices.append(matrix)

    if not map_epochs:
        raise RuntimeError(
            f"No TEC maps found in IONEX file: {path}"
        )

    order = np.argsort(
        np.asarray(
            [
                e.timestamp()
                for e in map_epochs
            ]
        )
    )

    map_epochs = [
        map_epochs[j]
        for j in order
    ]

    map_matrices = [
        map_matrices[j]
        for j in order
    ]

    return IonexGrid(
        map_epochs,
        master_lats,
        master_lons,
        np.stack(
            map_matrices,
            axis=0,
        ),
    )


# ============================================================================
# IONEX SPATIAL / TEMPORAL INTERPOLATION
# ============================================================================

def normalize_longitudes(lon, grid_lons):
    lon = np.asarray(lon, dtype=float)

    lo = float(grid_lons[0])
    hi = float(grid_lons[-1])
    span = hi - lo

    if span >= 359.0:
        out = ((lon - lo) % 360.0) + lo

        # Keep values inside [lo, hi].
        out = np.where(out > hi, out - 360.0, out)
        return out

    return lon


def bilinear_interpolate(grid2d, grid_lats, grid_lons, query_lats, query_lons):
    query_lats = np.asarray(query_lats, dtype=float)
    query_lons = normalize_longitudes(query_lons, grid_lons)

    out = np.full(query_lats.shape, np.nan, dtype=float)

    finite = np.isfinite(query_lats) & np.isfinite(query_lons)

    finite &= query_lats >= grid_lats[0]
    finite &= query_lats <= grid_lats[-1]
    finite &= query_lons >= grid_lons[0]
    finite &= query_lons <= grid_lons[-1]

    if not np.any(finite):
        return out

    latq = query_lats[finite]
    lonq = query_lons[finite]

    iy1 = np.searchsorted(grid_lats, latq, side="right")
    ix1 = np.searchsorted(grid_lons, lonq, side="right")

    iy1 = np.clip(iy1, 1, len(grid_lats) - 1)
    ix1 = np.clip(ix1, 1, len(grid_lons) - 1)

    iy0 = iy1 - 1
    ix0 = ix1 - 1

    y0 = grid_lats[iy0]
    y1 = grid_lats[iy1]
    x0 = grid_lons[ix0]
    x1 = grid_lons[ix1]

    fy = np.divide(
        latq - y0,
        y1 - y0,
        out=np.zeros_like(latq),
        where=(y1 != y0),
    )

    fx = np.divide(
        lonq - x0,
        x1 - x0,
        out=np.zeros_like(lonq),
        where=(x1 != x0),
    )

    q00 = grid2d[iy0, ix0]
    q01 = grid2d[iy0, ix1]
    q10 = grid2d[iy1, ix0]
    q11 = grid2d[iy1, ix1]

    weights = np.vstack([
        (1.0 - fy) * (1.0 - fx),
        (1.0 - fy) * fx,
        fy * (1.0 - fx),
        fy * fx,
    ])

    values = np.vstack([q00, q01, q10, q11])

    valid = np.isfinite(values)

    weighted_sum = np.nansum(
        np.where(valid, values * weights, 0.0),
        axis=0,
    )

    weight_sum = np.sum(
        np.where(valid, weights, 0.0),
        axis=0,
    )

    result = np.divide(
        weighted_sum,
        weight_sum,
        out=np.full_like(weighted_sum, np.nan),
        where=weight_sum > 0,
    )

    out[finite] = result
    return out


def interpolate_ionex(grid: IonexGrid, epochs, lats, lons):
    """
    Bilinear spatial + linear temporal interpolation.

    pandas compatibility note
    -------------------------
    pd.to_datetime() returns a Series when the input is a Series, but returns
    a DatetimeIndex for several other input types. Only DatetimeIndex exposes
    .asi8 directly. Convert explicitly to DatetimeIndex so this function works
    identically with both PyTECGg and pyOASIS pandas Series.
    """
    t_parsed = pd.to_datetime(
        epochs,
        utc=True,
        errors="coerce",
    )

    t = pd.DatetimeIndex(
        t_parsed
    )

    valid_time = ~t.isna()

    t_seconds = np.full(
        len(t),
        np.nan,
        dtype=float,
    )

    if np.any(valid_time):
        # IMPORTANT:
        # pandas preserves the native datetime resolution of Parquet data.
        # The production PyTECGg files use datetime64[us, UTC], whereas
        # timestamps generated from pyOASIS MJD can use datetime64[ns, UTC].
        #
        # DatetimeIndex.asi8 is expressed in the index's native resolution:
        #   datetime64[us] -> microseconds
        #   datetime64[ns] -> nanoseconds
        #
        # Therefore dividing asi8 by 1e9 without first fixing the unit makes
        # PyTECGg epochs about 1000 times too small and places them outside
        # the IONEX time span. Normalize explicitly to nanoseconds first.
        t_ns = t.as_unit(
            "ns"
        )

        t_seconds[valid_time] = (
            t_ns.asi8[valid_time].astype(np.float64) / 1e9
        )

    lats = np.asarray(lats, dtype=float)
    lons = np.asarray(lons, dtype=float)

    out = np.full(len(t), np.nan, dtype=float)

    valid = (
        np.isfinite(t_seconds)
        & np.isfinite(lats)
        & np.isfinite(lons)
    )

    valid &= t_seconds >= grid.time_seconds[0]
    valid &= t_seconds <= grid.time_seconds[-1]

    if not np.any(valid):
        return out

    idx = np.where(valid)[0]
    tv = t_seconds[idx]

    hi = np.searchsorted(
        grid.time_seconds,
        tv,
        side="right",
    )

    hi = np.clip(hi, 1, len(grid.time_seconds) - 1)
    lo = hi - 1

    # Group by temporal bracket for efficient vectorized spatial interpolation.
    pairs = np.column_stack((lo, hi))
    unique_pairs = np.unique(pairs, axis=0)

    for lo_i, hi_i in unique_pairs:
        group_mask = (lo == lo_i) & (hi == hi_i)
        target_idx = idx[group_mask]

        v0 = bilinear_interpolate(
            grid.tec[lo_i],
            grid.latitudes,
            grid.longitudes,
            lats[target_idx],
            lons[target_idx],
        )

        v1 = bilinear_interpolate(
            grid.tec[hi_i],
            grid.latitudes,
            grid.longitudes,
            lats[target_idx],
            lons[target_idx],
        )

        t0 = grid.time_seconds[lo_i]
        t1 = grid.time_seconds[hi_i]

        alpha = (t_seconds[target_idx] - t0) / (t1 - t0)

        both = np.isfinite(v0) & np.isfinite(v1)
        only0 = np.isfinite(v0) & ~np.isfinite(v1)
        only1 = ~np.isfinite(v0) & np.isfinite(v1)

        vals = np.full(len(target_idx), np.nan, dtype=float)

        vals[both] = (
            v0[both] * (1.0 - alpha[both])
            + v1[both] * alpha[both]
        )

        vals[only0] = v0[only0]
        vals[only1] = v1[only1]

        out[target_idx] = vals

    return out


# ============================================================================
# DATA READERS
# ============================================================================

def polars_to_pandas(df: pl.DataFrame) -> pd.DataFrame:
    try:
        return df.to_pandas()
    except Exception:
        return pd.DataFrame(df.to_dicts())


def load_pytecgg(path: Path):
    df_pl = pl.read_parquet(
        path,
        columns=[
            "epoch",
            "sv",
            "lat_ipp",
            "lon_ipp",
            "ele",
            "vtec",
            "veq",
        ],
    )

    df = polars_to_pandas(df_pl)

    df["epoch"] = pd.to_datetime(
        df["epoch"],
        utc=True,
        errors="coerce",
    )

    for col in ["lat_ipp", "lon_ipp", "ele", "vtec", "veq"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def load_pyoasis(path: Path):
    names = [
        "sat",
        "mjd",
        "stec",
        "lon_ipp",
        "lat_ipp",
        "elevation",
        "azimuth",
        "vtec",
    ]

    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=names,
        engine="python",
    )

    for col in names:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["epoch"] = pd.to_datetime(
        df["mjd"],
        unit="D",
        origin="1858-11-17",
        utc=True,
        errors="coerce",
    )

    return df


# ============================================================================
# STATISTICS
# ============================================================================

def clean_pair(method, reference):
    method = np.asarray(method, dtype=float)
    reference = np.asarray(reference, dtype=float)

    good = np.isfinite(method) & np.isfinite(reference)

    return method[good], reference[good]


def calculate_statistics(method, reference):
    method, reference = clean_pair(method, reference)

    n = len(method)

    result = {
        "n": int(n),
        "bias_tecu": np.nan,
        "mae_tecu": np.nan,
        "rmse_tecu": np.nan,
        "std_residual_tecu": np.nan,
        "median_residual_tecu": np.nan,
        "p05_residual_tecu": np.nan,
        "p95_residual_tecu": np.nan,
        "pearson_r": np.nan,
    }

    if n == 0:
        return result

    resid = method - reference

    result["bias_tecu"] = float(np.mean(resid))
    result["mae_tecu"] = float(np.mean(np.abs(resid)))
    result["rmse_tecu"] = float(np.sqrt(np.mean(resid ** 2)))
    result["std_residual_tecu"] = float(np.std(resid))
    result["median_residual_tecu"] = float(np.median(resid))
    result["p05_residual_tecu"] = float(np.quantile(resid, 0.05))
    result["p95_residual_tecu"] = float(np.quantile(resid, 0.95))

    if (
        n >= 2
        and np.std(method) > 0
        and np.std(reference) > 0
    ):
        result["pearson_r"] = float(
            np.corrcoef(method, reference)[0, 1]
        )

    return result


class GlobalAccumulator:
    """
    Incremental exact sums for mean/RMSE/correlation plus residual chunks for
    global median and percentiles.
    """

    def __init__(self):
        self.n = 0

        self.sum_m = 0.0
        self.sum_r = 0.0
        self.sum_m2 = 0.0
        self.sum_r2 = 0.0
        self.sum_mr = 0.0

        self.sum_res = 0.0
        self.sum_abs_res = 0.0
        self.sum_res2 = 0.0

        self.residual_chunks = []

    def add(self, method, reference):
        method, reference = clean_pair(method, reference)

        if len(method) == 0:
            return

        resid = method - reference

        self.n += len(method)

        self.sum_m += float(method.sum())
        self.sum_r += float(reference.sum())
        self.sum_m2 += float(np.dot(method, method))
        self.sum_r2 += float(np.dot(reference, reference))
        self.sum_mr += float(np.dot(method, reference))

        self.sum_res += float(resid.sum())
        self.sum_abs_res += float(np.abs(resid).sum())
        self.sum_res2 += float(np.dot(resid, resid))

        self.residual_chunks.append(
            resid.astype(np.float32, copy=False)
        )

    def statistics(self):
        if self.n == 0:
            return calculate_statistics([], [])

        n = float(self.n)

        bias = self.sum_res / n
        mae = self.sum_abs_res / n
        rmse = math.sqrt(self.sum_res2 / n)

        var_res = max(
            0.0,
            self.sum_res2 / n - bias * bias,
        )

        std_res = math.sqrt(var_res)

        numerator = (
            self.sum_mr - self.sum_m * self.sum_r / n
        )

        denom_m = (
            self.sum_m2 - self.sum_m * self.sum_m / n
        )

        denom_r = (
            self.sum_r2 - self.sum_r * self.sum_r / n
        )

        if denom_m > 0 and denom_r > 0:
            corr = numerator / math.sqrt(denom_m * denom_r)
        else:
            corr = np.nan

        residuals = np.concatenate(self.residual_chunks)

        return {
            "n": int(self.n),
            "bias_tecu": bias,
            "mae_tecu": mae,
            "rmse_tecu": rmse,
            "std_residual_tecu": std_res,
            "median_residual_tecu": float(np.median(residuals)),
            "p05_residual_tecu": float(np.quantile(residuals, 0.05)),
            "p95_residual_tecu": float(np.quantile(residuals, 0.95)),
            "pearson_r": float(corr) if np.isfinite(corr) else np.nan,
        }


# ============================================================================
# OUTPUT HELPERS
# ============================================================================

def save_dataframe(df: pd.DataFrame, path: Path):
    if not SAVE_DAILY_PARQUET:
        return

    try:
        df.to_parquet(
            path,
            index=False,
        )
    except Exception as exc:
        fallback = path.with_suffix(".csv.gz")

        print(
            f"  Parquet save failed ({exc}); "
            f"saving gzip CSV instead."
        )

        df.to_csv(
            fallback,
            index=False,
            compression="gzip",
        )


def add_daily_stat(rows, date, comparison, method, reference):
    method_clean, reference_clean = clean_pair(
        method,
        reference,
    )

    stat = calculate_statistics(
        method_clean,
        reference_clean,
    )

    rows.append({
        "date":
            date.date().isoformat(),
        "year":
            int(date.year),
        "month":
            int(date.month),
        "doy":
            date.timetuple().tm_yday,
        "comparison":
            comparison,
        "method_mean_tecu":
            float(np.mean(method_clean))
            if len(method_clean)
            else np.nan,
        "method_median_tecu":
            float(np.median(method_clean))
            if len(method_clean)
            else np.nan,
        "igs_mean_tecu":
            float(np.mean(reference_clean))
            if len(reference_clean)
            else np.nan,
        "igs_median_tecu":
            float(np.median(reference_clean))
            if len(reference_clean)
            else np.nan,
        **stat,
    })


def sample_for_scatter(df, method_col, reference_col):
    x = df[[method_col, reference_col]].dropna()

    if len(x) > SCATTER_SAMPLE_PER_DAY:
        x = x.sample(
            SCATTER_SAMPLE_PER_DAY,
            random_state=42,
        )

    return x


# ============================================================================
# PLOTS
# ============================================================================

def plot_station_timeseries(modelscale: pd.DataFrame):
    if modelscale.empty:
        return

    path = OUTPUT_ROOT / f"{STATION}_{YEAR}_PyTECGg_VEq_vs_IGS_2H_timeseries.png"

    fig, ax = plt.subplots(figsize=(15, 7))

    ax.plot(
        modelscale["epoch"],
        modelscale["pytecgg_veq"],
        linewidth=1.2,
        label="PyTECGg VEq",
    )

    ax.plot(
        modelscale["epoch"],
        modelscale["igs_vtec"],
        linewidth=1.5,
        label="IGS Final GIM",
    )

    ax.set_title(
        f"{STATION} {YEAR}: PyTECGg VEq vs IGS Final GIM ({MODEL_SCALE} scale)"
    )

    ax.set_xlabel("UTC date")
    ax.set_ylabel("VTEC (TECU)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax.xaxis.set_major_locator(
        mdates.AutoDateLocator()
    )

    ax.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(
            ax.xaxis.get_major_locator()
        )
    )

    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)

    print("Plot:", path)


def plot_station_residuals(modelscale: pd.DataFrame):
    if modelscale.empty:
        return

    path = OUTPUT_ROOT / f"{STATION}_{YEAR}_PyTECGg_minus_IGS_2H_residuals.png"

    resid = (
        modelscale["pytecgg_veq"]
        - modelscale["igs_vtec"]
    )

    fig, ax = plt.subplots(figsize=(15, 6))

    ax.plot(
        modelscale["epoch"],
        resid,
        linewidth=1.0,
    )

    ax.axhline(
        0.0,
        linewidth=1.0,
    )

    ax.set_title(
        f"{STATION} {YEAR}: PyTECGg VEq - IGS Final GIM"
    )

    ax.set_xlabel("UTC date")
    ax.set_ylabel("Residual (TECU)")
    ax.grid(True, alpha=0.3)

    ax.xaxis.set_major_locator(
        mdates.AutoDateLocator()
    )

    ax.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(
            ax.xaxis.get_major_locator()
        )
    )

    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)

    print("Plot:", path)


def plot_scatter(sample_df, method_col, ref_col, title, filename):
    if sample_df.empty:
        return

    x = pd.to_numeric(sample_df[ref_col], errors="coerce")
    y = pd.to_numeric(sample_df[method_col], errors="coerce")

    good = np.isfinite(x) & np.isfinite(y)

    x = np.asarray(x[good], dtype=float)
    y = np.asarray(y[good], dtype=float)

    if len(x) == 0:
        return

    mn = float(min(x.min(), y.min()))
    mx = float(max(x.max(), y.max()))

    pad = max((mx - mn) * 0.05, 1.0)
    mn -= pad
    mx += pad

    path = OUTPUT_ROOT / filename

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(
        x,
        y,
        s=8,
        alpha=0.35,
    )

    ax.plot(
        [mn, mx],
        [mn, mx],
        linewidth=1.2,
        label="1:1",
    )

    ax.set_xlim(mn, mx)
    ax.set_ylim(mn, mx)

    ax.set_xlabel("IGS Final GIM VTEC (TECU)")
    ax.set_ylabel("GNSS-derived VTEC (TECU)")
    ax.set_title(title)

    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)

    print("Plot:", path)


def plot_daily_rmse(daily_stats: pd.DataFrame):
    if daily_stats.empty:
        return

    path = OUTPUT_ROOT / f"{STATION}_{YEAR}_daily_RMSE_vs_IGS.png"

    fig, ax = plt.subplots(figsize=(15, 7))

    for comparison, group in daily_stats.groupby("comparison"):
        x = pd.to_datetime(group["date"])
        y = group["rmse_tecu"]

        ax.plot(
            x,
            y,
            marker="o",
            markersize=3,
            linewidth=1.0,
            label=comparison,
        )

    ax.set_title(
        f"{STATION} {YEAR}: daily TEC RMSE relative to IGS Final GIM"
    )

    ax.set_xlabel("UTC date")
    ax.set_ylabel("RMSE (TECU)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax.xaxis.set_major_locator(
        mdates.AutoDateLocator()
    )

    ax.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(
            ax.xaxis.get_major_locator()
        )
    )

    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)

    print("Plot:", path)


# ============================================================================
# MAIN
# ============================================================================

def validate_current_year():
    """
    Validate the currently configured YEAR.

    Returns
    -------
    summary_df, daily_stats_df, availability_row
    """
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    DAILY_OUTPUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    IONEX_CACHE.mkdir(
        parents=True,
        exist_ok=True,
    )

    print()
    print("=" * 92)
    print(
        f"KOH2 {YEAR} TEC VALIDATION AGAINST IGS FINAL GIM"
    )
    print("=" * 92)
    print("ROOT               :", ROOT)
    print("Station            :", STATION)
    print(f"KOH2 latitude      : {KOH2_LAT:.8f} deg")
    print(f"KOH2 longitude     : {KOH2_LON:.8f} deg")
    print(f"KOH2 height        : {KOH2_H:.3f} m")
    print("Elevation cutoff   :", MIN_ELEVATION_DEG, "deg")
    print("IGS temporal grid  : 2 hours")
    print("=" * 92)

    if not ROOT.is_dir():
        print(
            f"[NOTE] Year root does not exist: {ROOT}"
        )

        availability = {
            "year":
                YEAR,
            "year_root_exists":
                False,
            "pytecgg_days":
                0,
            "pyoasis_days":
                0,
            "union_days":
                0,
            "successful_ionex_days":
                0,
        }

        return (
            pd.DataFrame(),
            pd.DataFrame(),
            availability,
        )

    pytecgg_files, pyoasis_files, days = discover_data_files()

    if SELECTED_DOY is not None:
        pytecgg_files = {
            doy: path
            for doy, path in pytecgg_files.items()
            if doy == SELECTED_DOY
        }
        pyoasis_files = {
            doy: path
            for doy, path in pyoasis_files.items()
            if doy == SELECTED_DOY
        }
        days = [
            doy
            for doy in days
            if doy == SELECTED_DOY
        ]

    print()
    print("PyTECGg days found :", len(pytecgg_files))
    print("pyOASIS days found :", len(pyoasis_files))
    print("Union of days      :", len(days))

    if not days:
        print(
            f"[NOTE] No {YEAR} PyTECGg or pyOASIS TEC files found."
        )

        availability = {
            "year":
                YEAR,
            "year_root_exists":
                True,
            "pytecgg_days":
                len(pytecgg_files),
            "pyoasis_days":
                len(pyoasis_files),
            "union_days":
                0,
            "successful_ionex_days":
                0,
        }

        return (
            pd.DataFrame(),
            pd.DataFrame(),
            availability,
        )

    print()
    print(
        "DOYs:",
        ", ".join(
            f"{d:03d}"
            for d in days
        ),
    )

    accumulators = {
        "PyTECGg_VEq_vs_IGS_station":
            GlobalAccumulator(),
        "PyTECGg_VTEC_vs_IGS_IPP":
            GlobalAccumulator(),
        "pyOASIS_VTEC_vs_IGS_IPP":
            GlobalAccumulator(),
    }

    daily_stats_rows = []
    modelscale_parts = []

    scatter_samples = {
        "PyTECGg_VEq_vs_IGS_station": [],
        "PyTECGg_VTEC_vs_IGS_IPP": [],
        "pyOASIS_VTEC_vs_IGS_IPP": [],
    }

    successful_ionex = 0
    failed_ionex = 0

    for k, doy in enumerate(
        days,
        start=1,
    ):
        date = date_from_year_doy(
            YEAR,
            doy,
        )

        print()
        print("-" * 92)
        print(
            f"{YEAR} DAY {k}/{len(days)}: "
            f"{date:%Y-%m-%d}  DOY {doy:03d}"
        )
        print("-" * 92)

        try:
            ionex_path = download_ionex(
                YEAR,
                doy,
            )

            ionex = read_ionex(
                ionex_path
            )

            print(
                f"  IONEX maps: {len(ionex.epochs)}  "
                f"{ionex.epochs[0]} -> {ionex.epochs[-1]}"
            )

            successful_ionex += 1

        except Exception as exc:
            failed_ionex += 1

            print(
                "[ERROR] IGS Final IONEX failed:",
                exc,
            )

            print(
                "        Skipping this day."
            )

            continue

        # --------------------------------------------------------------------
        # PyTECGg
        # --------------------------------------------------------------------

        pytecgg_path = pytecgg_files.get(
            doy
        )

        if pytecgg_path is not None:
            print()
            print(
                "  PyTECGg:",
                pytecgg_path,
            )

            try:
                pt = load_pytecgg(
                    pytecgg_path
                )

                # A) Station VEq vs IGS at KOH2
                station = (
                    pt[
                        [
                            "epoch",
                            "veq",
                        ]
                    ]
                    .dropna()
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
                    "igs_vtec"
                ] = interpolate_ionex(
                    ionex,
                    station[
                        "epoch"
                    ],
                    np.full(
                        len(station),
                        KOH2_LAT,
                    ),
                    np.full(
                        len(station),
                        KOH2_LON,
                    ),
                )

                station = station.rename(
                    columns={
                        "veq":
                            "pytecgg_veq",
                    }
                )

                station[
                    "residual_tecu"
                ] = (
                    station[
                        "pytecgg_veq"
                    ]
                    - station[
                        "igs_vtec"
                    ]
                )

                station = station.dropna(
                    subset=[
                        "pytecgg_veq",
                        "igs_vtec",
                    ]
                )

                print(
                    "    PyTECGg VEq station matches:",
                    len(station),
                )

                accumulators[
                    "PyTECGg_VEq_vs_IGS_station"
                ].add(
                    station[
                        "pytecgg_veq"
                    ],
                    station[
                        "igs_vtec"
                    ],
                )

                add_daily_stat(
                    daily_stats_rows,
                    date,
                    "PyTECGg_VEq_vs_IGS_station",
                    station[
                        "pytecgg_veq"
                    ],
                    station[
                        "igs_vtec"
                    ],
                )

                save_dataframe(
                    station,
                    DAILY_OUTPUT
                    / (
                        f"{STATION}_{YEAR}{doy:03d}_"
                        "PyTECGg_VEq_vs_IGS_station.parquet"
                    ),
                )

                scatter_samples[
                    "PyTECGg_VEq_vs_IGS_station"
                ].append(
                    sample_for_scatter(
                        station,
                        "pytecgg_veq",
                        "igs_vtec",
                    )
                )

                if not station.empty:
                    tmp = (
                        station
                        .set_index(
                            "epoch"
                        )[
                            [
                                "pytecgg_veq",
                                "igs_vtec",
                            ]
                        ]
                        .resample(
                            MODEL_SCALE
                        )
                        .median()
                        .dropna()
                        .reset_index()
                    )

                    modelscale_parts.append(
                        tmp
                    )

                # B) PyTECGg VTEC at IPP vs IGS at same IPP
                ipp = (
                    pt.loc[
                        pt[
                            "ele"
                        ]
                        >= MIN_ELEVATION_DEG,
                        [
                            "epoch",
                            "sv",
                            "lat_ipp",
                            "lon_ipp",
                            "ele",
                            "vtec",
                        ],
                    ]
                    .dropna(
                        subset=[
                            "epoch",
                            "lat_ipp",
                            "lon_ipp",
                            "vtec",
                        ]
                    )
                    .copy()
                )

                ipp[
                    "igs_vtec"
                ] = interpolate_ionex(
                    ionex,
                    ipp[
                        "epoch"
                    ],
                    ipp[
                        "lat_ipp"
                    ],
                    ipp[
                        "lon_ipp"
                    ],
                )

                ipp[
                    "residual_tecu"
                ] = (
                    ipp[
                        "vtec"
                    ]
                    - ipp[
                        "igs_vtec"
                    ]
                )

                ipp = ipp.dropna(
                    subset=[
                        "vtec",
                        "igs_vtec",
                    ]
                )

                print(
                    "    PyTECGg IPP matches        :",
                    len(ipp),
                )

                accumulators[
                    "PyTECGg_VTEC_vs_IGS_IPP"
                ].add(
                    ipp[
                        "vtec"
                    ],
                    ipp[
                        "igs_vtec"
                    ],
                )

                add_daily_stat(
                    daily_stats_rows,
                    date,
                    "PyTECGg_VTEC_vs_IGS_IPP",
                    ipp[
                        "vtec"
                    ],
                    ipp[
                        "igs_vtec"
                    ],
                )

                save_dataframe(
                    ipp,
                    DAILY_OUTPUT
                    / (
                        f"{STATION}_{YEAR}{doy:03d}_"
                        "PyTECGg_VTEC_vs_IGS_IPP.parquet"
                    ),
                )

                scatter_samples[
                    "PyTECGg_VTEC_vs_IGS_IPP"
                ].append(
                    sample_for_scatter(
                        ipp,
                        "vtec",
                        "igs_vtec",
                    )
                )

            except Exception as exc:
                print(
                    "    [ERROR] PyTECGg validation failed:",
                    exc,
                )

        else:
            print(
                "  PyTECGg: no file for this day."
            )

        # --------------------------------------------------------------------
        # pyOASIS
        # --------------------------------------------------------------------

        pyoasis_path = pyoasis_files.get(
            doy
        )

        if pyoasis_path is not None:
            print()
            print(
                "  pyOASIS:",
                pyoasis_path,
            )

            try:
                po = load_pyoasis(
                    pyoasis_path
                )

                po = (
                    po.loc[
                        po[
                            "elevation"
                        ].abs()
                        >= MIN_ELEVATION_DEG,
                        [
                            "epoch",
                            "sat",
                            "lat_ipp",
                            "lon_ipp",
                            "elevation",
                            "vtec",
                        ],
                    ]
                    .dropna(
                        subset=[
                            "epoch",
                            "lat_ipp",
                            "lon_ipp",
                            "vtec",
                        ]
                    )
                    .copy()
                )

                po[
                    "igs_vtec"
                ] = interpolate_ionex(
                    ionex,
                    po[
                        "epoch"
                    ],
                    po[
                        "lat_ipp"
                    ],
                    po[
                        "lon_ipp"
                    ],
                )

                po[
                    "residual_tecu"
                ] = (
                    po[
                        "vtec"
                    ]
                    - po[
                        "igs_vtec"
                    ]
                )

                po = po.dropna(
                    subset=[
                        "vtec",
                        "igs_vtec",
                    ]
                )

                print(
                    "    pyOASIS IPP matches         :",
                    len(po),
                )

                accumulators[
                    "pyOASIS_VTEC_vs_IGS_IPP"
                ].add(
                    po[
                        "vtec"
                    ],
                    po[
                        "igs_vtec"
                    ],
                )

                add_daily_stat(
                    daily_stats_rows,
                    date,
                    "pyOASIS_VTEC_vs_IGS_IPP",
                    po[
                        "vtec"
                    ],
                    po[
                        "igs_vtec"
                    ],
                )

                save_dataframe(
                    po,
                    DAILY_OUTPUT
                    / (
                        f"{STATION}_{YEAR}{doy:03d}_"
                        "pyOASIS_VTEC_vs_IGS_IPP.parquet"
                    ),
                )

                scatter_samples[
                    "pyOASIS_VTEC_vs_IGS_IPP"
                ].append(
                    sample_for_scatter(
                        po,
                        "vtec",
                        "igs_vtec",
                    )
                )

            except Exception as exc:
                print(
                    "    [ERROR] pyOASIS validation failed:",
                    exc,
                )

        else:
            print(
                "  pyOASIS: no L1-L2 TEC file for this day."
            )

    # ------------------------------------------------------------------------
    # SUMMARY TABLES
    # ------------------------------------------------------------------------

    summary_rows = []

    for comparison, acc in accumulators.items():
        summary_rows.append({
            "year":
                YEAR,
            "comparison":
                comparison,
            **acc.statistics(),
        })

    summary_df = pd.DataFrame(
        summary_rows
    )

    daily_stats_df = pd.DataFrame(
        daily_stats_rows
    )

    summary_csv = (
        OUTPUT_ROOT
        / f"{STATION}_{YEAR}_IGS_validation_summary.csv"
    )

    daily_csv = (
        OUTPUT_ROOT
        / f"{STATION}_{YEAR}_IGS_validation_daily_statistics.csv"
    )

    summary_df.to_csv(
        summary_csv,
        index=False,
    )

    daily_stats_df.to_csv(
        daily_csv,
        index=False,
    )

    # Model-scale station data.
    if modelscale_parts:
        modelscale = (
            pd.concat(
                modelscale_parts,
                ignore_index=True,
            )
            .drop_duplicates(
                subset=[
                    "epoch",
                ]
            )
            .sort_values(
                "epoch"
            )
            .reset_index(
                drop=True
            )
        )

        modelscale[
            "residual_tecu"
        ] = (
            modelscale[
                "pytecgg_veq"
            ]
            - modelscale[
                "igs_vtec"
            ]
        )

        modelscale_csv = (
            OUTPUT_ROOT
            / f"{STATION}_{YEAR}_PyTECGg_VEq_vs_IGS_2H.csv"
        )

        modelscale.to_csv(
            modelscale_csv,
            index=False,
        )

    else:
        modelscale = pd.DataFrame()

    # ------------------------------------------------------------------------
    # PLOTS
    # ------------------------------------------------------------------------

    plot_station_timeseries(
        modelscale
    )

    plot_station_residuals(
        modelscale
    )

    plot_daily_rmse(
        daily_stats_df
    )

    if scatter_samples[
        "PyTECGg_VEq_vs_IGS_station"
    ]:
        x = pd.concat(
            scatter_samples[
                "PyTECGg_VEq_vs_IGS_station"
            ],
            ignore_index=True,
        )

        plot_scatter(
            x,
            "pytecgg_veq",
            "igs_vtec",
            (
                f"{STATION} {YEAR}: "
                "PyTECGg VEq vs IGS Final GIM"
            ),
            (
                f"{STATION}_{YEAR}_"
                "PyTECGg_VEq_vs_IGS_scatter.png"
            ),
        )

    if scatter_samples[
        "PyTECGg_VTEC_vs_IGS_IPP"
    ]:
        x = pd.concat(
            scatter_samples[
                "PyTECGg_VTEC_vs_IGS_IPP"
            ],
            ignore_index=True,
        )

        plot_scatter(
            x,
            "vtec",
            "igs_vtec",
            (
                f"{STATION} {YEAR}: "
                "PyTECGg IPP VTEC vs IGS Final GIM"
            ),
            (
                f"{STATION}_{YEAR}_"
                "PyTECGg_IPP_VTEC_vs_IGS_scatter.png"
            ),
        )

    if scatter_samples[
        "pyOASIS_VTEC_vs_IGS_IPP"
    ]:
        x = pd.concat(
            scatter_samples[
                "pyOASIS_VTEC_vs_IGS_IPP"
            ],
            ignore_index=True,
        )

        plot_scatter(
            x,
            "vtec",
            "igs_vtec",
            (
                f"{STATION} {YEAR}: "
                "pyOASIS IPP VTEC vs IGS Final GIM"
            ),
            (
                f"{STATION}_{YEAR}_"
                "pyOASIS_IPP_VTEC_vs_IGS_scatter.png"
            ),
        )

    # ------------------------------------------------------------------------
    # REPORT
    # ------------------------------------------------------------------------

    report_path = (
        OUTPUT_ROOT
        / f"{STATION}_{YEAR}_IGS_validation_report.txt"
    )

    with open(
        report_path,
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            "KOH2 TEC validation against IGS Final GIM\n"
        )

        f.write(
            "=" * 72
            + "\n"
        )

        f.write(
            f"Year: {YEAR}\n"
        )

        f.write(
            f"Discovered days: {len(days)}\n"
        )

        f.write(
            f"PyTECGg days: {len(pytecgg_files)}\n"
        )

        f.write(
            f"pyOASIS days: {len(pyoasis_files)}\n"
        )

        f.write(
            f"Successful IONEX days: {successful_ionex}\n"
        )

        f.write(
            f"Failed/unavailable IONEX days: {failed_ionex}\n"
        )

        f.write(
            f"Elevation cutoff: {MIN_ELEVATION_DEG:.1f} deg\n"
        )

        f.write(
            f"KOH2: lat={KOH2_LAT:.8f}, "
            f"lon={KOH2_LON:.8f}, "
            f"h={KOH2_H:.3f} m\n"
        )

        f.write(
            "\n"
        )

        f.write(
            "Residual convention: GNSS-derived VTEC minus IGS GIM VTEC\n"
        )

        f.write(
            "\n"
        )

        f.write(
            summary_df.to_string(
                index=False
            )
        )

        f.write(
            "\n"
        )

    print()
    print(
        f"{YEAR} SUMMARY"
    )

    print(
        "-" * 92
    )

    print(
        summary_df.to_string(
            index=False
        )
    )

    print(
        "-" * 92
    )

    availability = {
        "year":
            YEAR,
        "year_root_exists":
            True,
        "pytecgg_days":
            len(pytecgg_files),
        "pyoasis_days":
            len(pyoasis_files),
        "union_days":
            len(days),
        "successful_ionex_days":
            successful_ionex,
        "failed_or_unavailable_ionex_days":
            failed_ionex,
        "output":
            str(
                OUTPUT_ROOT
            ),
    }

    return (
        summary_df,
        daily_stats_df,
        availability,
    )


def create_equal_day_summaries(
    daily_all: pd.DataFrame,
):
    if daily_all.empty:
        return (
            pd.DataFrame(),
            pd.DataFrame(),
        )

    rows = []

    for (
        year,
        comparison,
    ), g in daily_all.groupby(
        [
            "year",
            "comparison",
        ],
        dropna=False,
    ):
        rows.append({
            "year":
                int(year),
            "comparison":
                comparison,
            "n_days":
                int(
                    g[
                        "date"
                    ].nunique()
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
            "mean_daily_bias_tecu":
                g[
                    "bias_tecu"
                ].mean(),
            "median_daily_bias_tecu":
                g[
                    "bias_tecu"
                ].median(),
            "mean_daily_mae_tecu":
                g[
                    "mae_tecu"
                ].mean(),
            "median_daily_mae_tecu":
                g[
                    "mae_tecu"
                ].median(),
            "mean_daily_rmse_tecu":
                g[
                    "rmse_tecu"
                ].mean(),
            "median_daily_rmse_tecu":
                g[
                    "rmse_tecu"
                ].median(),
            "median_daily_residual_std_tecu":
                g[
                    "std_residual_tecu"
                ].median(),
            "median_daily_pearson_r":
                g[
                    "pearson_r"
                ].median(),
            "median_of_daily_method_median_tecu":
                g[
                    "method_median_tecu"
                ].median(),
            "median_of_daily_igs_median_tecu":
                g[
                    "igs_median_tecu"
                ].median(),
        })

    annual = pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "comparison",
        ]
    )

    month_rows = []

    for (
        year,
        month,
        comparison,
    ), g in daily_all.groupby(
        [
            "year",
            "month",
            "comparison",
        ],
        dropna=False,
    ):
        month_rows.append({
            "year":
                int(year),
            "month":
                int(month),
            "comparison":
                comparison,
            "n_days":
                int(
                    g[
                        "date"
                    ].nunique()
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
            "median_daily_bias_tecu":
                g[
                    "bias_tecu"
                ].median(),
            "median_daily_rmse_tecu":
                g[
                    "rmse_tecu"
                ].median(),
            "median_daily_pearson_r":
                g[
                    "pearson_r"
                ].median(),
            "median_of_daily_method_median_tecu":
                g[
                    "method_median_tecu"
                ].median(),
            "median_of_daily_igs_median_tecu":
                g[
                    "igs_median_tecu"
                ].median(),
        })

    monthly = pd.DataFrame(
        month_rows
    ).sort_values(
        [
            "year",
            "month",
            "comparison",
        ]
    )

    return (
        annual,
        monthly,
    )


def plot_multiyear_metric(
    yearly_summary: pd.DataFrame,
    metric: str,
    ylabel: str,
    filename: str,
):
    if yearly_summary.empty:
        return

    fig, ax = plt.subplots(
        figsize=(
            11,
            6,
        )
    )

    for comparison, g in yearly_summary.groupby(
        "comparison"
    ):
        g = g.sort_values(
            "year"
        )

        ax.plot(
            g[
                "year"
            ],
            g[
                metric
            ],
            marker="o",
            linewidth=1.2,
            label=comparison,
        )

    if metric == "bias_tecu":
        ax.axhline(
            0.0,
            linewidth=1.0,
        )

    ax.set_xlabel(
        "Year"
    )

    ax.set_ylabel(
        ylabel
    )

    ax.set_title(
        f"{STATION} 2019-2026: {ylabel} relative to IGS Final GIM"
    )

    ax.grid(
        True,
        alpha=0.3,
    )

    ax.legend(
        fontsize=8,
    )

    fig.tight_layout()

    fig.savefig(
        MASTER_OUTPUT
        / filename,
        dpi=DPI,
    )

    plt.close(
        fig
    )


def write_master_report(
    yearly_summary,
    equal_day,
    availability,
):
    path = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "IGS_validation_master_report.txt"
        )
    )

    with open(
        path,
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            "KOH2 2019-2026 TEC VALIDATION AGAINST IGS FINAL GIM\n"
        )

        f.write(
            "=" * 90
            + "\n\n"
        )

        f.write(
            "Residual convention: method - IGS\n"
        )

        f.write(
            f"Elevation cutoff for IPP comparisons: "
            f"{MIN_ELEVATION_DEG:.1f} deg\n"
        )

        f.write(
            "Only available KOH2 observation days are processed; "
            "the script does not assume continuous yearly coverage.\n"
        )

        f.write(
            "IGS Final GIM is used consistently. Recent dates for which a "
            "Final GIM is not yet available are reported as unavailable rather "
            "than mixing Rapid and Final products.\n\n"
        )

        f.write(
            "DATA AVAILABILITY\n"
        )

        f.write(
            "-" * 90
            + "\n"
        )

        f.write(
            availability.to_string(
                index=False
            )
        )

        f.write(
            "\n\nPOINT-WEIGHTED YEARLY VALIDATION\n"
        )

        f.write(
            "-" * 90
            + "\n"
        )

        f.write(
            yearly_summary.to_string(
                index=False
            )
        )

        f.write(
            "\n\nEQUAL-DAY YEARLY SUMMARY\n"
        )

        f.write(
            "-" * 90
            + "\n"
        )

        f.write(
            equal_day.to_string(
                index=False
            )
        )

        f.write(
            "\n"
        )

    return path


def main():
    args = parse_args()
    configure_runtime(args)

    MASTER_OUTPUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    CACHE_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    BY_YEAR_OUTPUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "=" * 100
    )

    print(
        "KOH2 2019-2026 MULTI-YEAR TEC VALIDATION AGAINST IGS FINAL GIM"
    )

    print(
        "=" * 100
    )

    print(
        "Python:",
        sys.executable,
    )

    print(
        "Years:",
        ", ".join(
            str(y)
            for y in YEARS
        ),
    )

    print(
        "Base root:",
        BASE_ROOT,
    )

    print(
        "Master output:",
        MASTER_OUTPUT,
    )

    print(
        "IONEX cache root:",
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
        "=" * 100
    )

    all_year_summaries = []
    all_daily = []
    availability_rows = []

    for year in YEARS:
        configure_year(
            year
        )

        try:
            summary_df, daily_df, availability = validate_current_year()

        except Exception as exc:
            print()
            print(
                f"[FATAL YEAR ERROR] {year}:",
                repr(
                    exc
                ),
            )

            summary_df = pd.DataFrame()
            daily_df = pd.DataFrame()

            availability = {
                "year":
                    year,
                "year_root_exists":
                    ROOT.is_dir(),
                "pytecgg_days":
                    np.nan,
                "pyoasis_days":
                    np.nan,
                "union_days":
                    np.nan,
                "successful_ionex_days":
                    np.nan,
                "status":
                    f"ERROR: {exc!r}",
            }

        availability_rows.append(
            availability
        )

        if not summary_df.empty:
            all_year_summaries.append(
                summary_df
            )

        if not daily_df.empty:
            all_daily.append(
                daily_df
            )

    yearly_summary = (
        pd.concat(
            all_year_summaries,
            ignore_index=True,
        )
        if all_year_summaries
        else pd.DataFrame()
    )

    daily_all = (
        pd.concat(
            all_daily,
            ignore_index=True,
        )
        if all_daily
        else pd.DataFrame()
    )

    availability_df = pd.DataFrame(
        availability_rows
    )

    (
        equal_day_summary,
        monthly_summary,
    ) = create_equal_day_summaries(
        daily_all
    )

    yearly_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "IGS_validation_yearly_summary.csv"
        )
    )

    daily_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "IGS_validation_daily_statistics.csv"
        )
    )

    equal_day_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "IGS_validation_equal_day_yearly_summary.csv"
        )
    )

    monthly_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "IGS_validation_monthly_summary.csv"
        )
    )

    availability_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "validation_data_availability.csv"
        )
    )

    yearly_summary.to_csv(
        yearly_file,
        index=False,
    )

    daily_all.to_csv(
        daily_file,
        index=False,
    )

    equal_day_summary.to_csv(
        equal_day_file,
        index=False,
    )

    monthly_summary.to_csv(
        monthly_file,
        index=False,
    )

    availability_df.to_csv(
        availability_file,
        index=False,
    )

    if not yearly_summary.empty:
        plot_multiyear_metric(
            yearly_summary,
            "bias_tecu",
            "Bias (TECU)",
            (
                f"{STATION}_2019_2026_"
                "annual_bias_vs_IGS.png"
            ),
        )

        plot_multiyear_metric(
            yearly_summary,
            "rmse_tecu",
            "RMSE (TECU)",
            (
                f"{STATION}_2019_2026_"
                "annual_RMSE_vs_IGS.png"
            ),
        )

        plot_multiyear_metric(
            yearly_summary,
            "pearson_r",
            "Pearson correlation r",
            (
                f"{STATION}_2019_2026_"
                "annual_correlation_vs_IGS.png"
            ),
        )

    report_file = write_master_report(
        yearly_summary,
        equal_day_summary,
        availability_df,
    )

    print()
    print(
        "=" * 100
    )

    print(
        "MULTI-YEAR VALIDATION COMPLETE"
    )

    print(
        "=" * 100
    )

    print(
        "Availability:",
        availability_file,
    )

    print(
        "Yearly summary:",
        yearly_file,
    )

    print(
        "Equal-day summary:",
        equal_day_file,
    )

    print(
        "Monthly summary:",
        monthly_file,
    )

    print(
        "Daily statistics:",
        daily_file,
    )

    print(
        "Report:",
        report_file,
    )

    print()
    print(
        "YEARLY POINT-WEIGHTED SUMMARY"
    )

    print(
        "-" * 100
    )

    if yearly_summary.empty:
        print(
            "No validation results."
        )
    else:
        print(
            yearly_summary.to_string(
                index=False
            )
        )

    print(
        "-" * 100
    )


if __name__ == "__main__":
    main()
