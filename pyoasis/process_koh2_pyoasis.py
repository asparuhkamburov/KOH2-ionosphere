#!/usr/bin/env python3
from __future__ import annotations

r"""
KOH2 pyOASIS processing workflow
================================

Publication-oriented version of the processing script used for KOH2.

The program:
1. Traverses a year/month/day KOH2 archive.
2. Finds a suitable observation RINEX file for each day.
3. Downloads and extracts the corresponding GFZ MGEX rapid SP3 orbit.
4. Prepares a short-name RINEX 2.x staging file for pyOASIS using GFZRNX.
5. Runs the pyOASIS processing sequence:
       SP3intp -> RNXclean -> RNXlevelling
       -> ROTIcalc -> DTECcalc -> SIDXcalc -> TECcalc
6. Skips days whose final TEC product already exists unless --force is used.

The script contains no machine-specific absolute paths. Paths and year are
provided on the command line or through environment variables.

Example
-------
python process_koh2_pyoasis.py ^
    --year 2025 ^
    --date 2025-01-01 ^
    --data-root "<PYOASIS_DATA_ROOT>" ^
    --gfzrnx "D:\GNSS\gfzrnx_2.2.0_win10_64.exe"

The --data-root argument may point either to:
    <PYOASIS_DATA_ROOT>
or directly to:
    <PYOASIS_DATA_ROOT>\2023

Notes
-----
- CDDIS/Earthdata authentication is expected to be configured outside this
  script (for example through the user's .netrc/_netrc configuration).
- GFZRNX is an external executable and is not distributed with this repository.
- pyOASIS is an external dependency and is not redistributed here.
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import requests

# GeoRinex/xarray may emit FutureWarnings during pyOASIS processing.
warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    module=r"georinex(\..*)?$",
)


DEFAULT_CDDIS_BASE = "https://cddis.nasa.gov/archive/gnss/products"


@dataclass(frozen=True)
class Config:
    year: int
    year_root: Path
    gfzrnx: Path
    station: str
    cddis_base: str
    keep_gz: bool
    timeout: int
    force: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Process one year of KOH2 GNSS observations with GFZRNX and pyOASIS."
        )
    )

    parser.add_argument(
        "--year",
        type=int,
        required=True,
        help="Calendar year to process, e.g. 2025.",
    )

    parser.add_argument(
        "--date",
        type=str,
        default=None,
        help=(
            "Optional single UTC calendar date to process in YYYY-MM-DD format, "
            "e.g. 2025-01-01. The date must belong to --year."
        ),
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=os.environ.get("KOH2_DATA_ROOT"),
        help=(
            "Root containing annual folders, or the annual folder itself. "
            "Can also be set with KOH2_DATA_ROOT."
        ),
    )

    parser.add_argument(
        "--gfzrnx",
        type=Path,
        default=os.environ.get("GFZRNX_PATH"),
        help=(
            "Path to the GFZRNX executable. "
            "Can also be set with GFZRNX_PATH."
        ),
    )

    parser.add_argument(
        "--station",
        default="KOH2",
        help="Four-character station code. Default: KOH2.",
    )

    parser.add_argument(
        "--cddis-base",
        default=DEFAULT_CDDIS_BASE,
        help=f"CDDIS products base URL. Default: {DEFAULT_CDDIS_BASE}",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=90,
        help="HTTP download timeout in seconds. Default: 90.",
    )

    parser.add_argument(
        "--delete-gz",
        action="store_true",
        help="Delete downloaded .gz orbit files after successful extraction.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Reprocess days even when the expected final TEC file exists.",
    )

    return parser.parse_args()


def resolve_year_root(data_root: Path, year: int) -> Path:
    """Accept either a parent data root or the year directory itself."""
    data_root = data_root.expanduser().resolve()

    if data_root.is_dir() and data_root.name == str(year):
        return data_root

    year_root = data_root / str(year)

    if year_root.is_dir():
        return year_root

    raise FileNotFoundError(
        "Could not find the requested year directory.\n"
        f"Data root supplied: {data_root}\n"
        f"Tried: {data_root}\n"
        f"Tried: {year_root}"
    )


def parse_selected_date(
    date_text: str | None,
    year: int,
) -> datetime | None:
    """Parse an optional single-date restriction."""
    if date_text is None:
        return None

    try:
        selected = datetime.strptime(
            date_text,
            "%Y-%m-%d",
        )
    except ValueError as exc:
        raise ValueError(
            "--date must use YYYY-MM-DD format, e.g. 2025-01-01."
        ) from exc

    if selected.year != year:
        raise ValueError(
            f"--date {date_text} does not belong to --year {year}."
        )

    return selected


def validate_inputs(args: argparse.Namespace) -> Config:
    if args.data_root is None:
        raise SystemExit(
            "Missing --data-root. Supply it on the command line or set "
            "KOH2_DATA_ROOT."
        )

    if args.gfzrnx is None:
        raise SystemExit(
            "Missing --gfzrnx. Supply it on the command line or set "
            "GFZRNX_PATH."
        )

    gfzrnx = Path(args.gfzrnx).expanduser().resolve()

    if not gfzrnx.is_file():
        raise FileNotFoundError(
            f"GFZRNX executable not found:\n    {gfzrnx}"
        )

    year_root = resolve_year_root(
        Path(args.data_root),
        args.year,
    )

    station = args.station.strip().upper()

    if len(station) != 4:
        raise ValueError(
            f"Station code must contain four characters: {station!r}"
        )

    return Config(
        year=args.year,
        year_root=year_root,
        gfzrnx=gfzrnx,
        station=station,
        cddis_base=args.cddis_base.rstrip("/"),
        keep_gz=not args.delete_gz,
        timeout=args.timeout,
        force=args.force,
    )


def get_rinex_version(filepath: Path) -> str | None:
    """Return the RINEX version written in the first header line."""
    try:
        with filepath.open(
            "r",
            encoding="ascii",
            errors="ignore",
        ) as stream:
            first_line = stream.readline()

        return first_line[:9].strip()

    except OSError:
        return None


def gps_week(date: datetime) -> int:
    """Return the GPS week number for a calendar date."""
    gps_epoch = datetime(1980, 1, 6)
    return (date - gps_epoch).days // 7


def extract_gzip(gz_path: Path, keep_gz: bool = True) -> Path:
    """Extract a .gz file and return the uncompressed path."""
    output_path = gz_path.with_suffix("")

    if output_path.is_file():
        print("SP3 already extracted:")
        print(" ", output_path)
        return output_path

    print("Extracting:")
    print(" ", gz_path.name)

    with gzip.open(gz_path, "rb") as source:
        with output_path.open("wb") as target:
            shutil.copyfileobj(source, target)

    if not keep_gz:
        gz_path.unlink(missing_ok=True)

    print("Extracted SP3:")
    print(" ", output_path)

    return output_path


def download_sp3(
    date: datetime,
    sp3_folder: Path,
    config: Config,
) -> Path | None:
    """
    Download the GFZ MGEX rapid precise orbit used by this workflow.

    Expected filename:
        GFZ0MGXRAP_YYYYDDD0000_01D_05M_ORB.SP3.gz
    """
    year = date.year
    doy = int(date.strftime("%j"))
    week = gps_week(date)

    filename = (
        f"GFZ0MGXRAP_"
        f"{year}{doy:03d}0000_"
        f"01D_05M_ORB.SP3.gz"
    )

    sp3_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    gz_path = sp3_folder / filename
    sp3_path = gz_path.with_suffix("")

    print()
    print("GPS WEEK:", week)
    print("Expected SP3:")
    print(" ", filename)

    if sp3_path.is_file():
        print("SP3 already exists:")
        print(" ", sp3_path)
        return sp3_path

    if not gz_path.is_file():
        url = f"{config.cddis_base}/{week:04d}/{filename}"

        print("URL:")
        print(" ", url)

        temp_path = gz_path.with_suffix(
            gz_path.suffix + ".part"
        )

        try:
            with requests.get(
                url,
                stream=True,
                timeout=config.timeout,
                allow_redirects=True,
            ) as response:

                if response.status_code == 404:
                    print()
                    print("SP3 FILE NOT FOUND:")
                    print(" ", filename)
                    return None

                final_url = response.url.lower()
                content_type = response.headers.get(
                    "Content-Type",
                    "",
                ).lower()

                if (
                    "urs.earthdata.nasa.gov" in final_url
                    or "text/html" in content_type
                ):
                    print()
                    print("EARTHDATA LOGIN ERROR")
                    print(
                        "Check .netrc/_netrc and Earthdata authorization."
                    )
                    return None

                response.raise_for_status()

                with temp_path.open("wb") as stream:
                    for chunk in response.iter_content(
                        chunk_size=1024 * 1024
                    ):
                        if chunk:
                            stream.write(chunk)

                os.replace(
                    temp_path,
                    gz_path,
                )

                print("Downloaded:")
                print(" ", gz_path)

        except requests.RequestException as exc:
            temp_path.unlink(missing_ok=True)
            print()
            print("DOWNLOAD ERROR:")
            print(" ", exc)
            return None

    else:
        print("Compressed SP3 already exists:")
        print(" ", gz_path)

    try:
        return extract_gzip(
            gz_path,
            keep_gz=config.keep_gz,
        )

    except (OSError, gzip.BadGzipFile) as exc:
        print()
        print("SP3 EXTRACTION FAILED:")
        print(" ", exc)
        return None


def iter_day_folders(
    year_root: Path,
    year: int,
):
    """Yield valid (date, day_path) pairs from YEAR/MONTH/DAY folders."""
    for month_path in sorted(year_root.iterdir()):
        if not month_path.is_dir():
            continue

        if not month_path.name.isdigit():
            continue

        for day_path in sorted(month_path.iterdir()):
            if not day_path.is_dir():
                continue

            if not day_path.name.isdigit():
                continue

            try:
                date = datetime(
                    year,
                    int(month_path.name),
                    int(day_path.name),
                )

            except ValueError:
                print(
                    "Invalid date folder:",
                    day_path,
                )
                continue

            yield date, day_path


def find_observation_files(
    rinex_folder: Path,
) -> list[Path]:
    """Return supported observation files from a day's RINEX folder."""
    if not rinex_folder.is_dir():
        return []

    supported = []

    for path in rinex_folder.iterdir():
        if not path.is_file():
            continue

        name = path.name.lower()

        if (
            name.endswith(".rnx")
            or name.endswith(".crx")
            or name.endswith(".rnx.gz")
            or name.endswith(".crx.gz")
            or name.endswith("o")
        ):
            supported.append(path)

    return sorted(supported)


def select_observation_file(
    observation_files: list[Path],
) -> Path | None:
    """
    Reproduce the original pyOASIS input priority:

    1. Full-day 30-second mixed RINEX.
    2. Any available 30-second mixed RINEX segment.
    3. Full-day 15-second mixed RINEX.
    """
    rules = [
        "_01D_30S_MO",
        "_30S_MO",
        "_01D_15S_MO",
    ]

    for token in rules:
        for path in observation_files:
            if token in path.name.upper():
                return path

    return None


def short_rinex_name(
    station: str,
    date: datetime,
) -> str:
    """Return the short RINEX observation name expected by pyOASIS."""
    doy = int(date.strftime("%j"))
    yy = date.strftime("%y")

    return (
        f"{station.lower()}"
        f"{doy:03d}"
        f"0."
        f"{yy}"
        f"o"
    )


def prepare_rinex2(
    selected_obs_path: Path,
    staging_path: Path,
    gfzrnx: Path,
) -> Path:
    """
    Create or reuse the RINEX 2.x short-name staging file required by pyOASIS.
    """
    staging_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if staging_path.is_file():
        existing_version = get_rinex_version(
            staging_path
        )

        print()
        print(
            "Existing pyOASIS staging RINEX version:",
            existing_version,
        )

        if (
            existing_version is not None
            and existing_version.startswith("2.")
        ):
            print(
                "Existing RINEX 2 file is suitable."
            )
            return staging_path

        print(
            "Existing staging file is not RINEX 2; removing it."
        )
        staging_path.unlink()

    print()
    print("----------------------------------------")
    print("GFZRNX: RINEX -> RINEX 2.x")
    print("----------------------------------------")
    print("INPUT:")
    print(" ", selected_obs_path)
    print("OUTPUT:")
    print(" ", staging_path)

    command = [
        str(gfzrnx),
        "-finp",
        str(selected_obs_path),
        "-fout",
        str(staging_path),
        "-vosc",
        "2",
    ]

    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )

        if result.stdout:
            print(result.stdout)

        if result.stderr:
            print(result.stderr)

    except subprocess.CalledProcessError as exc:
        print()
        print("GFZRNX CONVERSION FAILED")
        print("Return code:", exc.returncode)
        print("STDOUT:")
        print(exc.stdout)
        print("STDERR:")
        print(exc.stderr)
        raise

    converted_version = get_rinex_version(
        staging_path
    )

    print()
    print(
        "pyOASIS RINEX version:",
        converted_version,
    )

    if not (
        converted_version
        and converted_version.startswith("2.")
    ):
        raise RuntimeError(
            "GFZRNX output is not RINEX 2.x"
        )

    print(
        "RINEX 2 preparation successful:"
    )
    print(" ", staging_path)

    return staging_path


def expected_tec_path(
    day_path: Path,
    station: str,
    date: datetime,
) -> Path:
    doy = int(date.strftime("%j"))

    return (
        day_path
        / "pyOASIS_OUTPUT"
        / "INDICES"
        / "TEC"
        / f"{station}_{doy:03d}_{date.year}_L1L2.TEC"
    )


def run_pyoasis_day(
    date: datetime,
    day_path: Path,
    selected_obs_path: Path,
    sp3_path: Path,
    config: Config,
) -> None:
    """Execute the complete pyOASIS processing chain for one day."""
    try:
        import pyOASIS
    except ImportError as exc:
        raise RuntimeError(
            "pyOASIS could not be imported. Activate the Python environment "
            "in which pyOASIS is installed."
        ) from exc
    doy = int(date.strftime("%j"))
    doy_str = f"{doy:03d}"
    year_str = str(date.year)
    station = config.station

    output_folder = (
        day_path
        / "pyOASIS_OUTPUT"
    )

    oasis_input_rinex = (
        day_path
        / "pyOASIS_INPUT"
        / "RINEX"
    )

    oasis_orbit_output = (
        output_folder
        / "ORBITS"
    )

    oasis_rinex_output = (
        output_folder
        / "RINEX"
    )

    sp3_folder = (
        day_path
        / "PRODUCTS"
        / "SP3"
    )

    for folder in [
        output_folder,
        oasis_input_rinex,
        oasis_orbit_output,
        oasis_rinex_output,
    ]:
        folder.mkdir(
            parents=True,
            exist_ok=True,
        )

    staging_path = (
        oasis_input_rinex
        / short_rinex_name(
            station,
            date,
        )
    )

    prepare_rinex2(
        selected_obs_path,
        staging_path,
        config.gfzrnx,
    )

    print()
    print("========================================")
    print("STARTING pyOASIS")
    print("========================================")
    print("DATE :", date.strftime("%Y-%m-%d"))
    print("OBS  :", selected_obs_path)
    print("SP3  :", sp3_path)

    # STEP 1 — SP3 interpolation
    print()
    print("----------------------------------------")
    print("pyOASIS STEP 1: SP3intp")
    print("----------------------------------------")

    pyOASIS.SP3intp(
        year_str,
        doy_str,
        str(sp3_folder),
        str(oasis_orbit_output),
    )

    orbit_processed = (
        oasis_orbit_output
        / f"ORBITS_{year_str}_{doy_str}.SP3"
    )

    if not orbit_processed.is_file():
        raise RuntimeError(
            "pyOASIS did not generate the expected interpolated orbit:\n"
            f"    {orbit_processed}"
        )

    print(
        "Orbit preparation successful:"
    )
    print(" ", orbit_processed)

    # STEP 2 — RINEX cleaning
    print()
    print("----------------------------------------")
    print("pyOASIS STEP 2: RNXclean")
    print("----------------------------------------")

    pyOASIS.RNXclean(
        station,
        doy_str,
        year_str,
        str(oasis_input_rinex),
        str(oasis_orbit_output),
        str(oasis_rinex_output),
    )

    # STEP 3 — RINEX levelling
    print()
    print("----------------------------------------")
    print("pyOASIS STEP 3: RNXlevelling")
    print("----------------------------------------")

    original_stderr = sys.stderr

    try:
        pyOASIS.RNXlevelling(
            station,
            str(oasis_rinex_output),
            show_plot=False,
        )

    finally:
        # RNXlevelling may redirect sys.stderr internally.
        sys.stderr = original_stderr

    rnx3_files = [
        path
        for path in oasis_rinex_output.iterdir()
        if path.is_file()
        and path.name.upper().endswith(".RNX3")
    ]

    print()
    print(
        "RNX3 FILES GENERATED:",
        len(rnx3_files),
    )

    if not rnx3_files:
        raise RuntimeError(
            "RNXlevelling did not generate any RNX3 files."
        )

    # Index directories
    indices_root = (
        output_folder
        / "INDICES"
    )

    roti_output = indices_root / "ROTI"
    dtec_output = indices_root / "DTEC"
    sidx_output = indices_root / "SIDX"
    tec_output = indices_root / "TEC"

    for folder in [
        roti_output,
        dtec_output,
        sidx_output,
        tec_output,
    ]:
        folder.mkdir(
            parents=True,
            exist_ok=True,
        )

    # STEP 4 — ROTI
    print()
    print("----------------------------------------")
    print("pyOASIS STEP 4: ROTIcalc")
    print("----------------------------------------")

    pyOASIS.ROTIcalc(
        station,
        doy_str,
        year_str,
        str(oasis_rinex_output),
        str(roti_output),
        show_plot=False,
    )

    # STEP 5 — DTEC
    print()
    print("----------------------------------------")
    print("pyOASIS STEP 5: DTECcalc")
    print("----------------------------------------")

    pyOASIS.DTECcalc(
        station,
        doy_str,
        year_str,
        str(oasis_rinex_output),
        str(dtec_output),
        show_plot=False,
    )

    # STEP 6 — SIDX
    print()
    print("----------------------------------------")
    print("pyOASIS STEP 6: SIDXcalc")
    print("----------------------------------------")

    pyOASIS.SIDXcalc(
        station,
        doy_str,
        year_str,
        str(oasis_rinex_output),
        str(sidx_output),
        show_plot=False,
    )

    # STEP 7 — TEC
    print()
    print("----------------------------------------")
    print("pyOASIS STEP 7: TECcalc")
    print("----------------------------------------")

    pyOASIS.TECcalc(
        station,
        doy_str,
        year_str,
        str(oasis_rinex_output),
        str(tec_output),
        show_plot=False,
    )

    final_tec = expected_tec_path(
        day_path,
        station,
        date,
    )

    print()
    print("========================================")
    print("FULL pyOASIS DAY FINISHED")
    print("========================================")
    print("DATE :", date.strftime("%Y-%m-%d"))
    print("RNX3 :", len(rnx3_files))
    print("ROTI :", roti_output)
    print("DTEC :", dtec_output)
    print("SIDX :", sidx_output)
    print("TEC  :", tec_output)

    if final_tec.is_file():
        print("FINAL TEC PRODUCT:")
        print(" ", final_tec)
    else:
        print(
            "WARNING: expected final TEC filename was not found:"
        )
        print(" ", final_tec)

    print("========================================")


def process_day(
    date: datetime,
    day_path: Path,
    config: Config,
) -> str:
    """Process one day and return a status string."""
    doy = int(date.strftime("%j"))

    print()
    print("========================================")
    print("DATE :", date.strftime("%Y-%m-%d"))
    print("DOY  :", f"{doy:03d}")
    print("========================================")

    final_tec = expected_tec_path(
        day_path,
        config.station,
        date,
    )

    if final_tec.is_file() and not config.force:
        print(
            "pyOASIS already completed - SKIPPING"
        )
        print("TEC:", final_tec)
        return "already-complete"

    rinex_folder = day_path / "RINEX"

    observation_files = find_observation_files(
        rinex_folder
    )

    if not observation_files:
        print(
            "NO RINEX OBSERVATION FILES"
        )
        return "no-rinex"

    selected_obs_path = select_observation_file(
        observation_files
    )

    if selected_obs_path is None:
        print(
            "NO SUITABLE 15 s / 30 s RINEX AVAILABLE"
        )
        print("Available files:")

        for path in observation_files:
            print(" ", path.name)

        print(
            "=> preprocessing/subsampling required"
        )
        print(
            "=> SP3 download skipped for this day"
        )

        return "no-suitable-rinex"

    print("pyOASIS INPUT:")
    print(" ", selected_obs_path)

    products_folder = (
        day_path
        / "PRODUCTS"
    )

    sp3_folder = (
        products_folder
        / "SP3"
    )

    # Retained for compatibility with the existing archive layout.
    (products_folder / "CLK").mkdir(
        parents=True,
        exist_ok=True,
    )
    (products_folder / "NAV").mkdir(
        parents=True,
        exist_ok=True,
    )

    sp3_path = download_sp3(
        date,
        sp3_folder,
        config,
    )

    if sp3_path is None:
        print(
            "SP3 NOT AVAILABLE / DOWNLOAD FAILED"
        )
        return "no-sp3"

    run_pyoasis_day(
        date,
        day_path,
        selected_obs_path,
        sp3_path,
        config,
    )

    return "processed"


def main() -> int:
    args = parse_args()
    config = validate_inputs(args)
    selected_date = parse_selected_date(
        args.date,
        config.year,
    )

    print("=" * 80)
    print("KOH2 pyOASIS PROCESSING WORKFLOW")
    print("=" * 80)
    print("Year       :", config.year)
    print("Year root  :", config.year_root)
    print("Station    :", config.station)
    print("GFZRNX     :", config.gfzrnx)
    print("CDDIS base :", config.cddis_base)
    print("Keep .gz   :", config.keep_gz)
    print("Force      :", config.force)
    print(
        "Date only  :",
        selected_date.strftime("%Y-%m-%d")
        if selected_date is not None
        else "ALL DAYS",
    )

    counters = {
        "processed": 0,
        "already-complete": 0,
        "no-rinex": 0,
        "no-suitable-rinex": 0,
        "no-sp3": 0,
        "failed": 0,
    }

    day_count = 0

    matched_selected_date = False

    for date, day_path in iter_day_folders(
        config.year_root,
        config.year,
    ):
        if (
            selected_date is not None
            and date.date() != selected_date.date()
        ):
            continue

        matched_selected_date = True
        day_count += 1

        try:
            status = process_day(
                date,
                day_path,
                config,
            )

            counters[status] += 1

        except Exception as exc:
            counters["failed"] += 1

            print()
            print("PROCESSING FAILED")
            print("DATE:", date.strftime("%Y-%m-%d"))
            print("DAY :", day_path)
            print(
                f"{type(exc).__name__}: {exc}"
            )
            print()

    if (
        selected_date is not None
        and not matched_selected_date
    ):
        raise FileNotFoundError(
            "No matching YEAR/MONTH/DAY directory was found for "
            f"{selected_date.strftime('%Y-%m-%d')} under:\n"
            f"    {config.year_root}"
        )

    print()
    print("=" * 80)
    print(
        "DATE SUMMARY"
        if selected_date is not None
        else "YEAR SUMMARY"
    )
    print("=" * 80)
    print("Day folders examined  :", day_count)

    for key, value in counters.items():
        print(
            f"{key:20s}: {value}"
        )

    if counters["failed"]:
        print()
        print(
            "One or more days failed. Review the messages above."
        )
        return 1

    print()
    print(
        "Processing completed without unhandled day-level failures."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
