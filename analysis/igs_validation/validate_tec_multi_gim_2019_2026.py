from __future__ import annotations

r"""
KOH2 2019-2026 MULTI-GIM COMPARISON
===================================

References:
    CODE final GIM  (CODG / COD0OPSFIN)
    ESA final GIM   (ESAG / ESA0OPSFIN)
    JPL final GIM   (JPLG / JPL0OPSFIN)
    UPC final GIM   (UPCG / UPC0OPSFIN)

This script is intended to be run after:
    validate_tec_igs_2019_2026_E_drive_V4.py

It deliberately reuses the exact same:
    - KOH2 coordinates
    - PyTECGg reader
    - pyOASIS reader
    - IONEX parser
    - spatial/temporal interpolation
    - elevation cutoff
    - residual/statistics definitions

from that validated IGS script, so the comparison methodology remains
consistent across all GIM references.

IMPORTANT SCIENTIFIC NOTE
-------------------------
These individual GIMs are products of IGS ionospheric analysis centers.
The IGS combined final GIM is not statistically independent of the IAAC
products. Therefore call this a "multi-GIM" or "inter-product comparison",
not four independent validations.

Residual convention:
    method VTEC - reference GIM VTEC

Root production data:
    E:\KOH2data\2019
    ...
    E:\KOH2data\2026

Outputs:
    E:\KOH2data\TEC_VALIDATION_MULTI_GIM_2019_2026
"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
from urllib.parse import urlparse
import importlib.util
import math
import netrc
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

MASTER_OUTPUT = (
    BASE_ROOT
    / "TEC_VALIDATION_MULTI_GIM_2019_2026"
)

CACHE_ROOT = (
    MASTER_OUTPUT
    / "_IONEX"
)

BY_YEAR_OUTPUT = (
    MASTER_OUTPUT
    / "BY_YEAR"
)

CDDIS_BASE = (
    "https://cddis.nasa.gov/archive/gnss/products/ionex"
)

MIN_ELEVATION_DEG = 30.0

SAVE_DAILY_MATCHED = False

DOWNLOAD_ATTEMPTS_PER_CANDIDATE = 3
DOWNLOAD_CONNECT_TIMEOUT_S = 30
DOWNLOAD_READ_TIMEOUT_S = 60
DOWNLOAD_RETRY_WAIT_S = 5
DOWNLOAD_CHUNK_BYTES = 256 * 1024


# Candidate temporal intervals are deliberately permissive.  Historical
# products use the short name, which does not encode sampling.  Long-format
# products do encode sampling and some centers/products have changed submission
# conventions.  Trying both 01H and 02H avoids hard-coding a false assumption.
GIMS = {
    "CODE": {
        "long_ac": "COD",
        "legacy": "codg",
        "sampling_candidates": ["01H", "02H"],
    },
    "ESA": {
        "long_ac": "ESA",
        "legacy": "esag",
        "sampling_candidates": ["02H", "01H"],
    },
    "JPL": {
        "long_ac": "JPL",
        "legacy": "jplg",
        "sampling_candidates": ["02H", "01H"],
    },
    "UPC": {
        "long_ac": "UPC",
        "legacy": "upcg",
        "sampling_candidates": ["02H", "01H"],
    },
}


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare KOH2 PyTECGg and pyOASIS VTEC with final "
            "CODE/ESA/JPL/UPC GIM products."
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
            "<data-root>/TEC_VALIDATION_MULTI_GIM_2019_2026"
        ),
    )
    parser.add_argument(
        "--gim-cache-root",
        type=Path,
        default=None,
        help=(
            "Optional existing multi-GIM IONEX cache root. "
            "Expected layout: <root>/<GIM>/<YEAR>/..."
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
    global BASE_ROOT, YEARS, MASTER_OUTPUT, CACHE_ROOT, BY_YEAR_OUTPUT
    global SELECTED_DOY

    BASE_ROOT = args.data_root.resolve()

    MASTER_OUTPUT = (
        args.output_root.resolve()
        if args.output_root is not None
        else BASE_ROOT / "TEC_VALIDATION_MULTI_GIM_2019_2026"
    )

    CACHE_ROOT = (
        args.gim_cache_root.resolve()
        if args.gim_cache_root is not None
        else MASTER_OUTPUT / "_IONEX"
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


# =============================================================================
# LOAD THE ALREADY-VALIDATED PUBLICATION IGS CORE
# =============================================================================

def load_validation_core():
    """
    Import the validated publication IGS core without executing its main().
    """
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
            "\nCould not find the validated publication IGS core script:\n"
            f"    {filename}\n\n"
            "Place this multi-GIM script in the same folder as V4, or run it "
            "from the folder containing V4."
        )

    spec = importlib.util.spec_from_file_location(
        "koh2_igs_v4_core",
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

# Reuse exact validated constants/functions.
KOH2_LAT = CORE.KOH2_LAT
KOH2_LON = CORE.KOH2_LON
KOH2_H = CORE.KOH2_H

load_pytecgg = CORE.load_pytecgg
load_pyoasis = CORE.load_pyoasis
read_ionex = CORE.read_ionex
interpolate_ionex = CORE.interpolate_ionex
calculate_statistics = CORE.calculate_statistics
GlobalAccumulator = CORE.GlobalAccumulator
decompress_ionex = CORE.decompress_ionex
get_cddis_session = CORE.get_cddis_session


# =============================================================================
# FILE DISCOVERY
# =============================================================================

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
        m.group(1)
    )

    y = int(
        m.group(2)
    )

    if y != year:
        return None

    return doy


def discover_data_files(
    year: int,
):
    root = (
        BASE_ROOT
        / str(year)
    )

    pytecgg = {}
    pyoasis = {}

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


# =============================================================================
# GIM DOWNLOAD
# =============================================================================

def gim_candidate_names(
    model: str,
    year: int,
    doy: int,
):
    cfg = GIMS[
        model
    ]

    ac = cfg[
        "long_ac"
    ]

    legacy = cfg[
        "legacy"
    ]

    yy = (
        year
        % 100
    )

    names = []

    for sampling in cfg[
        "sampling_candidates"
    ]:
        names.extend([
            (
                f"{ac}0OPSFIN_"
                f"{year}{doy:03d}0000_"
                f"01D_{sampling}_GIM.INX.gz"
            ),
            (
                f"{ac}0OPSFIN_"
                f"{year}{doy:03d}0000_"
                f"01D_{sampling}_ION.IOX.gz"
            ),
        ])

    # Historical IGS short product filename.
    names.extend([
        f"{legacy}{doy:03d}0.{yy:02d}i.Z",
        f"{legacy.upper()}{doy:03d}0.{yy:02d}I.Z",
    ])

    # Preserve order but remove duplicate candidates.
    return list(
        dict.fromkeys(
            names
        )
    )


def gim_url(
    year: int,
    doy: int,
    name: str,
):
    return (
        f"{CDDIS_BASE}/"
        f"{year}/"
        f"{doy:03d}/"
        f"{name}"
    )


def decompressed_path_for(
    path: Path,
):
    low = (
        path.name.lower()
    )

    if low.endswith(
        ".gz"
    ):
        return path.with_suffix(
            ""
        )

    if path.name.endswith(
        ".Z"
    ) or path.name.endswith(
        ".z"
    ):
        return path.with_suffix(
            ""
        )

    return path


def try_existing_gim(
    model: str,
    year: int,
    doy: int,
    names,
):
    cache = (
        CACHE_ROOT
        / model
        / str(year)
    )

    for name in names:
        compressed = (
            cache
            / name
        )

        output = decompressed_path_for(
            compressed
        )

        if (
            output.is_file()
            and output.stat().st_size > 1000
        ):
            return output

        if (
            compressed.is_file()
            and compressed.stat().st_size > 1000
        ):
            try:
                return decompress_ionex(
                    compressed,
                    output,
                )
            except Exception:
                pass

    return None


def download_gim(
    model: str,
    year: int,
    doy: int,
):
    cache = (
        CACHE_ROOT
        / model
        / str(year)
    )

    cache.mkdir(
        parents=True,
        exist_ok=True,
    )

    names = gim_candidate_names(
        model,
        year,
        doy,
    )

    existing = try_existing_gim(
        model,
        year,
        doy,
        names,
    )

    if existing is not None:
        return existing

    session = get_cddis_session()

    errors = []

    for name in names:
        url = gim_url(
            year,
            doy,
            name,
        )

        compressed = (
            cache
            / name
        )

        output = decompressed_path_for(
            compressed
        )

        for attempt in range(
            1,
            DOWNLOAD_ATTEMPTS_PER_CANDIDATE
            + 1,
        ):
            tmp = compressed.with_suffix(
                compressed.suffix
                + ".part"
            )

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
                    break

                response.raise_for_status()

                content_type = response.headers.get(
                    "Content-Type",
                    "",
                ).lower()

                first_chunk = True
                bytes_written = 0

                with open(
                    tmp,
                    "wb",
                ) as f:
                    for chunk in response.iter_content(
                        chunk_size=DOWNLOAD_CHUNK_BYTES
                    ):
                        if not chunk:
                            continue

                        if first_chunk:
                            head = chunk[
                                :500
                            ].lower()

                            if (
                                b"<html" in head
                                or b"earthdata login" in head
                                or "text/html" in content_type
                            ):
                                raise RuntimeError(
                                    "CDDIS returned HTML instead of IONEX."
                                )

                            first_chunk = False

                        f.write(
                            chunk
                        )

                        bytes_written += len(
                            chunk
                        )

                if bytes_written < 1000:
                    raise RuntimeError(
                        f"Downloaded file too small: "
                        f"{bytes_written} bytes"
                    )

                tmp.replace(
                    compressed
                )

                print(
                    f"      {model}: downloaded "
                    f"{compressed.name}"
                )

                return decompress_ionex(
                    compressed,
                    output,
                )

            except KeyboardInterrupt:
                if tmp.exists():
                    try:
                        tmp.unlink()
                    except Exception:
                        pass
                raise

            except Exception as exc:
                errors.append(
                    (
                        f"{name} attempt {attempt}/"
                        f"{DOWNLOAD_ATTEMPTS_PER_CANDIDATE}: "
                        f"{exc!r}"
                    )
                )

                if tmp.exists():
                    try:
                        tmp.unlink()
                    except Exception:
                        pass

                if (
                    attempt
                    < DOWNLOAD_ATTEMPTS_PER_CANDIDATE
                ):
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
        (
            f"No {model} final GIM could be obtained "
            f"for {year} DOY {doy:03d}."
        )
    )


# =============================================================================
# DAILY STATISTICS
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


def daily_statistics_row(
    year,
    doy,
    model,
    comparison,
    method,
    reference,
):
    stats = calculate_statistics(
        method,
        reference,
    )

    method_clean = np.asarray(
        method,
        dtype=float,
    )

    ref_clean = np.asarray(
        reference,
        dtype=float,
    )

    good = (
        np.isfinite(
            method_clean
        )
        & np.isfinite(
            ref_clean
        )
    )

    method_clean = method_clean[
        good
    ]

    ref_clean = ref_clean[
        good
    ]

    date = date_from_year_doy(
        year,
        doy,
    )

    return {
        "year":
            year,
        "month":
            date.month,
        "date":
            date.date().isoformat(),
        "doy":
            doy,
        "reference_gim":
            model,
        "comparison":
            comparison,
        **stats,
        "method_median_tecu":
            (
                float(
                    np.median(
                        method_clean
                    )
                )
                if len(
                    method_clean
                )
                else np.nan
            ),
        "reference_median_tecu":
            (
                float(
                    np.median(
                        ref_clean
                    )
                )
                if len(
                    ref_clean
                )
                else np.nan
            ),
    }


# =============================================================================
# METHOD-SPECIFIC VALIDATION
# =============================================================================

def validate_pytecgg(
    pt: pd.DataFrame,
    grid,
    year: int,
    doy: int,
    model: str,
    accumulators,
    daily_rows,
    daily_output: Path | None,
):
    # A) PyTECGg VEq at KOH2 station position.
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
        "reference_vtec"
    ] = interpolate_ionex(
        grid,
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

    station = station.dropna(
        subset=[
            "veq",
            "reference_vtec",
        ]
    )

    comparison = (
        "PyTECGg_VEq_vs_GIM_station"
    )

    accumulators[
        comparison
    ].add(
        station[
            "veq"
        ],
        station[
            "reference_vtec"
        ],
    )

    daily_rows.append(
        daily_statistics_row(
            year,
            doy,
            model,
            comparison,
            station[
                "veq"
            ],
            station[
                "reference_vtec"
            ],
        )
    )

    # B) PyTECGg VTEC at PyTECGg IPP.
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
    ].copy()

    ipp[
        "reference_vtec"
    ] = interpolate_ionex(
        grid,
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

    ipp = ipp.dropna(
        subset=[
            "vtec",
            "reference_vtec",
        ]
    )

    comparison = (
        "PyTECGg_VTEC_vs_GIM_IPP"
    )

    accumulators[
        comparison
    ].add(
        ipp[
            "vtec"
        ],
        ipp[
            "reference_vtec"
        ],
    )

    daily_rows.append(
        daily_statistics_row(
            year,
            doy,
            model,
            comparison,
            ipp[
                "vtec"
            ],
            ipp[
                "reference_vtec"
            ],
        )
    )

    if (
        SAVE_DAILY_MATCHED
        and daily_output is not None
    ):
        daily_output.mkdir(
            parents=True,
            exist_ok=True,
        )

        ipp.to_parquet(
            daily_output
            / (
                f"{STATION}_{year}{doy:03d}_"
                f"PyTECGg_vs_{model}.parquet"
            ),
            index=False,
        )

    return (
        len(
            station
        ),
        len(
            ipp
        ),
    )


def validate_pyoasis(
    po: pd.DataFrame,
    grid,
    year: int,
    doy: int,
    model: str,
    accumulators,
    daily_rows,
    daily_output: Path | None,
):
    po = po[
        [
            "epoch",
            "sat",
            "lat_ipp",
            "lon_ipp",
            "elevation",
            "vtec",
        ]
    ].copy()

    po = po[
        pd.to_numeric(
            po[
                "elevation"
            ],
            errors="coerce",
        )
        >= MIN_ELEVATION_DEG
    ].copy()

    po[
        "reference_vtec"
    ] = interpolate_ionex(
        grid,
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

    po = po.dropna(
        subset=[
            "vtec",
            "reference_vtec",
        ]
    )

    comparison = (
        "pyOASIS_VTEC_vs_GIM_IPP"
    )

    accumulators[
        comparison
    ].add(
        po[
            "vtec"
        ],
        po[
            "reference_vtec"
        ],
    )

    daily_rows.append(
        daily_statistics_row(
            year,
            doy,
            model,
            comparison,
            po[
                "vtec"
            ],
            po[
                "reference_vtec"
            ],
        )
    )

    if (
        SAVE_DAILY_MATCHED
        and daily_output is not None
    ):
        daily_output.mkdir(
            parents=True,
            exist_ok=True,
        )

        po.to_parquet(
            daily_output
            / (
                f"{STATION}_{year}{doy:03d}_"
                f"pyOASIS_vs_{model}.parquet"
            ),
            index=False,
        )

    return len(
        po
    )


# =============================================================================
# SUMMARY FUNCTIONS
# =============================================================================

def point_weighted_summary(
    all_accumulators,
):
    rows = []

    for (
        year,
        model,
        comparison,
    ), accumulator in all_accumulators.items():
        rows.append({
            "year":
                year,
            "reference_gim":
                model,
            "comparison":
                comparison,
            **accumulator.statistics(),
        })

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "reference_gim",
            "comparison",
        ]
    )


def equal_day_summary(
    daily: pd.DataFrame,
):
    if daily.empty:
        return pd.DataFrame()

    rows = []

    for (
        year,
        model,
        comparison,
    ), g in daily.groupby(
        [
            "year",
            "reference_gim",
            "comparison",
        ],
        dropna=False,
    ):
        rows.append({
            "year":
                int(
                    year
                ),
            "reference_gim":
                model,
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
            "median_of_daily_reference_median_tecu":
                g[
                    "reference_median_tecu"
                ].median(),
        })

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "comparison",
            "reference_gim",
        ]
    )


def strict_common_day_summary(
    daily: pd.DataFrame,
):
    """
    For each year + comparison, retain only dates having valid statistics from
    ALL four GIM references.  This is the fairest inter-GIM comparison.
    """
    if daily.empty:
        return pd.DataFrame()

    rows = []

    expected_models = set(
        GIMS
    )

    for (
        year,
        comparison,
    ), g in daily.groupby(
        [
            "year",
            "comparison",
        ]
    ):
        valid = g[
            pd.to_numeric(
                g[
                    "n"
                ],
                errors="coerce",
            )
            > 0
        ].copy()

        dates = []

        for date, gd in valid.groupby(
            "date"
        ):
            models = set(
                gd[
                    "reference_gim"
                ].dropna()
            )

            if expected_models.issubset(
                models
            ):
                dates.append(
                    date
                )

        if not dates:
            continue

        common = valid[
            valid[
                "date"
            ].isin(
                dates
            )
        ].copy()

        for model, gm in common.groupby(
            "reference_gim"
        ):
            rows.append({
                "year":
                    int(
                        year
                    ),
                "comparison":
                    comparison,
                "reference_gim":
                    model,
                "n_common_days":
                    int(
                        len(
                            dates
                        )
                    ),
                "n_points_total":
                    int(
                        pd.to_numeric(
                            gm[
                                "n"
                            ],
                            errors="coerce",
                        ).fillna(
                            0
                        ).sum()
                    ),
                "mean_daily_bias_tecu":
                    gm[
                        "bias_tecu"
                    ].mean(),
                "median_daily_bias_tecu":
                    gm[
                        "bias_tecu"
                    ].median(),
                "mean_daily_rmse_tecu":
                    gm[
                        "rmse_tecu"
                    ].mean(),
                "median_daily_rmse_tecu":
                    gm[
                        "rmse_tecu"
                    ].median(),
                "median_daily_pearson_r":
                    gm[
                        "pearson_r"
                    ].median(),
                "median_of_daily_reference_median_tecu":
                    gm[
                        "reference_median_tecu"
                    ].median(),
            })

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "comparison",
            "reference_gim",
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
        model,
        comparison,
    ), g in daily.groupby(
        [
            "year",
            "month",
            "reference_gim",
            "comparison",
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
            "reference_gim":
                model,
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
            "median_of_daily_reference_median_tecu":
                g[
                    "reference_median_tecu"
                ].median(),
        })

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "year",
            "month",
            "comparison",
            "reference_gim",
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

    BY_YEAR_OUTPUT.mkdir(
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
        "KOH2 2019-2026 MULTI-GIM COMPARISON"
    )

    print(
        "=" * 118
    )

    print(
        "Root:",
        BASE_ROOT,
    )

    print(
        "Output:",
        MASTER_OUTPUT,
    )

    print(
        "GIM cache:",
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
        "References:",
        ", ".join(
            GIMS
        ),
    )

    print(
        "Elevation cutoff:",
        MIN_ELEVATION_DEG,
        "deg",
    )

    print(
        "Residual convention: method VTEC - reference GIM VTEC"
    )

    print(
        "=" * 118
    )

    all_accumulators = {}
    daily_rows = []
    availability_rows = []

    for year in YEARS:
        root = (
            BASE_ROOT
            / str(year)
        )

        print()
        print(
            "=" * 118
        )

        print(
            year
        )

        print(
            "=" * 118
        )

        if not root.is_dir():
            print(
                "Year root missing:",
                root,
            )
            continue

        (
            pytecgg_files,
            pyoasis_files,
            days,
        ) = discover_data_files(
            year
        )

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

        print(
            "PyTECGg days:",
            len(
                pytecgg_files
            ),
        )

        print(
            "pyOASIS days:",
            len(
                pyoasis_files
            ),
        )

        print(
            "Union days:",
            len(
                days
            ),
        )

        # Cache loaded method data once per day; do not reread for every GIM.
        for index, doy in enumerate(
            days,
            1,
        ):
            print()
            print(
                f"{index:3d}/{len(days):3d}  "
                f"{year} DOY {doy:03d}"
            )

            pt = None
            po = None

            if doy in pytecgg_files:
                try:
                    pt = load_pytecgg(
                        pytecgg_files[
                            doy
                        ]
                    )
                except Exception as exc:
                    print(
                        "    [ERROR] PyTECGg read:",
                        repr(
                            exc
                        ),
                    )

            if doy in pyoasis_files:
                try:
                    po = load_pyoasis(
                        pyoasis_files[
                            doy
                        ]
                    )
                except Exception as exc:
                    print(
                        "    [ERROR] pyOASIS read:",
                        repr(
                            exc
                        ),
                    )

            for model in GIMS:
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
                    "reference_gim":
                        model,
                    "pytecgg_present":
                        pt is not None,
                    "pyoasis_present":
                        po is not None,
                    "gim_status":
                        "FAILED",
                    "pytecgg_veq_matches":
                        0,
                    "pytecgg_vtec_matches":
                        0,
                    "pyoasis_vtec_matches":
                        0,
                }

                try:
                    gim_path = download_gim(
                        model,
                        year,
                        doy,
                    )

                    grid = read_ionex(
                        gim_path
                    )

                    status[
                        "gim_status"
                    ] = "OK"

                    status[
                        "gim_file"
                    ] = str(
                        gim_path
                    )

                    for comparison in [
                        "PyTECGg_VEq_vs_GIM_station",
                        "PyTECGg_VTEC_vs_GIM_IPP",
                        "pyOASIS_VTEC_vs_GIM_IPP",
                    ]:
                        key = (
                            year,
                            model,
                            comparison,
                        )

                        if key not in all_accumulators:
                            all_accumulators[
                                key
                            ] = GlobalAccumulator()

                    daily_output = (
                        BY_YEAR_OUTPUT
                        / str(year)
                        / "DAILY_MATCHED"
                        / model
                    )

                    if pt is not None:
                        (
                            n_veq,
                            n_pt,
                        ) = validate_pytecgg(
                            pt,
                            grid,
                            year,
                            doy,
                            model,
                            {
                                k[
                                    2
                                ]:
                                    v
                                for k, v in all_accumulators.items()
                                if (
                                    k[
                                        0
                                    ]
                                    == year
                                    and k[
                                        1
                                    ]
                                    == model
                                )
                            },
                            daily_rows,
                            daily_output,
                        )

                        status[
                            "pytecgg_veq_matches"
                        ] = n_veq

                        status[
                            "pytecgg_vtec_matches"
                        ] = n_pt

                    if po is not None:
                        n_po = validate_pyoasis(
                            po,
                            grid,
                            year,
                            doy,
                            model,
                            {
                                k[
                                    2
                                ]:
                                    v
                                for k, v in all_accumulators.items()
                                if (
                                    k[
                                        0
                                    ]
                                    == year
                                    and k[
                                        1
                                    ]
                                    == model
                                )
                            },
                            daily_rows,
                            daily_output,
                        )

                        status[
                            "pyoasis_vtec_matches"
                        ] = n_po

                    print(
                        f"    {model:4s}: "
                        f"PyVEq={status['pytecgg_veq_matches']:,}  "
                        f"PyVTEC={status['pytecgg_vtec_matches']:,}  "
                        f"pyOASIS={status['pyoasis_vtec_matches']:,}"
                    )

                except Exception as exc:
                    status[
                        "error"
                    ] = repr(
                        exc
                    )

                    print(
                        f"    {model:4s}: [UNAVAILABLE/ERROR] "
                        f"{exc}"
                    )

                availability_rows.append(
                    status
                )

        # Save a restart-friendly snapshot after each year.
        daily_snapshot = pd.DataFrame(
            daily_rows
        )

        availability_snapshot = pd.DataFrame(
            availability_rows
        )

        daily_snapshot.to_csv(
            MASTER_OUTPUT
            / (
                f"{STATION}_2019_2026_"
                "multiGIM_daily_statistics.csv"
            ),
            index=False,
        )

        availability_snapshot.to_csv(
            MASTER_OUTPUT
            / (
                f"{STATION}_2019_2026_"
                "multiGIM_availability.csv"
            ),
            index=False,
        )

    # -------------------------------------------------------------------------
    # FINAL TABLES
    # -------------------------------------------------------------------------

    point = point_weighted_summary(
        all_accumulators
    )

    daily = pd.DataFrame(
        daily_rows
    )

    availability = pd.DataFrame(
        availability_rows
    )

    equal = equal_day_summary(
        daily
    )

    strict = strict_common_day_summary(
        daily
    )

    monthly = monthly_summary(
        daily
    )

    point_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "multiGIM_yearly_point_weighted_summary.csv"
        )
    )

    equal_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "multiGIM_equal_day_yearly_summary.csv"
        )
    )

    strict_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "multiGIM_strict_common_day_yearly_summary.csv"
        )
    )

    monthly_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "multiGIM_monthly_summary.csv"
        )
    )

    daily_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "multiGIM_daily_statistics.csv"
        )
    )

    availability_file = (
        MASTER_OUTPUT
        / (
            f"{STATION}_2019_2026_"
            "multiGIM_availability.csv"
        )
    )

    point.to_csv(
        point_file,
        index=False,
    )

    equal.to_csv(
        equal_file,
        index=False,
    )

    strict.to_csv(
        strict_file,
        index=False,
    )

    monthly.to_csv(
        monthly_file,
        index=False,
    )

    daily.to_csv(
        daily_file,
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
        "MULTI-GIM COMPARISON COMPLETE"
    )

    print(
        "=" * 118
    )

    print(
        "Point-weighted:",
        point_file,
    )

    print(
        "Equal-day:",
        equal_file,
    )

    print(
        "STRICT common-day:",
        strict_file,
    )

    print(
        "Monthly:",
        monthly_file,
    )

    print(
        "Daily:",
        daily_file,
    )

    print(
        "Availability:",
        availability_file,
    )

    print()

    if not strict.empty:
        print(
            "STRICT COMMON-DAY YEARLY SUMMARY"
        )

        print(
            "-" * 118
        )

        print(
            strict.to_string(
                index=False
            )
        )

        print(
            "-" * 118
        )


if __name__ == "__main__":
    main()
