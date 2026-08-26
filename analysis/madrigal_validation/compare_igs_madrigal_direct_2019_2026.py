from __future__ import annotations

r"""
DIRECT IGS FINAL GIM <-> MADRIGAL GNSS VTEC COMPARISON, KOH2 REGION, 2019-2026
=============================================================================

Purpose
-------
Remove PyTECGg and pyOASIS completely from the reference-product comparison.

For each available Madrigal 5-minute / 1-degree regional grid bin:
    1. use the Madrigal bin epoch, latitude, and longitude;
    2. interpolate IGS Final GIM VTEC to that exact epoch/location;
    3. calculate residual:

            IGS VTEC - Madrigal VTEC

This directly tests whether the growing 2023-2025 difference inferred from
the two KOH2 processors is actually an inter-reference effect.

Data reuse
----------
Madrigal:
    <MADRIGAL_VALIDATION_ROOT>\
        _MADRIGAL_CACHE\VTEC_REGIONAL_LOCAL\YYYY\
        madrigal_regional_local_YYYY_DDD.parquet

IGS:
    Reuses the validated V4 IGS validator and its cached IONEX files:
    <IGS_IONEX_ROOT>

No PyTECGg or pyOASIS values are used in the statistics.

Fixed regional box
------------------
To avoid a day-dependent footprint, use a constant box centered on KOH2:

    latitude : KOH2 latitude +/- 5 degrees
    longitude: KOH2 longitude +/- 10 degrees

This is intentionally inside the regional Madrigal extracts produced during
the prior KOH2 comparison and represents the broad KOH2/IPP neighborhood.

Outputs
-------
<IGS_MADRIGAL_REFERENCE_ROOT>

    KOH2_2019_2026_IGS_vs_Madrigal_daily_statistics.csv
    KOH2_2019_2026_IGS_vs_Madrigal_yearly_point_weighted_summary.csv
    KOH2_2019_2026_IGS_vs_Madrigal_equal_day_yearly_summary.csv
    KOH2_2019_2026_IGS_vs_Madrigal_monthly_summary.csv
    KOH2_2019_2026_IGS_vs_Madrigal_availability.csv
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
IGS_CACHE_ROOT = None

MADRIGAL_CACHE_ROOT = (
    BASE_ROOT
    / "TEC_VALIDATION_MADRIGAL_2019_2026"
    / "_MADRIGAL_CACHE"
    / "VTEC_REGIONAL_LOCAL"
)

OUTPUT_ROOT = (
    BASE_ROOT
    / "TEC_REFERENCE_IGS_MADRIGAL_2019_2026"
)

# Fixed regional box relative to KOH2.
LAT_HALF_WIDTH_DEG = 5.0
LON_HALF_WIDTH_DEG = 10.0


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Directly compare IGS Final GIM with cached Madrigal MAPGPS VTEC "
            "in the fixed KOH2 regional box."
        )
    )
    parser.add_argument(
        "--data-root",
        required=True,
        type=Path,
        help="KOH2 project data root.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: "
            "<data-root>/TEC_REFERENCE_IGS_MADRIGAL_2019_2026"
        ),
    )
    parser.add_argument(
        "--madrigal-cache-root",
        type=Path,
        default=None,
        help=(
            "Root containing VTEC_REGIONAL_LOCAL/YEAR/"
            "madrigal_regional_local_YEAR_DDD.parquet."
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
    global MADRIGAL_CACHE_ROOT, OUTPUT_ROOT, IGS_CACHE_ROOT

    BASE_ROOT = args.data_root.resolve()

    OUTPUT_ROOT = (
        args.output_root.resolve()
        if args.output_root is not None
        else BASE_ROOT / "TEC_REFERENCE_IGS_MADRIGAL_2019_2026"
    )

    MADRIGAL_CACHE_ROOT = (
        args.madrigal_cache_root.resolve()
        if args.madrigal_cache_root is not None
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

KOH2_LAT = CORE.KOH2_LAT
KOH2_LON = CORE.KOH2_LON

configure_year = CORE.configure_year
download_ionex = CORE.download_ionex
read_ionex = CORE.read_ionex
interpolate_ionex = CORE.interpolate_ionex
calculate_statistics = CORE.calculate_statistics
GlobalAccumulator = CORE.GlobalAccumulator


LAT_MIN = KOH2_LAT - LAT_HALF_WIDTH_DEG
LAT_MAX = KOH2_LAT + LAT_HALF_WIDTH_DEG
LON_MIN = KOH2_LON - LON_HALF_WIDTH_DEG
LON_MAX = KOH2_LON + LON_HALF_WIDTH_DEG


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


def discover_madrigal_days(
    year: int,
):
    folder = (
        MADRIGAL_CACHE_ROOT
        / str(
            year
        )
    )

    result = {}

    if not folder.is_dir():
        return result

    for path in folder.glob(
        f"madrigal_regional_local_{year}_*.parquet"
    ):
        m = re.search(
            rf"madrigal_regional_local_{year}_(\d{{3}})\.parquet$",
            path.name,
        )

        if not m:
            continue

        doy = int(
            m.group(
                1
            )
        )

        result[
            doy
        ] = path

    return dict(
        sorted(
            result.items()
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


def load_madrigal_fixed_region(
    path: Path,
):
    df = pd.read_parquet(
        path
    )

    required = {
        "epoch",
        "gdlat",
        "glon",
        "tec",
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

    mask = (
        df[
            "epoch"
        ].notna()
        & np.isfinite(
            df[
                "gdlat"
            ]
        )
        & np.isfinite(
            df[
                "glon"
            ]
        )
        & np.isfinite(
            df[
                "tec"
            ]
        )
        & (
            df[
                "gdlat"
            ]
            >= LAT_MIN
        )
        & (
            df[
                "gdlat"
            ]
            <= LAT_MAX
        )
        & (
            df[
                "glon"
            ]
            >= LON_MIN
        )
        & (
            df[
                "glon"
            ]
            <= LON_MAX
        )
    )

    return (
        df[
            mask
        ]
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


def daily_row(
    year: int,
    doy: int,
    df: pd.DataFrame,
):
    stats = calculate_statistics(
        df[
            "igs_vtec"
        ],
        df[
            "tec"
        ],
    )

    date = date_from_year_doy(
        year,
        doy,
    )

    residual = (
        df[
            "igs_vtec"
        ]
        - df[
            "tec"
        ]
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
            "IGS_Final_GIM_vs_Madrigal",
        "residual_convention":
            "IGS_minus_Madrigal",
        **stats,
        "igs_median_tecu":
            (
                float(
                    np.nanmedian(
                        df[
                            "igs_vtec"
                        ]
                    )
                )
                if len(
                    df
                )
                else np.nan
            ),
        "madrigal_median_tecu":
            (
                float(
                    np.nanmedian(
                        df[
                            "tec"
                        ]
                    )
                )
                if len(
                    df
                )
                else np.nan
            ),
        "median_madrigal_dtec_tecu":
            (
                float(
                    np.nanmedian(
                        df[
                            "dtec"
                        ]
                    )
                )
                if (
                    len(
                        df
                    )
                    and np.any(
                        np.isfinite(
                            df[
                                "dtec"
                            ]
                        )
                    )
                )
                else np.nan
            ),
        "mean_igs_tecu":
            (
                float(
                    np.nanmean(
                        df[
                            "igs_vtec"
                        ]
                    )
                )
                if len(
                    df
                )
                else np.nan
            ),
        "mean_madrigal_tecu":
            (
                float(
                    np.nanmean(
                        df[
                            "tec"
                        ]
                    )
                )
                if len(
                    df
                )
                else np.nan
            ),
        "median_abs_reference_difference_tecu":
            (
                float(
                    np.nanmedian(
                        np.abs(
                            residual
                        )
                    )
                )
                if len(
                    df
                )
                else np.nan
            ),
    }


# =============================================================================
# SUMMARIES
# =============================================================================

def point_weighted_summary(
    accumulators,
):
    rows = []

    for year, acc in sorted(
        accumulators.items()
    ):
        rows.append({
            "year":
                year,
            "comparison":
                "IGS_Final_GIM_vs_Madrigal",
            "residual_convention":
                "IGS_minus_Madrigal",
            **acc.statistics(),
        })

    return pd.DataFrame(
        rows
    )


def equal_day_summary(
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
            "comparison":
                "IGS_Final_GIM_vs_Madrigal",
            "residual_convention":
                "IGS_minus_Madrigal",
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
            "median_of_daily_igs_median_tecu":
                g[
                    "igs_median_tecu"
                ].median(),
            "median_of_daily_madrigal_median_tecu":
                g[
                    "madrigal_median_tecu"
                ].median(),
            "median_daily_madrigal_dtec_tecu":
                g[
                    "median_madrigal_dtec_tecu"
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
            "median_of_daily_igs_median_tecu":
                g[
                    "igs_median_tecu"
                ].median(),
            "median_of_daily_madrigal_median_tecu":
                g[
                    "madrigal_median_tecu"
                ].median(),
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

    # The imported publication IGS core owns download_ionex()/configure_year().
    # Point its cache globals to the selected validated cache before use.
    CORE.BASE_ROOT = BASE_ROOT
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
        "DIRECT IGS FINAL GIM vs MADRIGAL GNSS VTEC"
    )

    print(
        "=" * 118
    )

    print(
        "Residual convention: IGS - Madrigal"
    )

    print(
        f"Fixed region: lat {LAT_MIN:.4f} .. {LAT_MAX:.4f}, "
        f"lon {LON_MIN:.4f} .. {LON_MAX:.4f}"
    )

    print(
        "Madrigal regional cache:",
        MADRIGAL_CACHE_ROOT,
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

    daily_rows = []
    availability_rows = []
    accumulators = {}

    for year in YEARS:
        files = discover_madrigal_days(
            year
        )

        if SELECTED_DOY is not None:
            files = {
                doy: path
                for doy, path in files.items()
                if doy == SELECTED_DOY
            }

        print()
        print(
            "=" * 118
        )

        print(
            f"{year}: Madrigal cached days = {len(files)}"
        )

        print(
            "=" * 118
        )

        configure_year(
            year
        )

        accumulators[
            year
        ] = GlobalAccumulator()

        for i, (
            doy,
            path,
        ) in enumerate(
            files.items(),
            1,
        ):
            print(
                f"{i:3d}/{len(files):3d}  "
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
                "madrigal_cache":
                    str(
                        path
                    ),
                "status":
                    "FAILED",
                "n_madrigal_fixed_region":
                    0,
                "n_igs_matched":
                    0,
            }

            try:
                madrigal = load_madrigal_fixed_region(
                    path
                )

                status[
                    "n_madrigal_fixed_region"
                ] = len(
                    madrigal
                )

                if madrigal.empty:
                    raise RuntimeError(
                        "No Madrigal bins inside the fixed regional box."
                    )

                ionex_path = download_ionex(
                    year,
                    doy,
                )

                grid = read_ionex(
                    ionex_path
                )

                madrigal[
                    "igs_vtec"
                ] = interpolate_ionex(
                    grid,
                    madrigal[
                        "epoch"
                    ],
                    madrigal[
                        "gdlat"
                    ],
                    madrigal[
                        "glon"
                    ],
                )

                matched = madrigal.dropna(
                    subset=[
                        "igs_vtec",
                        "tec",
                    ]
                ).copy()

                status[
                    "n_igs_matched"
                ] = len(
                    matched
                )

                if matched.empty:
                    raise RuntimeError(
                        "No valid IGS-Madrigal matched bins."
                    )

                row = daily_row(
                    year,
                    doy,
                    matched,
                )

                daily_rows.append(
                    row
                )

                accumulators[
                    year
                ].add(
                    matched[
                        "igs_vtec"
                    ],
                    matched[
                        "tec"
                    ],
                )

                status[
                    "status"
                ] = "OK"

                print(
                    f"    n={row['n']:,}  "
                    f"bias={row['bias_tecu']:+.3f}  "
                    f"RMSE={row['rmse_tecu']:.3f}  "
                    f"r={row['pearson_r']:.3f}"
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

        # Restart-friendly snapshots after each year.
        pd.DataFrame(
            daily_rows
        ).to_csv(
            OUTPUT_ROOT
            / (
                f"{STATION}_2019_2026_"
                "IGS_vs_Madrigal_daily_statistics.csv"
            ),
            index=False,
        )

        pd.DataFrame(
            availability_rows
        ).to_csv(
            OUTPUT_ROOT
            / (
                f"{STATION}_2019_2026_"
                "IGS_vs_Madrigal_availability.csv"
            ),
            index=False,
        )

    daily = pd.DataFrame(
        daily_rows
    )

    availability = pd.DataFrame(
        availability_rows
    )

    point = point_weighted_summary(
        accumulators
    )

    equal = equal_day_summary(
        daily
    )

    monthly = monthly_summary(
        daily
    )

    daily_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "IGS_vs_Madrigal_daily_statistics.csv"
        )
    )

    availability_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "IGS_vs_Madrigal_availability.csv"
        )
    )

    point_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "IGS_vs_Madrigal_yearly_point_weighted_summary.csv"
        )
    )

    equal_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "IGS_vs_Madrigal_equal_day_yearly_summary.csv"
        )
    )

    monthly_file = (
        OUTPUT_ROOT
        / (
            f"{STATION}_2019_2026_"
            "IGS_vs_Madrigal_monthly_summary.csv"
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

    point.to_csv(
        point_file,
        index=False,
    )

    equal.to_csv(
        equal_file,
        index=False,
    )

    monthly.to_csv(
        monthly_file,
        index=False,
    )

    print()
    print(
        "=" * 118
    )

    print(
        "DIRECT IGS-MADRIGAL COMPARISON COMPLETE"
    )

    print(
        "=" * 118
    )

    print(
        "Daily:",
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

    print(
        "Monthly:",
        monthly_file,
    )

    if not equal.empty:
        print()
        print(
            "EQUAL-DAY YEARLY SUMMARY"
        )

        print(
            "-" * 118
        )

        print(
            equal.to_string(
                index=False
            )
        )

        print(
            "-" * 118
        )


if __name__ == "__main__":
    main()
