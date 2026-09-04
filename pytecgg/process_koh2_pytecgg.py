# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Asparuh Kamburov
#
# Original KOH2-ionosphere workflow/orchestration code.
# PyTECGg is a separately installed external dependency and is not vendored
# in this repository. Upstream PyTECGg is licensed under GPL-3.0-or-later:
# https://github.com/viventriglia/PyTECGg

"""
PyTECGg batch processor for YEAR/MM/DD/RINEX archives.

Per day it processes ALL available 30-second RINEX observation segments,
concatenates them, removes duplicate observations, downloads BRDC navigation
from BKG, computes geometry, arcs, calibrated sTEC/vTEC and VEq, and writes
CSV/Parquet outputs.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime
import argparse
import traceback

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd

import polars as pl

from pytecgg import GNSSContext
from pytecgg.parsing import read_rinex_obs, read_rinex_nav
from pytecgg.utils import download_nav_bkg
from pytecgg.linear_combinations import calculate_linear_combinations
from pytecgg.satellites import prepare_ephemeris, satellite_coordinates, calculate_ipp
from pytecgg.tec_calibration import extract_arcs, calculate_tec, calculate_vertical_equivalent


# ============================================================
# SCIENTIFIC PROCESSING SETTINGS
# ============================================================
#
# These defaults reproduce the operational KOH2 workflow.
# Paths, year and optional single-date selection are supplied by CLI.

BASE: Path | None = None
YEAR: int | None = None

SAMPLING_TOKEN = "_30S_"
H_IPP = 350_000.0
MIN_ELEVATION = 30.0
SELECTION_MODE = "quality"
MIN_ARC_LENGTH = 30
THRESHOLD_ABS = 5.0
THRESHOLD_STD = 5.0
THRESHOLD_JUMP = 10.0
MAX_POLYNOMIAL_DEGREE = 3
BATCH_SIZE_EPOCHS = 30
SKIP_COMPLETED = True

SUPPORTED_OBS_ENDINGS = (".rnx", ".rnx.gz", ".crx", ".crx.gz")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Process KOH2 30-second RINEX observations with PyTECGg. "
            "The scientific defaults reproduce the operational workflow."
        )
    )
    parser.add_argument(
        "--data-root",
        required=True,
        type=Path,
        help="Parent of YEAR/MM/DD directories, or the YEAR directory itself.",
    )
    parser.add_argument(
        "--year",
        required=True,
        type=int,
        help="Four-digit processing year.",
    )
    parser.add_argument(
        "--date",
        type=str,
        default=None,
        help="Optional single date in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recalculate even when the expected TEC output already exists.",
    )
    return parser.parse_args()


def resolve_year_root(data_root: Path, year: int) -> Path:
    root = data_root.resolve()

    if root.name == str(year):
        return root

    return root / str(year)


def is_supported_obs(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in SUPPORTED_OBS_ENDINGS)


def find_30s_observation_files(rinex_dir: Path) -> list[Path]:
    if not rinex_dir.is_dir():
        return []
    return sorted(
        p for p in rinex_dir.iterdir()
        if p.is_file() and is_supported_obs(p) and SAMPLING_TOKEN in p.name.upper()
    )


def parse_all_observations(obs_files: list[Path]):
    frames = []
    receiver_pos = None
    rinex_version = None
    receiver_name = None

    for i, obs_path in enumerate(obs_files, start=1):
        print(f"  OBS {i:02d}/{len(obs_files):02d}: {obs_path.name}")
        df_obs, rec_pos, version = read_rinex_obs(obs_path)

        if df_obs.is_empty():
            print("      WARNING: empty observation table; skipping file")
            continue

        frames.append(df_obs)

        if receiver_pos is None:
            receiver_pos = rec_pos
            rinex_version = version
            receiver_name = obs_path.name[:4].lower()
        elif version != rinex_version:
            print(f"      WARNING: mixed RINEX versions: {rinex_version} and {version}")

    if not frames:
        raise RuntimeError("No usable observation records were parsed.")

    obs_all = pl.concat(frames, how="vertical_relaxed")
    rows_before = obs_all.height

    dedup_cols = [c for c in ("epoch", "sv", "observable") if c in obs_all.columns]
    if len(dedup_cols) == 3:
        obs_all = obs_all.unique(subset=dedup_cols, keep="first", maintain_order=True)

    if all(c in obs_all.columns for c in ("epoch", "sv", "observable")):
        obs_all = obs_all.sort(["epoch", "sv", "observable"])

    removed = rows_before - obs_all.height
    print(f"  Parsed observation records : {rows_before:,}")
    print(f"  Duplicate records removed  : {removed:,}")
    print(f"  Records entering PyTECGg   : {obs_all.height:,}")

    return obs_all, receiver_pos, rinex_version, receiver_name


def merge_navigation_files(nav_files: list[Path]) -> dict:
    full_nav: dict[str, pl.DataFrame] = {}
    parsed_any = False

    for nav_path in sorted(nav_files):
        if not nav_path.is_file() or nav_path.stat().st_size == 0:
            continue

        try:
            print("  NAV:", nav_path.name)
            day_nav = read_rinex_nav(nav_path)
        except Exception as exc:
            print(f"      Cannot parse as RINEX NAV: {exc}")
            continue

        if not day_nav:
            continue

        parsed_any = True
        for system, df_nav in day_nav.items():
            if df_nav is None or df_nav.is_empty():
                continue
            if system not in full_nav or full_nav[system].is_empty():
                full_nav[system] = df_nav
            else:
                full_nav[system] = pl.concat([full_nav[system], df_nav], how="vertical_relaxed")

    if not parsed_any or not full_nav:
        raise RuntimeError("No usable broadcast navigation data were parsed.")

    return full_nav


def systems_from_observations(df_obs: pl.DataFrame) -> list[str]:
    supported = {"G", "R", "E", "C"}
    systems = sorted({
        str(sv)[0].upper()
        for sv in df_obs["sv"].drop_nulls().to_list()
        if str(sv) and str(sv)[0].upper() in supported
    })
    if not systems:
        raise RuntimeError("No supported GNSS satellites found.")
    return systems


def save_manifest(path: Path, date: datetime, doy: int, obs_files: list[Path], systems: list[str], obs_rows: int):
    lines = [
        f"Date: {date:%Y-%m-%d}",
        f"DOY: {doy:03d}",
        f"Systems: {','.join(systems)}",
        f"Observation records after de-duplication: {obs_rows}",
        f"Elevation mask: {MIN_ELEVATION:.1f} deg",
        f"IPP height: {H_IPP:.1f} m",
        "",
        "30-second RINEX segments used:",
    ]
    lines.extend(f"  {p}" for p in obs_files)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ============================================================
# PLOTTING
# ============================================================

def save_pytecgg_plots(
    df_cal,
    df_veq,
    output_dir,
    station,
    date,
    doy,
):
    """
    Create daily PyTECGg time-series PNG plots.

    Outputs:
        *_PyTECGg_TEC.png
            vTEC observations from all satellites +
            station Vertical Equivalent (VEq)

        *_PyTECGg_STEC.png
            calibrated slant TEC

        *_PyTECGg_VEQ.png
            station-level Vertical Equivalent TEC
    """

    print()
    print("----------------------------------------")
    print("PyTECGg: creating daily plots")
    print("----------------------------------------")

   
# --------------------------------------------------------
    # Convert Polars -> pandas without requiring PyArrow
    # --------------------------------------------------------

    if isinstance(df_cal, pl.DataFrame):
        cal = pd.DataFrame(
            df_cal.to_dicts()
        )
    else:
        cal = pd.DataFrame(
            df_cal
        )

    if isinstance(df_veq, pl.DataFrame):
        veq = pd.DataFrame(
            df_veq.to_dicts()
        )
    else:
        veq = pd.DataFrame(
            df_veq
        )


    # --------------------------------------------------------
    # Convert epoch
    # --------------------------------------------------------

    if "epoch" in cal.columns:

        cal["epoch"] = pd.to_datetime(
            cal["epoch"],
            errors="coerce"
        )

    if "epoch" in veq.columns:

        veq["epoch"] = pd.to_datetime(
            veq["epoch"],
            errors="coerce"
        )


    # ========================================================
    # PLOT 1
    #
    # vTEC satellite observations + VEq
    # Similar concept to pyOASIS daily TEC plot
    # ========================================================

    if (
        "epoch" in cal.columns
        and "vtec" in cal.columns
    ):

        plot_df = cal[
            ["epoch", "sv", "vtec"]
        ].copy()

        plot_df["vtec"] = pd.to_numeric(
            plot_df["vtec"],
            errors="coerce"
        )

        plot_df = plot_df.dropna(
            subset=["epoch", "vtec"]
        )

        if not plot_df.empty:

            fig, ax = plt.subplots(
                figsize=(14, 7)
            )

            # ----------------------------------------------
            # Raw satellite vTEC observations
            # ----------------------------------------------

            ax.scatter(
                plot_df["epoch"],
                plot_df["vtec"],
                s=5,
                alpha=0.30,
                label="Satellite vTEC"
            )


            # ----------------------------------------------
            # VEq
            # ----------------------------------------------

            if (
                "epoch" in veq.columns
                and "veq" in veq.columns
            ):

                veq_plot = veq[
                    ["epoch", "veq"]
                ].copy()

                veq_plot["veq"] = pd.to_numeric(
                    veq_plot["veq"],
                    errors="coerce"
                )

                veq_plot = veq_plot.dropna(
                    subset=["epoch", "veq"]
                )

                # PyTECGg VEq is the same station value
                # repeated for satellites at an epoch.
                # Reduce to one value per epoch.

                veq_plot = (
                    veq_plot
                    .groupby(
                        "epoch",
                        as_index=False
                    )["veq"]
                    .median()
                    .sort_values("epoch")
                )

                if not veq_plot.empty:

                    ax.plot(
                        veq_plot["epoch"],
                        veq_plot["veq"],
                        linewidth=2.0,
                        label="VEq"
                    )


            # ----------------------------------------------
            # Formatting
            # ----------------------------------------------

            ax.set_title(
                f"{station}  |  "
                f"{date:%Y-%m-%d}  |  "
                f"DOY {doy:03d}\n"
                "PyTECGg calibrated vertical TEC"
            )

            ax.set_xlabel(
                "Time (UTC)"
            )

            ax.set_ylabel(
                "TEC [TECU]"
            )

            ax.grid(
                True,
                alpha=0.30
            )

            ax.xaxis.set_major_locator(
                mdates.HourLocator(
                    interval=2
                )
            )

            ax.xaxis.set_major_formatter(
                mdates.DateFormatter(
                    "%H:%M"
                )
            )

            ax.legend()

            fig.autofmt_xdate()

            plt.tight_layout()

            output_png = (
                output_dir
                / (
                    f"{station}_"
                    f"{doy:03d}_"
                    f"{date.year}_"
                    f"PyTECGg_TEC.png"
                )
            )

            plt.savefig(
                output_png,
                dpi=200,
                bbox_inches="tight"
            )

            plt.close(fig)

            print(
                "TEC plot :",
                output_png
            )


    # ========================================================
    # PLOT 2
    #
    # Calibrated slant TEC
    # ========================================================

    if (
        "epoch" in cal.columns
        and "stec" in cal.columns
    ):

        stec_df = cal[
            ["epoch", "sv", "stec"]
        ].copy()

        stec_df["stec"] = pd.to_numeric(
            stec_df["stec"],
            errors="coerce"
        )

        stec_df = stec_df.dropna(
            subset=["epoch", "stec"]
        )

        if not stec_df.empty:

            fig, ax = plt.subplots(
                figsize=(14, 7)
            )

            ax.scatter(
                stec_df["epoch"],
                stec_df["stec"],
                s=5,
                alpha=0.30
            )

            ax.set_title(
                f"{station}  |  "
                f"{date:%Y-%m-%d}  |  "
                f"DOY {doy:03d}\n"
                "PyTECGg calibrated slant TEC"
            )

            ax.set_xlabel(
                "Time (UTC)"
            )

            ax.set_ylabel(
                "sTEC [TECU]"
            )

            ax.grid(
                True,
                alpha=0.30
            )

            ax.xaxis.set_major_locator(
                mdates.HourLocator(
                    interval=2
                )
            )

            ax.xaxis.set_major_formatter(
                mdates.DateFormatter(
                    "%H:%M"
                )
            )

            fig.autofmt_xdate()

            plt.tight_layout()

            output_png = (
                output_dir
                / (
                    f"{station}_"
                    f"{doy:03d}_"
                    f"{date.year}_"
                    f"PyTECGg_STEC.png"
                )
            )

            plt.savefig(
                output_png,
                dpi=200,
                bbox_inches="tight"
            )

            plt.close(fig)

            print(
                "sTEC plot:",
                output_png
            )


    # ========================================================
    # PLOT 3
    #
    # VEq station time series
    # ========================================================

    if (
        "epoch" in veq.columns
        and "veq" in veq.columns
    ):

        veq_plot = veq[
            ["epoch", "veq"]
        ].copy()

        veq_plot["veq"] = pd.to_numeric(
            veq_plot["veq"],
            errors="coerce"
        )

        veq_plot = veq_plot.dropna(
            subset=["epoch", "veq"]
        )

        veq_plot = (
            veq_plot
            .groupby(
                "epoch",
                as_index=False
            )["veq"]
            .median()
            .sort_values("epoch")
        )

        if not veq_plot.empty:

            fig, ax = plt.subplots(
                figsize=(14, 6)
            )

            ax.plot(
                veq_plot["epoch"],
                veq_plot["veq"],
                linewidth=1.8
            )

            ax.set_title(
                f"{station}  |  "
                f"{date:%Y-%m-%d}  |  "
                f"DOY {doy:03d}\n"
                "PyTECGg Vertical Equivalent TEC"
            )

            ax.set_xlabel(
                "Time (UTC)"
            )

            ax.set_ylabel(
                "VEq [TECU]"
            )

            ax.grid(
                True,
                alpha=0.30
            )

            ax.xaxis.set_major_locator(
                mdates.HourLocator(
                    interval=2
                )
            )

            ax.xaxis.set_major_formatter(
                mdates.DateFormatter(
                    "%H:%M"
                )
            )

            fig.autofmt_xdate()

            plt.tight_layout()

            output_png = (
                output_dir
                / (
                    f"{station}_"
                    f"{doy:03d}_"
                    f"{date.year}_"
                    f"PyTECGg_VEQ.png"
                )
            )

            plt.savefig(
                output_png,
                dpi=200,
                bbox_inches="tight"
            )

            plt.close(fig)

            print(
                "VEq plot :",
                output_png
            )


    # Make absolutely sure batch processing
    # does not accumulate matplotlib figures.

    plt.close("all")

def process_day(day_dir: Path, date: datetime) -> bool:
    doy = int(date.strftime("%j"))
    rinex_dir = day_dir / "RINEX"
    nav_dir = day_dir / "PyTECGg_INPUT" / "NAV"
    output_dir = day_dir / "PyTECGg_OUTPUT"

    nav_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    obs_files = find_30s_observation_files(rinex_dir)

    print()
    print("=" * 72)
    print(f"DATE : {date:%Y-%m-%d}")
    print(f"DOY  : {doy:03d}")
    print(f"RINEX: {rinex_dir}")
    print(f"30S FILES: {len(obs_files)}")
    print("=" * 72)

    if not obs_files:
        print("NO 30-second RINEX observation segment available - SKIPPING")
        return False

    station = obs_files[0].name[:4].upper()
    tec_parquet = output_dir / f"{station}_{doy:03d}_{date.year}_PyTECGg_TEC.parquet"
    tec_csv = output_dir / f"{station}_{doy:03d}_{date.year}_PyTECGg_TEC.csv"
    veq_csv = output_dir / f"{station}_{doy:03d}_{date.year}_PyTECGg_VEQ.csv"

    # ============================================================
    # COMPLETED DAY / PLOT BACKFILL
    # ============================================================

    tec_plot_png = (
        output_dir
        / f"{station}_{doy:03d}_{date.year}_PyTECGg_TEC.png"
    )

    stec_plot_png = (
        output_dir
        / f"{station}_{doy:03d}_{date.year}_PyTECGg_STEC.png"
    )

    veq_plot_png = (
        output_dir
        / f"{station}_{doy:03d}_{date.year}_PyTECGg_VEQ.png"
    )

    veq_parquet_existing = (
        output_dir
        / f"{station}_{doy:03d}_{date.year}_PyTECGg_VEQ.parquet"
    )


    if SKIP_COMPLETED and tec_parquet.is_file():

        # --------------------------------------------------------
        # TEC exists, but plots may not exist yet
        # --------------------------------------------------------

        plots_complete = (
            tec_plot_png.is_file()
            and stec_plot_png.is_file()
            and veq_plot_png.is_file()
        )

        if plots_complete:

            print("PyTECGg already completed - SKIPPING")
            print("TEC :", tec_parquet)
            print("PNG : already generated")

            return True


        # --------------------------------------------------------
        # Backfill plots without recalculating PyTECGg
        # --------------------------------------------------------

        print("PyTECGg TEC already exists.")
        print("Missing PNG plot(s) - generating from existing results...")

        if not veq_parquet_existing.is_file():

            print("VEq Parquet is missing:")
            print(" ", veq_parquet_existing)
            print("Full processing will continue.")

        else:

            df_cal_existing = pl.read_parquet(
                tec_parquet
            )

            df_veq_existing = pl.read_parquet(
                veq_parquet_existing
            )

            save_pytecgg_plots(
                df_cal=df_cal_existing,
                df_veq=df_veq_existing,
                output_dir=output_dir,
                station=station,
                date=date,
                doy=doy,
            )

            print()
            print("Plot backfill finished.")
            print("No TEC recalculation required.")

            return True

    print("\n----------------------------------------")
    print("PyTECGg STEP 1: Parse all 30s RINEX segments")
    print("----------------------------------------")

    obs_all, rec_pos, rinex_version, receiver_name = parse_all_observations(obs_files)
    systems = systems_from_observations(obs_all)

    print("Receiver:", receiver_name.upper())
    print("RINEX version:", rinex_version)
    print("GNSS systems:", systems)
    print("ECEF receiver position:", rec_pos)

    print("\n----------------------------------------")
    print("PyTECGg STEP 2: Broadcast NAV (BKG)")
    print("----------------------------------------")

    existing_nav = [p for p in nav_dir.iterdir() if p.is_file() and p.stat().st_size > 0]
    if not existing_nav:
        print(f"Downloading BRDC navigation for {date.year} DOY {doy:03d} ...")
        download_nav_bkg(year=date.year, doys=[doy], output_path=nav_dir)
    else:
        print("Navigation file(s) already present - download skipped")

    nav_files = [p for p in nav_dir.iterdir() if p.is_file() and p.stat().st_size > 0]
    if not nav_files:
        raise RuntimeError(f"No navigation files found after download attempt: {nav_dir}")

    full_nav = merge_navigation_files(nav_files)

    print("\n----------------------------------------")
    print("PyTECGg STEP 3: Context + ephemerides")
    print("----------------------------------------")

    ctx = GNSSContext(
        receiver_pos=rec_pos,
        receiver_name=receiver_name,
        rinex_version=rinex_version,
        h_ipp=H_IPP,
        systems=systems,
    )

    ephem_dict = prepare_ephemeris(full_nav, ctx=ctx)
    if not ephem_dict:
        raise RuntimeError("prepare_ephemeris returned no usable satellites.")
    print("Prepared satellite ephemerides:", len(ephem_dict))

    print("\n----------------------------------------")
    print("PyTECGg STEP 4: Linear combinations")
    print("----------------------------------------")

    df_lc = calculate_linear_combinations(
        obs_all,
        ctx=ctx,
        combinations=["gflc_phase", "gflc_code", "mw"],
        selection_mode=SELECTION_MODE,
    )
    print("Linear-combination rows:", f"{df_lc.height:,}")

    print("\n----------------------------------------")
    print("PyTECGg STEP 5: Arc extraction / levelling")
    print("----------------------------------------")

    df_arcs = extract_arcs(
        df=df_lc,
        ctx=ctx,
        threshold_abs=THRESHOLD_ABS,
        threshold_std=THRESHOLD_STD,
        min_arc_length=MIN_ARC_LENGTH,
        threshold_jump=THRESHOLD_JUMP,
    )
    print("Arc-levelled rows:", f"{df_arcs.height:,}")

    print("\n----------------------------------------")
    print("PyTECGg STEP 6: Satellite geometry")
    print("----------------------------------------")

    df_coords = satellite_coordinates(
        sv_ids=df_arcs["sv"],
        epochs=df_arcs["epoch"],
        ephem_dict=ephem_dict,
    )
    df_geom_input = df_arcs.join(df_coords, on=["sv", "epoch"], how="left")

    print("\n----------------------------------------")
    print("PyTECGg STEP 7: IPP / azimuth / elevation")
    print("----------------------------------------")

    df_final = calculate_ipp(df_geom_input, ctx=ctx, min_elevation=MIN_ELEVATION)
    print(f"Rows after elevation >= {MIN_ELEVATION:.1f} deg:", f"{df_final.height:,}")

    if df_final.is_empty():
        raise RuntimeError("No observations remain after geometry/elevation filtering.")

    geom_parquet = output_dir / f"{station}_{doy:03d}_{date.year}_PyTECGg_GEOMETRY.parquet"
    df_final.write_parquet(geom_parquet)

    print("\n----------------------------------------")
    print("PyTECGg STEP 8: TEC calibration")
    print("----------------------------------------")

    df_cal = calculate_tec(
        df_final,
        ctx=ctx,
        max_polynomial_degree=MAX_POLYNOMIAL_DEGREE,
        batch_size_epochs=BATCH_SIZE_EPOCHS,
    )

    if df_cal.is_empty():
        raise RuntimeError("TEC calibration returned an empty DataFrame.")

    df_cal.write_parquet(tec_parquet)
    df_cal.write_csv(tec_csv)
    print("TEC Parquet:", tec_parquet)
    print("TEC CSV    :", tec_csv)

    print("\n----------------------------------------")
    print("PyTECGg STEP 9: Vertical Equivalent (VEq)")
    print("----------------------------------------")

    df_veq = calculate_vertical_equivalent(
        df_cal,
        ctx=ctx,
        max_polynomial_degree=MAX_POLYNOMIAL_DEGREE,
        batch_size_epochs=BATCH_SIZE_EPOCHS,
    )

    veq_parquet = output_dir / f"{station}_{doy:03d}_{date.year}_PyTECGg_VEQ.parquet"
    df_veq.write_parquet(veq_parquet)
    df_veq.write_csv(veq_csv)
    print("VEq Parquet:", veq_parquet)
    print("VEq CSV    :", veq_csv)
    
    # --------------------------------------------------------
    # STEP 10 - Daily PNG plots
    # --------------------------------------------------------

    print()
    print("----------------------------------------")
    print("PyTECGg STEP 10: Daily time-series plots")
    print("----------------------------------------")

    save_pytecgg_plots(
        df_cal=df_cal,
        df_veq=df_veq,
        output_dir=output_dir,
        station=station,
        date=date,
        doy=doy,
    )

    manifest = output_dir / f"{station}_{doy:03d}_{date.year}_PyTECGg_MANIFEST.txt"
    save_manifest(manifest, date, doy, obs_files, systems, obs_all.height)

    print("\n" + "=" * 72)
    print("FULL PyTECGg DAY FINISHED")
    print("DATE                :", f"{date:%Y-%m-%d}")
    print("DOY                 :", f"{doy:03d}")
    print("30S RINEX SEGMENTS  :", len(obs_files))
    print("TEC ROWS            :", df_cal.height)
    print("VEQ ROWS            :", df_veq.height)
    print("OUTPUT              :", output_dir)
    print("=" * 72)
    return True


def main():
    global BASE, YEAR, SKIP_COMPLETED

    args = parse_args()
    YEAR = args.year
    BASE = args.data_root.resolve()
    SKIP_COMPLETED = not args.force

    selected_date = None
    if args.date is not None:
        try:
            selected_date = datetime.strptime(args.date, "%Y-%m-%d")
        except ValueError as exc:
            raise ValueError("--date must use YYYY-MM-DD format") from exc

        if selected_date.year != YEAR:
            raise ValueError("--date year must match --year")

    year_root = resolve_year_root(BASE, YEAR)
    if not year_root.is_dir():
        raise FileNotFoundError(f"Year folder does not exist: {year_root}")

    print("=" * 72)
    print("KOH2 PyTECGg PROCESSING WORKFLOW")
    print("=" * 72)
    print("Year             :", YEAR)
    print("Year root        :", year_root)
    print("Date only        :", args.date if args.date else "ALL")
    print("Sampling token   :", SAMPLING_TOKEN)
    print("IPP height       :", f"{H_IPP:.1f} m")
    print("Elevation mask   :", f"{MIN_ELEVATION:.1f} deg")
    print("Selection mode   :", SELECTION_MODE)
    print("Minimum arc len. :", MIN_ARC_LENGTH)
    print("Threshold abs    :", THRESHOLD_ABS)
    print("Threshold std    :", THRESHOLD_STD)
    print("Threshold jump   :", THRESHOLD_JUMP)
    print("Polynomial degree:", MAX_POLYNOMIAL_DEGREE)
    print("Batch epochs     :", BATCH_SIZE_EPOCHS)
    print("Force            :", args.force)

    completed = 0
    no_data_or_failed = 0
    examined = 0

    for month_dir in sorted(year_root.iterdir()):
        if not month_dir.is_dir() or not month_dir.name.isdigit():
            continue
        month = int(month_dir.name)
        if not 1 <= month <= 12:
            continue

        for day_dir in sorted(month_dir.iterdir()):
            if not day_dir.is_dir() or not day_dir.name.isdigit():
                continue

            day = int(day_dir.name)
            try:
                date = datetime(YEAR, month, day)
            except ValueError:
                print("Invalid date folder:", day_dir)
                no_data_or_failed += 1
                continue

            if selected_date is not None and date.date() != selected_date.date():
                continue

            examined += 1

            try:
                if process_day(day_dir, date):
                    completed += 1
                else:
                    no_data_or_failed += 1
            except KeyboardInterrupt:
                raise
            except Exception as exc:
                no_data_or_failed += 1
                print("\n" + "!" * 72)
                print("ERROR PROCESSING DAY:", date.strftime("%Y-%m-%d"))
                print(type(exc).__name__ + ":", exc)
                print("The batch will continue with the next day.")
                print("!" * 72)
                traceback.print_exc()

    print("\n" + "=" * 72)
    print("PyTECGg BATCH FINISHED")
    print("YEAR           :", YEAR)
    print("DAYS EXAMINED  :", examined)
    print("COMPLETED      :", completed)
    print("NO DATA/FAILED :", no_data_or_failed)
    print("=" * 72)


if __name__ == "__main__":
    main()
