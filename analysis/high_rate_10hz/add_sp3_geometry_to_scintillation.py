#!/usr/bin/env python3
"""
Add precise-orbit satellite elevation/azimuth and analysis QC to KOH2
high-rate scintillation/proxy output.

Scientific role
---------------
This is a geometry/QC post-processor. It does NOT recompute SIGMA_PHI_RAD or
S4_CNO_PROXY. Satellite ECEF positions are interpolated from SP3 precise
ephemerides at each 60-s window midpoint and transformed to local ENU at KOH2.

Recommended analysis masks produced by this script:
    qc_phase_analysis
        sigma_phi_rad is finite
        AND qc_near_arc_edge_phase is False
        AND SP3 geometry is available
        AND elevation_deg >= threshold (default 30 deg)

    qc_s4_analysis
        s4_cno_proxy is finite
        AND qc_near_arc_edge_cn0 is False
        AND SP3 geometry is available
        AND elevation_deg >= threshold (default 30 deg)

Important
---------
- Input window times are matched directly to SP3 epoch labels. No UTC/GPS
  time-scale conversion is applied. This is appropriate when the RINEX epoch
  labels and SP3 epochs use the same GNSS time labeling, as expected for this
  workflow.
- No SP3 extrapolation is performed. Supply adjacent-day SP3 files if needed
  to cover windows close to the start/end of a daily SP3 file.
- SP3 positions are in km and are converted internally to metres.
- This script does not turn S4_CNO_PROXY into calibrated ISMR S4 and does not
  make SIGMA_PHI_RAD a reference-grade Phi60 product.

Dependencies
------------
numpy
pandas
scipy

Example
-------
python add_sp3_geometry_to_scintillation.py ^
  --input-csv KOH2_20250101_20250101_SCINT_1MIN_CONTINUOUS.csv ^
  --sp3 GFZ0MGXRAP_20250010000_01D_05M_ORB.SP3 ^
  --sp3 GFZ0MGXRAP_20250020000_01D_05M_ORB.SP3 ^
  --output-dir SCINT_CONTINUOUS_GEOMETRY ^
  --elevation-mask-deg 30
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline


# KOH2 approximate ECEF position from the station metadata used in this project.
DEFAULT_STATION_X_M = 1453335.2992
DEFAULT_STATION_Y_M = -2554570.1548
DEFAULT_STATION_Z_M = -5641700.7402
DEFAULT_ELEVATION_MASK_DEG = 30.0
DEFAULT_MAX_SP3_BRACKET_GAP_S = 1200.0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Add SP3-derived elevation/azimuth and QC flags to scintillation CSV."
    )
    p.add_argument("--input-csv", required=True, help="Continuous scintillation/proxy CSV.")
    p.add_argument(
        "--sp3",
        action="append",
        default=[],
        help="SP3/SP3.gz file. May be supplied multiple times.",
    )
    p.add_argument(
        "--sp3-dir",
        default=None,
        help="Optional directory searched for *.sp3, *.SP3, *.sp3.gz, *.SP3.gz.",
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--elevation-mask-deg",
        type=float,
        default=DEFAULT_ELEVATION_MASK_DEG,
    )
    p.add_argument("--station-x", type=float, default=DEFAULT_STATION_X_M)
    p.add_argument("--station-y", type=float, default=DEFAULT_STATION_Y_M)
    p.add_argument("--station-z", type=float, default=DEFAULT_STATION_Z_M)
    p.add_argument(
        "--max-sp3-bracket-gap-s",
        type=float,
        default=DEFAULT_MAX_SP3_BRACKET_GAP_S,
        help="Reject interpolation across a larger SP3 epoch gap.",
    )
    p.add_argument("--parquet", action="store_true")
    return p


def collect_sp3_files(explicit: Iterable[str], sp3_dir: str | None) -> list[Path]:
    files: list[Path] = [Path(x) for x in explicit]

    if sp3_dir:
        d = Path(sp3_dir)
        if not d.is_dir():
            raise FileNotFoundError(f"SP3 directory not found: {d}")
        patterns = ("*.sp3", "*.SP3", "*.sp3.gz", "*.SP3.gz")
        for pattern in patterns:
            files.extend(sorted(d.glob(pattern)))

    # Resolve duplicates while preserving deterministic order.
    unique: list[Path] = []
    seen: set[str] = set()
    for f in files:
        key = str(f.resolve()) if f.exists() else str(f)
        if key not in seen:
            seen.add(key)
            unique.append(f)

    if not unique:
        raise ValueError("No SP3 files supplied. Use --sp3 and/or --sp3-dir.")

    missing = [str(f) for f in unique if not f.is_file()]
    if missing:
        raise FileNotFoundError("SP3 file(s) not found:\n  " + "\n  ".join(missing))

    return unique


def open_text_auto(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("rt", encoding="ascii", errors="replace")


def parse_sp3(files: list[Path]) -> pd.DataFrame:
    """
    Parse SP3 position records.

    Returns columns:
        epoch, sv, x_m, y_m, z_m
    """
    records: list[tuple[pd.Timestamp, str, float, float, float]] = []

    for file_no, path in enumerate(files, 1):
        print(f"[SP3 {file_no}/{len(files)}] Parsing {path.name}")
        current_epoch: pd.Timestamp | None = None

        with open_text_auto(path) as fh:
            for line in fh:
                if not line:
                    continue

                if line.startswith("*"):
                    parts = line[1:].split()
                    if len(parts) < 6:
                        current_epoch = None
                        continue
                    year, month, day, hour, minute = map(int, parts[:5])
                    second = float(parts[5])
                    sec_int = int(math.floor(second))
                    micro = int(round((second - sec_int) * 1_000_000))
                    if micro >= 1_000_000:
                        sec_int += 1
                        micro -= 1_000_000
                    current_epoch = pd.Timestamp(
                        year=year,
                        month=month,
                        day=day,
                        hour=hour,
                        minute=minute,
                        second=sec_int,
                        microsecond=micro,
                    )
                    continue

                if current_epoch is None or not line.startswith("P"):
                    continue

                # SP3 P record: P + satellite id + X Y Z clock
                # Satellite id examples: G01, R24, E11.
                sv = line[1:4].strip()
                if len(sv) < 2:
                    continue

                try:
                    x_km = float(line[4:18])
                    y_km = float(line[18:32])
                    z_km = float(line[32:46])
                except ValueError:
                    # Fallback for unusually spaced SP3 variants.
                    parts = line[1:].split()
                    if len(parts) < 4:
                        continue
                    sv = parts[0]
                    try:
                        x_km, y_km, z_km = map(float, parts[1:4])
                    except ValueError:
                        continue

                xyz = np.array([x_km, y_km, z_km], dtype=float)
                if not np.all(np.isfinite(xyz)):
                    continue

                # SP3 missing-position sentinel is approximately 999999.999999 km.
                if np.any(np.abs(xyz) >= 999000.0):
                    continue

                records.append(
                    (
                        current_epoch,
                        sv,
                        x_km * 1000.0,
                        y_km * 1000.0,
                        z_km * 1000.0,
                    )
                )

    if not records:
        raise ValueError("No usable SP3 position records were parsed.")

    sp3 = pd.DataFrame(records, columns=["epoch", "sv", "x_m", "y_m", "z_m"])
    sp3 = (
        sp3.drop_duplicates(subset=["epoch", "sv"], keep="last")
        .sort_values(["sv", "epoch"])
        .reset_index(drop=True)
    )
    return sp3


def ecef_to_geodetic_wgs84(x: float, y: float, z: float) -> tuple[float, float, float]:
    """
    Convert ECEF metres to geodetic latitude, longitude (radians), height (m)
    on WGS84 using a stable iterative solution.
    """
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)

    lon = math.atan2(y, x)
    p = math.hypot(x, y)

    if p < 1e-12:
        lat = math.copysign(math.pi / 2.0, z)
        b = a * (1.0 - f)
        h = abs(z) - b
        return lat, lon, h

    lat = math.atan2(z, p * (1.0 - e2))

    for _ in range(20):
        sin_lat = math.sin(lat)
        n = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        h = p / math.cos(lat) - n
        lat_new = math.atan2(z, p * (1.0 - e2 * n / (n + h)))
        if abs(lat_new - lat) < 1e-14:
            lat = lat_new
            break
        lat = lat_new

    sin_lat = math.sin(lat)
    n = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    h = p / math.cos(lat) - n
    return lat, lon, h


def ecef_los_to_az_el(
    sat_xyz_m: np.ndarray,
    station_xyz_m: np.ndarray,
    lat_rad: float,
    lon_rad: float,
) -> tuple[np.ndarray, np.ndarray]:
    d = sat_xyz_m - station_xyz_m[None, :]
    dx, dy, dz = d[:, 0], d[:, 1], d[:, 2]

    slon, clon = math.sin(lon_rad), math.cos(lon_rad)
    slat, clat = math.sin(lat_rad), math.cos(lat_rad)

    east = -slon * dx + clon * dy
    north = -slat * clon * dx - slat * slon * dy + clat * dz
    up = clat * clon * dx + clat * slon * dy + slat * dz

    horiz = np.hypot(east, north)
    elev = np.degrees(np.arctan2(up, horiz))
    az = np.degrees(np.arctan2(east, north))
    az = np.mod(az, 360.0)
    return az, elev


def interpolate_satellite(
    sat: pd.DataFrame,
    query_times: pd.DatetimeIndex,
    max_bracket_gap_s: float,
) -> pd.DataFrame:
    """
    Cubic-spline interpolation inside SP3 coverage only.

    Geometry is rejected when the query lies outside the satellite's SP3
    coverage or when the bracketing SP3 epochs have a gap larger than the
    configured maximum. No extrapolation is used.
    """
    sat = sat.sort_values("epoch").drop_duplicates("epoch")
    if len(sat) < 2:
        return pd.DataFrame(index=query_times)

    t0 = sat["epoch"].iloc[0]
    t_sec = (sat["epoch"] - t0).dt.total_seconds().to_numpy(dtype=float)
    q_sec = (query_times - t0).total_seconds().to_numpy(dtype=float)

    xyz = sat[["x_m", "y_m", "z_m"]].to_numpy(dtype=float)

    if len(sat) >= 4:
        fx = CubicSpline(t_sec, xyz[:, 0], extrapolate=False)
        fy = CubicSpline(t_sec, xyz[:, 1], extrapolate=False)
        fz = CubicSpline(t_sec, xyz[:, 2], extrapolate=False)
        q_xyz = np.column_stack([fx(q_sec), fy(q_sec), fz(q_sec)])
        method = "cubic"
    else:
        q_xyz = np.column_stack(
            [
                np.interp(q_sec, t_sec, xyz[:, i], left=np.nan, right=np.nan)
                for i in range(3)
            ]
        )
        method = "linear"

    idx_right = np.searchsorted(t_sec, q_sec, side="left")
    valid_bracket = (idx_right > 0) & (idx_right < len(t_sec))

    # Exact SP3 epoch queries are valid even though searchsorted returns the
    # matching epoch as the right index; bracket using previous + exact.
    exact = np.zeros(len(q_sec), dtype=bool)
    inside_right = idx_right < len(t_sec)
    exact[inside_right] = np.isclose(
        t_sec[idx_right[inside_right]],
        q_sec[inside_right],
        atol=1e-6,
        rtol=0.0,
    )

    bracket_gap = np.full(len(q_sec), np.nan)
    nearest = np.full(len(q_sec), np.nan)

    vr = np.where(valid_bracket)[0]
    if len(vr):
        ir = idx_right[vr]
        left_t = t_sec[ir - 1]
        right_t = t_sec[ir]
        bracket_gap[vr] = right_t - left_t
        nearest[vr] = np.minimum(q_sec[vr] - left_t, right_t - q_sec[vr])

    # First exact epoch has no left bracket but is directly known.
    first_exact = exact & (idx_right == 0)
    if np.any(first_exact):
        bracket_gap[first_exact] = 0.0
        nearest[first_exact] = 0.0
        q_xyz[first_exact] = xyz[0]

    # Other exact epochs are direct known positions too.
    other_exact = exact & (idx_right > 0)
    if np.any(other_exact):
        nearest[other_exact] = 0.0
        q_xyz[other_exact] = xyz[idx_right[other_exact]]

    geometry_ok = np.all(np.isfinite(q_xyz), axis=1)
    geometry_ok &= (
        exact
        | (
            valid_bracket
            & np.isfinite(bracket_gap)
            & (bracket_gap <= max_bracket_gap_s)
        )
    )

    q_xyz[~geometry_ok, :] = np.nan

    out = pd.DataFrame(
        {
            "sat_x_m": q_xyz[:, 0],
            "sat_y_m": q_xyz[:, 1],
            "sat_z_m": q_xyz[:, 2],
            "sp3_nearest_epoch_distance_s": nearest,
            "sp3_bracket_gap_s": bracket_gap,
            "sp3_interpolation_method": method,
            "qc_sp3_geometry_available": geometry_ok,
        },
        index=query_times,
    )
    return out


def bool_series(df: pd.DataFrame, name: str) -> pd.Series:
    if name not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    return df[name].astype("boolean").fillna(False).astype(bool)


def main() -> int:
    args = build_parser().parse_args()

    input_csv = Path(args.input_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_csv.is_file():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    sp3_files = collect_sp3_files(args.sp3, args.sp3_dir)

    print("=" * 78)
    print("KOH2 SP3 ELEVATION / AZIMUTH QC")
    print("=" * 78)
    print(f"Input CSV: {input_csv}")
    print(f"SP3 files: {len(sp3_files)}")
    print(f"Elevation mask: {args.elevation_mask_deg:.1f} deg")
    print()

    df = pd.read_csv(input_csv)
    required = {"window_mid", "sv"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV missing required columns: {sorted(missing)}")

    df["window_mid"] = pd.to_datetime(df["window_mid"], errors="coerce")
    if df["window_mid"].isna().any():
        raise ValueError("Some window_mid values could not be parsed.")

    sp3 = parse_sp3(sp3_files)
    print(f"SP3 position records: {len(sp3):,}")
    print(f"SP3 satellites: {sp3['sv'].nunique()}")
    print(
        "SP3 epoch range: "
        f"{sp3['epoch'].min()} -> {sp3['epoch'].max()}"
    )

    station_xyz = np.array(
        [args.station_x, args.station_y, args.station_z], dtype=float
    )
    lat_rad, lon_rad, h_m = ecef_to_geodetic_wgs84(*station_xyz)

    print(
        "Station geodetic: "
        f"lat={math.degrees(lat_rad):.8f} deg, "
        f"lon={math.degrees(lon_rad):.8f} deg, "
        f"h={h_m:.3f} m"
    )

    # Compute geometry only once for each unique (window_mid, sv), then merge
    # back to all signal/frequency rows.
    keys = (
        df[["window_mid", "sv"]]
        .drop_duplicates()
        .sort_values(["sv", "window_mid"])
        .reset_index(drop=True)
    )

    geom_parts: list[pd.DataFrame] = []
    sp3_groups = {sv: g for sv, g in sp3.groupby("sv", sort=False)}

    for sv, q in keys.groupby("sv", sort=False):
        q_times = pd.DatetimeIndex(q["window_mid"])
        sat = sp3_groups.get(sv)

        if sat is None:
            part = q.copy()
            part["sat_x_m"] = np.nan
            part["sat_y_m"] = np.nan
            part["sat_z_m"] = np.nan
            part["sp3_nearest_epoch_distance_s"] = np.nan
            part["sp3_bracket_gap_s"] = np.nan
            part["sp3_interpolation_method"] = ""
            part["qc_sp3_geometry_available"] = False
            part["azimuth_deg"] = np.nan
            part["elevation_deg"] = np.nan
            geom_parts.append(part)
            continue

        interp = interpolate_satellite(
            sat,
            q_times,
            max_bracket_gap_s=args.max_sp3_bracket_gap_s,
        )

        part = q.copy()
        for col in interp.columns:
            part[col] = interp[col].to_numpy()

        sat_xyz = part[["sat_x_m", "sat_y_m", "sat_z_m"]].to_numpy(dtype=float)
        ok = part["qc_sp3_geometry_available"].to_numpy(dtype=bool)

        az = np.full(len(part), np.nan)
        el = np.full(len(part), np.nan)
        if np.any(ok):
            az_ok, el_ok = ecef_los_to_az_el(
                sat_xyz[ok],
                station_xyz,
                lat_rad,
                lon_rad,
            )
            az[ok] = az_ok
            el[ok] = el_ok

        part["azimuth_deg"] = az
        part["elevation_deg"] = el
        geom_parts.append(part)

    geom = pd.concat(geom_parts, ignore_index=True)

    # Satellite coordinates are useful for validation but unnecessarily large
    # in the publication-facing output, so retain only geometry diagnostics.
    geom_keep = geom[
        [
            "window_mid",
            "sv",
            "azimuth_deg",
            "elevation_deg",
            "sp3_nearest_epoch_distance_s",
            "sp3_bracket_gap_s",
            "sp3_interpolation_method",
            "qc_sp3_geometry_available",
        ]
    ]

    out = df.merge(geom_keep, on=["window_mid", "sv"], how="left", validate="many_to_one")
    out["qc_sp3_geometry_available"] = (
        out["qc_sp3_geometry_available"]
        .astype("boolean")
        .fillna(False)
        .astype(bool)
    )
    out["qc_elevation_ge_30"] = (
        out["qc_sp3_geometry_available"]
        & out["elevation_deg"].ge(args.elevation_mask_deg)
    )

    phase_edge_bad = bool_series(out, "qc_near_arc_edge_phase")
    cn0_edge_bad = bool_series(out, "qc_near_arc_edge_cn0")

    sigma_valid = (
        pd.to_numeric(out.get("sigma_phi_rad"), errors="coerce").notna()
        if "sigma_phi_rad" in out.columns
        else pd.Series(False, index=out.index)
    )
    s4_valid = (
        pd.to_numeric(out.get("s4_cno_proxy"), errors="coerce").notna()
        if "s4_cno_proxy" in out.columns
        else pd.Series(False, index=out.index)
    )

    out["qc_phase_analysis"] = (
        sigma_valid
        & (~phase_edge_bad)
        & out["qc_elevation_ge_30"]
    )
    out["qc_s4_analysis"] = (
        s4_valid
        & (~cn0_edge_bad)
        & out["qc_elevation_ge_30"]
    )

    stem = input_csv.stem
    out_csv = output_dir / f"{stem}_SP3_GEOMETRY_QC.csv"
    out.to_csv(out_csv, index=False)

    out_parquet = None
    if args.parquet:
        out_parquet = output_dir / f"{stem}_SP3_GEOMETRY_QC.parquet"
        out.to_parquet(out_parquet, index=False)

    n_geom = int(out["qc_sp3_geometry_available"].sum())
    n_elev = int(out["qc_elevation_ge_30"].sum())
    n_phase = int(out["qc_phase_analysis"].sum())
    n_s4 = int(out["qc_s4_analysis"].sum())

    available_sv = set(sp3["sv"].unique())
    requested_sv = set(out["sv"].dropna().astype(str).unique())
    missing_sv = sorted(requested_sv - available_sv)

    manifest = output_dir / f"{stem}_SP3_GEOMETRY_QC_MANIFEST.txt"
    with manifest.open("w", encoding="utf-8") as fh:
        fh.write("KOH2 SP3 ELEVATION / AZIMUTH QC\n")
        fh.write("=" * 48 + "\n")
        fh.write(f"Input CSV: {input_csv}\n")
        fh.write(f"Rows: {len(out)}\n")
        fh.write(f"SP3 files: {len(sp3_files)}\n")
        for p in sp3_files:
            fh.write(f"  {p}\n")
        fh.write(f"Station X/Y/Z m: {args.station_x:.4f} {args.station_y:.4f} {args.station_z:.4f}\n")
        fh.write(
            "Station geodetic deg/m: "
            f"{math.degrees(lat_rad):.8f} "
            f"{math.degrees(lon_rad):.8f} "
            f"{h_m:.3f}\n"
        )
        fh.write(f"Elevation mask: >= {args.elevation_mask_deg:.1f} deg\n")
        fh.write(f"Rows with SP3 geometry: {n_geom}\n")
        fh.write(f"Rows elevation >= mask: {n_elev}\n")
        fh.write(f"Rows qc_phase_analysis: {n_phase}\n")
        fh.write(f"Rows qc_s4_analysis: {n_s4}\n")
        fh.write(f"Missing requested satellites in SP3: {len(missing_sv)}\n")
        if missing_sv:
            fh.write("  " + " ".join(missing_sv) + "\n")
        fh.write("\n")
        fh.write("qc_phase_analysis definition:\n")
        fh.write("  finite sigma_phi_rad\n")
        fh.write("  AND qc_near_arc_edge_phase == False\n")
        fh.write("  AND SP3 geometry available\n")
        fh.write(f"  AND elevation_deg >= {args.elevation_mask_deg:.1f}\n")
        fh.write("\n")
        fh.write("qc_s4_analysis definition:\n")
        fh.write("  finite s4_cno_proxy\n")
        fh.write("  AND qc_near_arc_edge_cn0 == False\n")
        fh.write("  AND SP3 geometry available\n")
        fh.write(f"  AND elevation_deg >= {args.elevation_mask_deg:.1f}\n")
        fh.write("\n")
        fh.write("Scientific status:\n")
        fh.write("  Geometry/QC post-processing only.\n")
        fh.write("  S4_CNO_PROXY remains an uncalibrated C/N0-derived proxy.\n")
        fh.write("  SIGMA_PHI_RAD remains an experimental phase-fluctuation estimate.\n")

    print()
    print(f"Rows: {len(out):,}")
    print(f"Rows with SP3 geometry: {n_geom:,}")
    print(f"Rows elevation >= {args.elevation_mask_deg:.1f} deg: {n_elev:,}")
    print(f"Rows qc_phase_analysis: {n_phase:,}")
    print(f"Rows qc_s4_analysis: {n_s4:,}")
    if missing_sv:
        print("Satellites absent from SP3:", " ".join(missing_sv))
    print(f"CSV: {out_csv}")
    if out_parquet:
        print(f"Parquet: {out_parquet}")
    print(f"Manifest: {manifest}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
