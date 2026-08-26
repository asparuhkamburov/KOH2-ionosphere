#!/usr/bin/env python3
r"""
Continuous dual-frequency geometry-free phase-fluctuation processor for
high-rate KOH2 RINEX 3 observations.

Scientific role
---------------
The earlier single-frequency SIGMA_PHI_RAD diagnostic was found to contain a
strong nondispersive/common-mode component. This processor instead forms a
dual-frequency geometry-free carrier-phase combination before high-pass
filtering.

Default carrier pairs:
    GPS       L1C - L2W/L2X, selected automatically per satellite
    GLONASS   L1C - L2C

Processing:
    1. Stream all hourly RINEX 3 files in chronological order.
    2. Read carrier phase and RINEX LLI digits.
    3. Convert carrier cycles to metres using constellation-specific
       frequencies (including GLONASS FDMA channel numbers).
    4. Form geometry-free phase L_a*lambda_a - L_b*lambda_b.
    5. Split continuous arcs at:
         - observation gaps,
         - RINEX epoch flags,
         - LLI bit 0 (loss of lock),
         - LLI bit 1 (half-cycle ambiguity),
         - robust geometry-free jumps.
    6. Linearly detrend each continuous GF segment.
    7. Apply a 6th-order zero-phase Butterworth high-pass filter at 0.1 Hz.
    8. Convert the filtered GF metres to an equivalent ionospheric phase
       fluctuation on carrier A.
    9. Compute a 60-s population standard deviation.
   10. Apply 120-s segment-edge QC.
   11. Merge SP3-derived elevation/azimuth from the previous geometry-QC CSV
       and provide an elevation >=30 deg analysis mask.

Scientific status
-----------------
The output SIGMA_PHI_GF_EQUIV_RAD is an EXPERIMENTAL geometry-free
phase-fluctuation proxy. It is not asserted to be a reference-grade ISMR
Phi60 index.

Dependencies
------------
numpy
pandas
scipy

Example
-------
python process_koh2_gf_scintillation_continuous.py ^
  --input-dir "...\RINEX\10Z_01H" ^
  --geometry-csv "...SCINT_1MIN_CONTINUOUS_SP3_GEOMETRY_QC.csv" ^
  --output-dir "...SCINT_GF_CONTINUOUS" ^
  --file-token "_01H_10Z_"
"""

from __future__ import annotations

import argparse
import gzip
import math
import time
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

C_MPS = 299792458.0

DEFAULT_WINDOW_SECONDS = 60.0
DEFAULT_HP_CUTOFF_HZ = 0.1
DEFAULT_FILTER_ORDER = 6
DEFAULT_MIN_COMPLETENESS = 0.80
DEFAULT_GAP_FACTOR = 2.5
DEFAULT_EDGE_GUARD_SECONDS = 120.0
DEFAULT_JUMP_MAD_FACTOR = 12.0
DEFAULT_JUMP_ABS_FLOOR_M = 0.02
DEFAULT_MIN_SAMPLE_RATE_HZ = 5.0
DEFAULT_ELEVATION_MASK_DEG = 30.0

GPS_FREQUENCY_HZ = {
    "1": 1575.42e6,
    "2": 1227.60e6,
    "5": 1176.45e6,
}

# Preferred publication-facing pairs for the present KOH2 receiver.
DEFAULT_GPS_PHASE_A = "L1C"
DEFAULT_GPS_PHASE_B = "auto"
DEFAULT_GLO_PHASE_A = "L1C"
DEFAULT_GLO_PHASE_B = "L2C"


@dataclass
class SeriesBuffer:
    system: str
    sv: str
    phase_a: str
    phase_b: str
    frequency_a_hz: float
    frequency_b_hz: float
    time_unix_s: array = field(default_factory=lambda: array("d"))
    gf_m: array = field(default_factory=lambda: array("d"))
    lli_a: bytearray = field(default_factory=bytearray)
    lli_b: bytearray = field(default_factory=bytearray)
    epoch_flag: bytearray = field(default_factory=bytearray)

    def append(
        self,
        t_unix_s: float,
        gf_m: float,
        lli_a: int,
        lli_b: int,
        epoch_flag: int,
    ) -> None:
        self.time_unix_s.append(float(t_unix_s))
        self.gf_m.append(float(gf_m))
        self.lli_a.append(int(lli_a) & 0xFF)
        self.lli_b.append(int(lli_b) & 0xFF)
        self.epoch_flag.append(int(epoch_flag) & 0xFF)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Continuous L1/L2 geometry-free high-rate phase-fluctuation "
            "processor with LLI, arc-edge, and elevation QC."
        )
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--geometry-csv", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--file-token", default="_01H_10Z_")
    p.add_argument("--gps-phase-a", default=DEFAULT_GPS_PHASE_A)
    p.add_argument(
        "--gps-phase-b",
        default=DEFAULT_GPS_PHASE_B,
        help=(
            "GPS second carrier. Default 'auto' selects L2W or L2X "
            "per satellite from the actually populated RINEX fields."
        ),
    )
    p.add_argument("--glo-phase-a", default=DEFAULT_GLO_PHASE_A)
    p.add_argument("--glo-phase-b", default=DEFAULT_GLO_PHASE_B)
    p.add_argument("--window-seconds", type=float, default=DEFAULT_WINDOW_SECONDS)
    p.add_argument("--hp-cutoff-hz", type=float, default=DEFAULT_HP_CUTOFF_HZ)
    p.add_argument("--filter-order", type=int, default=DEFAULT_FILTER_ORDER)
    p.add_argument(
        "--min-window-completeness",
        type=float,
        default=DEFAULT_MIN_COMPLETENESS,
    )
    p.add_argument("--gap-factor", type=float, default=DEFAULT_GAP_FACTOR)
    p.add_argument(
        "--edge-guard-seconds",
        type=float,
        default=DEFAULT_EDGE_GUARD_SECONDS,
    )
    p.add_argument(
        "--jump-mad-factor",
        type=float,
        default=DEFAULT_JUMP_MAD_FACTOR,
    )
    p.add_argument(
        "--jump-abs-floor-m",
        type=float,
        default=DEFAULT_JUMP_ABS_FLOOR_M,
    )
    p.add_argument(
        "--min-sample-rate-hz",
        type=float,
        default=DEFAULT_MIN_SAMPLE_RATE_HZ,
    )
    p.add_argument(
        "--elevation-mask-deg",
        type=float,
        default=DEFAULT_ELEVATION_MASK_DEG,
    )
    p.add_argument("--parquet", action="store_true")
    return p


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="ignore", newline="")
    return path.open("r", encoding="ascii", errors="ignore", newline="")


def parse_header(stream: TextIO) -> tuple[dict[str, list[str]], dict[str, int]]:
    obs_types: dict[str, list[str]] = {}
    expected: dict[str, int] = {}
    current_system: str | None = None
    glo_channels: dict[str, int] = {}

    while True:
        line = stream.readline()
        if not line:
            raise EOFError("END OF HEADER not found")

        label = line[60:80].strip() if len(line) >= 60 else ""

        if label == "SYS / # / OBS TYPES":
            if line and line[0] != " ":
                current_system = line[0]
                try:
                    expected[current_system] = int(line[3:6])
                except ValueError:
                    expected[current_system] = 0
                obs_types[current_system] = []

            if current_system is not None:
                obs_types[current_system].extend(line[7:60].split())

        elif label == "GLONASS SLOT / FRQ #":
            tokens = line[:60].split()
            i = 1 if tokens and tokens[0].lstrip("+-").isdigit() else 0
            while i + 1 < len(tokens):
                sv = tokens[i]
                if sv.startswith("R"):
                    try:
                        glo_channels[sv] = int(tokens[i + 1])
                    except ValueError:
                        pass
                    i += 2
                else:
                    i += 1

        elif label == "END OF HEADER":
            break

    for system, n in expected.items():
        if n > 0:
            obs_types[system] = obs_types[system][:n]

    return obs_types, glo_channels


def collect_files(input_dir: Path, token: str) -> list[Path]:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    files = [
        p for p in input_dir.iterdir()
        if p.is_file()
        and token.lower() in p.name.lower()
        and (
            p.suffix.lower() in {".rnx", ".obs", ".gz"}
            or ".rnx." in p.name.lower()
        )
    ]
    files.sort(key=lambda p: p.name)

    if not files:
        raise FileNotFoundError(
            f"No RINEX files containing token {token!r} found in {input_dir}"
        )
    return files


def scan_headers(
    files: list[Path],
) -> tuple[dict[str, int], dict[str, list[str]]]:
    """
    Read only headers first. This ensures GLONASS channel numbers are known
    before any observation samples are converted to metres.
    """
    glo_channels: dict[str, int] = {}
    representative_obs: dict[str, list[str]] = {}

    for path in files:
        with open_text(path) as f:
            obs_types, channels = parse_header(f)

        for sv, k in channels.items():
            if sv in glo_channels and glo_channels[sv] != k:
                raise ValueError(
                    f"Conflicting GLONASS channel for {sv}: "
                    f"{glo_channels[sv]} vs {k}"
                )
            glo_channels[sv] = k

        for system, names in obs_types.items():
            if system not in representative_obs:
                representative_obs[system] = names

    return glo_channels, representative_obs


def glonass_frequency_hz(sv: str, phase_obs: str, channels: dict[str, int]) -> float | None:
    k = channels.get(sv)
    if k is None:
        return None

    band = phase_obs[1]
    if band == "1":
        return (1602.0 + 0.5625 * k) * 1e6
    if band == "2":
        return (1246.0 + 0.4375 * k) * 1e6
    if band == "3":
        return 1202.025e6
    return None


def carrier_frequency_hz(
    system: str,
    sv: str,
    phase_obs: str,
    glo_channels: dict[str, int],
) -> float | None:
    band = phase_obs[1]

    if system == "G":
        return GPS_FREQUENCY_HZ.get(band)

    if system == "R":
        return glonass_frequency_hz(sv, phase_obs, glo_channels)

    return None


def parse_epoch_header(line: str) -> tuple[pd.Timestamp, int, int]:
    parts = line[1:].split()
    if len(parts) < 8:
        raise ValueError(f"Malformed RINEX epoch line: {line.rstrip()}")

    year, month, day, hour, minute = map(int, parts[:5])
    sec = float(parts[5])
    flag = int(parts[6])
    nsat = int(parts[7])

    sec_int = int(math.floor(sec))
    microsecond = int(round((sec - sec_int) * 1_000_000))
    if microsecond >= 1_000_000:
        sec_int += 1
        microsecond -= 1_000_000

    epoch = pd.Timestamp(
        year=year,
        month=month,
        day=day,
        hour=hour,
        minute=minute,
        second=sec_int,
        microsecond=microsecond,
    )
    return epoch, flag, nsat


def parse_obs_field(record: str, obs_index: int) -> tuple[float | None, int]:
    start = 3 + 16 * obs_index
    field = record[start:start + 16]
    if len(field) < 14:
        return None, 0

    value_text = field[:14].strip()
    if not value_text:
        return None, 0

    try:
        value = float(value_text)
    except ValueError:
        return None, 0

    lli_text = field[14:15].strip()
    try:
        lli = int(lli_text) if lli_text else 0
    except ValueError:
        lli = 0

    return value, lli


def stream_rinex_files(
    files: list[Path],
    glo_channels: dict[str, int],
    args: argparse.Namespace,
) -> tuple[dict[str, SeriesBuffer], dict[str, int]]:
    series: dict[str, SeriesBuffer] = {}

    stats = {
        "epoch_records": 0,
        "satellite_records": 0,
        "paired_samples": 0,
        "missing_phase_pair_records": 0,
        "missing_glonass_channel_records": 0,
        "gps_l2w_samples": 0,
        "gps_l2x_samples": 0,
        "gps_pair_conflict_records": 0,
    }

    for file_number, path in enumerate(files, 1):
        print(f"[{file_number}/{len(files)}] Parsing {path.name}")

        with open_text(path) as f:
            obs_types, file_channels = parse_header(f)

            # Header scan has already collected channels, but accept additional
            # consistent entries if present.
            for sv, k in file_channels.items():
                if sv in glo_channels and glo_channels[sv] != k:
                    raise ValueError(
                        f"Conflicting GLONASS channel for {sv}: "
                        f"{glo_channels[sv]} vs {k}"
                    )
                glo_channels[sv] = k

            # Build observation indices for the selected carriers.
            g_names = obs_types.get("G", [])
            r_names = obs_types.get("R", [])

            gps_ia = (
                g_names.index(args.gps_phase_a)
                if args.gps_phase_a in g_names
                else None
            )

            gps_b_candidates: list[tuple[str, int]] = []
            if str(args.gps_phase_b).lower() == "auto":
                # In this KOH2 data set the GPS satellites divide naturally
                # between L2W and L2X. Both are on the same L2 frequency.
                for obs in ("L2W", "L2X"):
                    if obs in g_names:
                        gps_b_candidates.append((obs, g_names.index(obs)))
            elif args.gps_phase_b in g_names:
                gps_b_candidates.append(
                    (args.gps_phase_b, g_names.index(args.gps_phase_b))
                )

            glo_ia = (
                r_names.index(args.glo_phase_a)
                if args.glo_phase_a in r_names
                else None
            )
            glo_ib = (
                r_names.index(args.glo_phase_b)
                if args.glo_phase_b in r_names
                else None
            )

            while True:
                line = f.readline()
                if not line:
                    break
                if not line.startswith(">"):
                    continue

                epoch, epoch_flag, nsat = parse_epoch_header(line)
                stats["epoch_records"] += 1
                epoch_unix = epoch.timestamp()

                for _ in range(nsat):
                    record = f.readline()
                    if not record:
                        break

                    stats["satellite_records"] += 1
                    sv = record[:3].strip()
                    if len(sv) < 2:
                        continue

                    system = sv[0]

                    if system == "G":
                        if gps_ia is None or not gps_b_candidates:
                            continue

                        phase_a = args.gps_phase_a
                        value_a, lli_a = parse_obs_field(record, gps_ia)
                        if value_a is None:
                            stats["missing_phase_pair_records"] += 1
                            continue

                        # Once a satellite's L2 tracking type has been learned,
                        # keep it fixed. This avoids mixing ambiguity states.
                        existing_phase_b = (
                            series[sv].phase_b if sv in series else None
                        )

                        selected_b = None
                        value_b = None
                        lli_b = 0

                        if existing_phase_b is not None:
                            for obs_b, ib in gps_b_candidates:
                                if obs_b == existing_phase_b:
                                    vb, lb = parse_obs_field(record, ib)
                                    if vb is not None:
                                        selected_b = obs_b
                                        value_b = vb
                                        lli_b = lb
                                    break
                        else:
                            for obs_b, ib in gps_b_candidates:
                                vb, lb = parse_obs_field(record, ib)
                                if vb is not None:
                                    selected_b = obs_b
                                    value_b = vb
                                    lli_b = lb
                                    break

                        if selected_b is None or value_b is None:
                            # If a previously selected GPS signal disappears
                            # while another L2 observable is populated, do not
                            # silently switch ambiguity states.
                            if existing_phase_b is not None:
                                for obs_b, ib in gps_b_candidates:
                                    if obs_b == existing_phase_b:
                                        continue
                                    vb, _ = parse_obs_field(record, ib)
                                    if vb is not None:
                                        stats["gps_pair_conflict_records"] += 1
                                        break
                            stats["missing_phase_pair_records"] += 1
                            continue

                        phase_b = selected_b
                        if phase_b == "L2W":
                            stats["gps_l2w_samples"] += 1
                        elif phase_b == "L2X":
                            stats["gps_l2x_samples"] += 1

                    elif system == "R":
                        if glo_ia is None or glo_ib is None:
                            continue

                        phase_a = args.glo_phase_a
                        phase_b = args.glo_phase_b
                        value_a, lli_a = parse_obs_field(record, glo_ia)
                        value_b, lli_b = parse_obs_field(record, glo_ib)

                        if value_a is None or value_b is None:
                            stats["missing_phase_pair_records"] += 1
                            continue

                    else:
                        continue

                    fa = carrier_frequency_hz(
                        system, sv, phase_a, glo_channels
                    )
                    fb = carrier_frequency_hz(
                        system, sv, phase_b, glo_channels
                    )

                    if fa is None or fb is None:
                        if system == "R":
                            stats["missing_glonass_channel_records"] += 1
                        continue

                    lambda_a = C_MPS / fa
                    lambda_b = C_MPS / fb
                    gf_m = value_a * lambda_a - value_b * lambda_b

                    if not np.isfinite(gf_m):
                        continue

                    if sv not in series:
                        series[sv] = SeriesBuffer(
                            system=system,
                            sv=sv,
                            phase_a=phase_a,
                            phase_b=phase_b,
                            frequency_a_hz=fa,
                            frequency_b_hz=fb,
                        )
                    else:
                        s = series[sv]
                        if (
                            abs(s.frequency_a_hz - fa) > 1.0
                            or abs(s.frequency_b_hz - fb) > 1.0
                            or s.phase_a != phase_a
                            or s.phase_b != phase_b
                        ):
                            raise ValueError(
                                f"Inconsistent carrier definition for {sv}: "
                                f"{s.phase_a}-{s.phase_b} vs "
                                f"{phase_a}-{phase_b}"
                            )

                    series[sv].append(
                        epoch_unix,
                        gf_m,
                        lli_a,
                        lli_b,
                        epoch_flag,
                    )
                    stats["paired_samples"] += 1

    return series, stats


def robust_jump_flags(
    t: np.ndarray,
    gf: np.ndarray,
    gap_threshold_s: float,
    mad_factor: float,
    abs_floor_m: float,
) -> tuple[np.ndarray, float]:
    n = len(gf)
    flags = np.zeros(n, dtype=bool)

    if n < 3:
        return flags, abs_floor_m

    dt = np.diff(t)
    dgf = np.diff(gf)

    valid = (
        np.isfinite(dt)
        & np.isfinite(dgf)
        & (dt > 0)
        & (dt <= gap_threshold_s)
    )

    if valid.sum() < 10:
        threshold = abs_floor_m
    else:
        x = dgf[valid]
        median = float(np.median(x))
        mad = float(np.median(np.abs(x - median)))
        robust_sigma = 1.4826 * mad
        threshold = max(abs_floor_m, mad_factor * robust_sigma)

    if valid.any():
        x = dgf[valid]
        center = float(np.median(x))
    else:
        center = 0.0

    jump = np.zeros(n - 1, dtype=bool)
    jump[valid] = np.abs(dgf[valid] - center) > threshold
    flags[1:] = jump

    return flags, threshold


def phase_a_equivalent_scale_rad_per_m(fa: float, fb: float) -> float:
    """
    Convert geometry-free carrier range metres to equivalent first-order
    ionospheric carrier-A phase radians.

    Overall sign is immaterial for a standard-deviation metric.
    """
    gf_to_a_m = (fb * fb) / (fa * fa - fb * fb)
    return abs(gf_to_a_m * (2.0 * math.pi * fa / C_MPS))


def process_series(
    s: SeriesBuffer,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, dict[str, float | int]]:
    t = np.frombuffer(s.time_unix_s, dtype=np.float64)
    gf = np.frombuffer(s.gf_m, dtype=np.float64)
    lli_a = np.frombuffer(s.lli_a, dtype=np.uint8)
    lli_b = np.frombuffer(s.lli_b, dtype=np.uint8)
    epoch_flag = np.frombuffer(s.epoch_flag, dtype=np.uint8)

    n = len(t)
    if n < 20:
        return pd.DataFrame(), {
            "samples": n,
            "segments": 0,
            "jump_flags": 0,
            "lli_bit0_samples": 0,
            "lli_bit1_samples": 0,
            "epoch_flag_samples": 0,
        }

    # Remove duplicate epochs, preserving the last record.
    # Duplicate epochs are not expected but must not enter filter design.
    keep = np.ones(n, dtype=bool)
    if n > 1:
        keep[:-1] = t[:-1] != t[1:]

    t = t[keep]
    gf = gf[keep]
    lli_a = lli_a[keep]
    lli_b = lli_b[keep]
    epoch_flag = epoch_flag[keep]

    dt = np.diff(t)
    positive_dt = dt[np.isfinite(dt) & (dt > 0)]
    if len(positive_dt) == 0:
        return pd.DataFrame(), {
            "samples": len(t),
            "segments": 0,
            "jump_flags": 0,
            "lli_bit0_samples": 0,
            "lli_bit1_samples": 0,
            "epoch_flag_samples": 0,
        }

    sample_dt = float(np.median(positive_dt))
    fs = 1.0 / sample_dt

    if fs < args.min_sample_rate_hz:
        print(
            f"  {s.sv}: skipped, sample rate {fs:.3f} Hz "
            f"< {args.min_sample_rate_hz:.3f} Hz"
        )
        return pd.DataFrame(), {
            "samples": len(t),
            "segments": 0,
            "jump_flags": 0,
            "lli_bit0_samples": 0,
            "lli_bit1_samples": 0,
            "epoch_flag_samples": 0,
        }

    if not (0 < args.hp_cutoff_hz < fs / 2.0):
        raise ValueError(
            f"{s.sv}: invalid HP cutoff {args.hp_cutoff_hz} Hz "
            f"for fs={fs:.6f} Hz"
        )

    gap_threshold_s = args.gap_factor * sample_dt

    jump_flags, jump_threshold_m = robust_jump_flags(
        t,
        gf,
        gap_threshold_s,
        args.jump_mad_factor,
        args.jump_abs_floor_m,
    )

    lli_bit0 = ((lli_a & 1) != 0) | ((lli_b & 1) != 0)
    lli_bit1 = ((lli_a & 2) != 0) | ((lli_b & 2) != 0)
    bad_epoch_flag = epoch_flag != 0

    new_segment = np.zeros(len(t), dtype=bool)
    new_segment[0] = True
    if len(t) > 1:
        new_segment[1:] |= (
            (~np.isfinite(dt))
            | (dt <= 0)
            | (dt > gap_threshold_s)
        )
    new_segment |= lli_bit0
    new_segment |= lli_bit1
    new_segment |= bad_epoch_flag
    new_segment |= jump_flags

    segment_id = np.cumsum(new_segment) - 1
    n_segments = int(segment_id[-1] + 1)

    sos = butter(
        args.filter_order,
        args.hp_cutoff_hz,
        btype="highpass",
        fs=fs,
        output="sos",
    )

    hp_gf = np.full(len(t), np.nan, dtype=float)
    edge_seconds = np.full(len(t), np.nan, dtype=float)
    segment_samples = np.zeros(len(t), dtype=np.int32)

    # Conservative minimum length. sosfiltfilt will impose its own padding
    # requirement too; failed short segments are simply left invalid.
    min_filter_samples = max(
        30,
        int(math.ceil(fs * 3.0 / args.hp_cutoff_hz)),
    )

    for seg in range(n_segments):
        idx = np.flatnonzero(segment_id == seg)
        if len(idx) < min_filter_samples:
            continue

        tt = t[idx] - t[idx[0]]
        yy = gf[idx]

        if not np.all(np.isfinite(yy)):
            continue

        # Center time before linear detrending for numerical stability.
        tc = tt - float(np.mean(tt))
        coef = np.polyfit(tc, yy, 1)
        detrended = yy - np.polyval(coef, tc)

        try:
            filtered = sosfiltfilt(sos, detrended)
        except ValueError:
            continue

        hp_gf[idx] = filtered
        segment_samples[idx] = len(idx)
        edge_seconds[idx] = np.minimum(tt, tt[-1] - tt)

    scale_rad_per_m = phase_a_equivalent_scale_rad_per_m(
        s.frequency_a_hz,
        s.frequency_b_hz,
    )
    hp_equiv_rad = hp_gf * scale_rad_per_m

    # Align all windows to UTC minute boundaries.
    t_ns = np.rint(t * 1e9).astype(np.int64)
    epoch = pd.to_datetime(t_ns, unit="ns")
    window_ns = int(round(args.window_seconds * 1e9))
    window_start_ns = (t_ns // window_ns) * window_ns

    # Work with numpy grouping to avoid a 15M-row DataFrame.
    unique_windows, starts = np.unique(window_start_ns, return_index=True)
    ends = np.r_[starts[1:], len(t)]

    expected_samples = int(round(fs * args.window_seconds))
    min_samples = int(
        math.ceil(expected_samples * args.min_window_completeness)
    )

    output_rows: list[dict[str, object]] = []

    for ws_ns, i0, i1 in zip(unique_windows, starts, ends):
        x = hp_equiv_rad[i0:i1]
        x_gf = hp_gf[i0:i1]
        finite = np.isfinite(x)
        n_valid = int(finite.sum())

        sigma = (
            float(np.std(x[finite], ddof=0))
            if n_valid >= min_samples
            else np.nan
        )
        sigma_gf_m = (
            float(np.std(x_gf[np.isfinite(x_gf)], ddof=0))
            if int(np.isfinite(x_gf).sum()) >= min_samples
            else np.nan
        )

        edge_slice = edge_seconds[i0:i1]
        finite_edge = np.isfinite(edge_slice)
        min_edge = (
            float(np.min(edge_slice[finite_edge]))
            if finite_edge.any()
            else np.nan
        )

        n_lli0 = int(lli_bit0[i0:i1].sum())
        n_lli1 = int(lli_bit1[i0:i1].sum())
        n_jump = int(jump_flags[i0:i1].sum())
        n_epoch_flag = int(bad_epoch_flag[i0:i1].sum())

        near_edge = (
            True
            if not np.isfinite(min_edge)
            else min_edge <= args.edge_guard_seconds
        )

        ws = pd.Timestamp(ws_ns)
        wm = ws + pd.to_timedelta(args.window_seconds / 2.0, unit="s")

        output_rows.append({
            "window_start": ws,
            "window_mid": wm,
            "sv": s.sv,
            "system": s.system,
            "phase_a": s.phase_a,
            "phase_b": s.phase_b,
            "frequency_a_hz": s.frequency_a_hz,
            "frequency_b_hz": s.frequency_b_hz,
            "sample_rate_hz": fs,
            "n_samples": int(i1 - i0),
            "n_filtered_samples": n_valid,
            "expected_samples": expected_samples,
            "completeness": n_valid / expected_samples
                if expected_samples else np.nan,
            "sigma_gf_m": sigma_gf_m,
            "sigma_phi_gf_equiv_rad": sigma,
            "gf_to_phase_a_scale_rad_per_m": scale_rad_per_m,
            "min_edge_seconds_gf": min_edge,
            "qc_near_gf_segment_edge": near_edge,
            "lli_loss_lock_count": n_lli0,
            "lli_half_cycle_count": n_lli1,
            "gf_jump_count": n_jump,
            "rinex_epoch_flag_count": n_epoch_flag,
            "gf_jump_threshold_m": jump_threshold_m,
        })

    summary = {
        "samples": len(t),
        "segments": n_segments,
        "jump_flags": int(jump_flags.sum()),
        "lli_bit0_samples": int(lli_bit0.sum()),
        "lli_bit1_samples": int(lli_bit1.sum()),
        "epoch_flag_samples": int(bad_epoch_flag.sum()),
        "sample_rate_hz": fs,
        "jump_threshold_m": jump_threshold_m,
    }

    return pd.DataFrame(output_rows), summary


def load_geometry(
    path: Path,
    elevation_mask_deg: float,
) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Geometry CSV not found: {path}")

    g = pd.read_csv(path)

    required = {
        "window_mid",
        "sv",
        "elevation_deg",
        "azimuth_deg",
        "qc_sp3_geometry_available",
    }
    missing = required - set(g.columns)
    if missing:
        raise ValueError(
            f"Geometry CSV missing required columns: {sorted(missing)}"
        )

    g["window_mid"] = pd.to_datetime(g["window_mid"], errors="coerce")
    g["elevation_deg"] = pd.to_numeric(
        g["elevation_deg"], errors="coerce"
    )
    g["azimuth_deg"] = pd.to_numeric(
        g["azimuth_deg"], errors="coerce"
    )
    g["qc_sp3_geometry_available"] = (
        g["qc_sp3_geometry_available"]
        .astype("boolean")
        .fillna(False)
        .astype(bool)
    )

    # Existing geometry file contains multiple signal rows for the same
    # satellite/minute. SP3 geometry is signal independent; collapse it.
    geom = (
        g[
            [
                "window_mid",
                "sv",
                "elevation_deg",
                "azimuth_deg",
                "qc_sp3_geometry_available",
            ]
        ]
        .drop_duplicates(subset=["window_mid", "sv"])
        .copy()
    )

    geom["qc_elevation_ge_mask"] = (
        geom["qc_sp3_geometry_available"]
        & geom["elevation_deg"].ge(elevation_mask_deg)
    )

    return geom


def add_old_single_frequency_diagnostic(
    out: pd.DataFrame,
    geometry_csv_path: Path,
) -> pd.DataFrame:
    """
    Add the old phase-A single-frequency sigma only as a diagnostic comparison.
    It is not part of qc_gf_analysis.
    """
    g = pd.read_csv(
        geometry_csv_path,
        usecols=lambda c: c in {
            "window_start",
            "sv",
            "system",
            "signal",
            "sigma_phi_rad",
            "s4_cno_proxy",
        },
    )

    if not {
        "window_start", "sv", "signal", "sigma_phi_rad"
    }.issubset(g.columns):
        return out

    g["window_start"] = pd.to_datetime(
        g["window_start"], errors="coerce"
    )

    phase_a_signal = out["phase_a"].str[1:]
    out = out.copy()
    out["_phase_a_signal"] = phase_a_signal

    old = (
        g[
            [
                "window_start",
                "sv",
                "signal",
                "sigma_phi_rad",
                "s4_cno_proxy",
            ]
        ]
        .rename(
            columns={
                "signal": "_phase_a_signal",
                "sigma_phi_rad": "old_single_sigma_phase_a_rad",
                "s4_cno_proxy": "s4_cno_proxy_phase_a",
            }
        )
        .drop_duplicates(
            subset=["window_start", "sv", "_phase_a_signal"]
        )
    )

    out = out.merge(
        old,
        on=["window_start", "sv", "_phase_a_signal"],
        how="left",
        validate="one_to_one",
    )
    return out.drop(columns=["_phase_a_signal"])


def main() -> int:
    args = build_parser().parse_args()

    input_dir = Path(args.input_dir)
    geometry_csv = Path(args.geometry_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = collect_files(input_dir, args.file_token)

    print("=" * 84)
    print("KOH2 CONTINUOUS GEOMETRY-FREE HIGH-RATE PHASE PROCESSOR")
    print("=" * 84)
    print(f"Files: {len(files)}")
    print(
        "GPS pair: "
        f"{args.gps_phase_a} - "
        f"{'L2W/L2X auto per satellite' if str(args.gps_phase_b).lower() == 'auto' else args.gps_phase_b}"
    )
    print(f"GLONASS pair: {args.glo_phase_a} - {args.glo_phase_b}")
    print(f"Window: {args.window_seconds:.1f} s")
    print(f"High-pass cutoff: {args.hp_cutoff_hz:.3f} Hz")
    print(f"Filter order: {args.filter_order}")
    print(f"Edge guard: {args.edge_guard_seconds:.1f} s")
    print(f"Elevation mask: {args.elevation_mask_deg:.1f} deg")
    print()

    t0 = time.perf_counter()

    glo_channels, representative_obs = scan_headers(files)
    print(f"GLONASS frequency channels found: {len(glo_channels)}")
    if "G" in representative_obs:
        print(
            "GPS phase observables: "
            + " ".join(x for x in representative_obs["G"] if x.startswith("L"))
        )
    if "R" in representative_obs:
        print(
            "GLONASS phase observables: "
            + " ".join(x for x in representative_obs["R"] if x.startswith("L"))
        )
    print()

    t_parse0 = time.perf_counter()
    series, parse_stats = stream_rinex_files(
        files,
        glo_channels,
        args,
    )
    parse_seconds = time.perf_counter() - t_parse0

    print()
    print(f"Streaming parse: {parse_seconds:.3f} s")
    print(f"Epoch records: {parse_stats['epoch_records']:,}")
    print(
        f"Satellite records: {parse_stats['satellite_records']:,}"
    )
    print(f"Paired GF samples: {parse_stats['paired_samples']:,}")
    if str(args.gps_phase_b).lower() == "auto":
        print(f"GPS L1C-L2W samples: {parse_stats['gps_l2w_samples']:,}")
        print(f"GPS L1C-L2X samples: {parse_stats['gps_l2x_samples']:,}")
        print(
            "GPS pair-switch conflicts skipped: "
            f"{parse_stats['gps_pair_conflict_records']:,}"
        )
    print(f"Satellite GF series: {len(series)}")
    print(
        "Records missing selected phase pair: "
        f"{parse_stats['missing_phase_pair_records']:,}"
    )
    print(
        "GLONASS records skipped for missing channel: "
        f"{parse_stats['missing_glonass_channel_records']:,}"
    )

    t_metric0 = time.perf_counter()
    frames: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []

    for number, sv in enumerate(sorted(series), 1):
        s = series[sv]
        print(
            f"[GF {number}/{len(series)}] {sv} "
            f"{s.phase_a}-{s.phase_b}: {len(s.time_unix_s):,} samples"
        )
        frame, summary = process_series(s, args)
        if not frame.empty:
            frames.append(frame)

        summary_row = {
            "sv": sv,
            "system": s.system,
            "phase_a": s.phase_a,
            "phase_b": s.phase_b,
            "frequency_a_hz": s.frequency_a_hz,
            "frequency_b_hz": s.frequency_b_hz,
            **summary,
        }
        summaries.append(summary_row)

    if not frames:
        raise RuntimeError("No valid GF minute outputs were produced.")

    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(
        ["window_start", "system", "sv"]
    ).reset_index(drop=True)

    metric_seconds = time.perf_counter() - t_metric0

    geom = load_geometry(
        geometry_csv,
        args.elevation_mask_deg,
    )

    out = out.merge(
        geom,
        on=["window_mid", "sv"],
        how="left",
        validate="many_to_one",
    )

    out["qc_sp3_geometry_available"] = (
        out["qc_sp3_geometry_available"]
        .astype("boolean")
        .fillna(False)
        .astype(bool)
    )
    out["qc_elevation_ge_mask"] = (
        out["qc_elevation_ge_mask"]
        .astype("boolean")
        .fillna(False)
        .astype(bool)
    )

    out["qc_gf_analysis"] = (
        out["sigma_phi_gf_equiv_rad"].notna()
        & (~out["qc_near_gf_segment_edge"])
        & (out["lli_loss_lock_count"] == 0)
        & (out["lli_half_cycle_count"] == 0)
        & (out["gf_jump_count"] == 0)
        & (out["rinex_epoch_flag_count"] == 0)
        & out["qc_sp3_geometry_available"]
        & out["qc_elevation_ge_mask"]
    )

    out = add_old_single_frequency_diagnostic(
        out,
        geometry_csv,
    )

    # Useful diagnostic ratio. Do not interpret where the old metric is tiny.
    out["gf_over_old_single_sigma_ratio"] = np.where(
        pd.to_numeric(
            out["old_single_sigma_phase_a_rad"], errors="coerce"
        ) > 1e-6,
        out["sigma_phi_gf_equiv_rad"]
        / out["old_single_sigma_phase_a_rad"],
        np.nan,
    )

    date_min = out["window_start"].min()
    date_max = out["window_start"].max()
    if pd.isna(date_min) or pd.isna(date_max):
        tag = "UNKNOWN_DATE"
    else:
        tag = pd.Timestamp(date_min).strftime("%Y%m%d")

    out_csv = output_dir / f"KOH2_{tag}_GF_PHASE_1MIN_CONTINUOUS_QC.csv"
    out.to_csv(out_csv, index=False)

    out_parquet = None
    if args.parquet:
        out_parquet = (
            output_dir
            / f"KOH2_{tag}_GF_PHASE_1MIN_CONTINUOUS_QC.parquet"
        )
        out.to_parquet(out_parquet, index=False)

    summary_df = pd.DataFrame(summaries)
    summary_csv = output_dir / f"KOH2_{tag}_GF_SERIES_SUMMARY.csv"
    summary_df.to_csv(summary_csv, index=False)

    valid = out[out["qc_gf_analysis"]].copy()
    valid_g = valid[valid["system"] == "G"]
    valid_r = valid[valid["system"] == "R"]

    manifest = output_dir / f"KOH2_{tag}_GF_PHASE_MANIFEST.txt"
    with manifest.open("w", encoding="utf-8") as fh:
        fh.write("KOH2 CONTINUOUS GEOMETRY-FREE HIGH-RATE PHASE PROCESSING\n")
        fh.write("=" * 68 + "\n")
        fh.write(f"Files processed: {len(files)}\n")
        fh.write(
            "GPS pair: "
            f"{args.gps_phase_a} - "
            f"{'L2W/L2X auto per satellite' if str(args.gps_phase_b).lower() == 'auto' else args.gps_phase_b}\n"
        )
        fh.write(f"GLONASS pair: {args.glo_phase_a} - {args.glo_phase_b}\n")
        fh.write(f"Window: {args.window_seconds:.1f} s\n")
        fh.write(f"High-pass cutoff: {args.hp_cutoff_hz:.3f} Hz\n")
        fh.write(f"Butterworth order: {args.filter_order}\n")
        fh.write(
            f"Minimum window completeness: "
            f"{args.min_window_completeness:.2f}\n"
        )
        fh.write(
            f"GF jump MAD factor: {args.jump_mad_factor:.2f}\n"
        )
        fh.write(
            f"GF jump absolute floor: {args.jump_abs_floor_m:.4f} m\n"
        )
        fh.write(
            f"Segment edge guard: {args.edge_guard_seconds:.1f} s\n"
        )
        fh.write(
            f"Elevation mask: >= {args.elevation_mask_deg:.1f} deg\n"
        )
        fh.write("\n")
        fh.write(f"Epoch records: {parse_stats['epoch_records']}\n")
        fh.write(
            f"Satellite records: {parse_stats['satellite_records']}\n"
        )
        fh.write(
            f"Paired geometry-free samples: "
            f"{parse_stats['paired_samples']}\n"
        )
        fh.write(f"Satellite GF series: {len(series)}\n")
        fh.write(f"Output minute rows: {len(out)}\n")
        fh.write(f"QC-valid GF rows: {len(valid)}\n")
        fh.write(f"QC-valid GPS rows: {len(valid_g)}\n")
        fh.write(f"QC-valid GLONASS rows: {len(valid_r)}\n")
        fh.write("\n")
        fh.write("qc_gf_analysis definition:\n")
        fh.write("  finite sigma_phi_gf_equiv_rad\n")
        fh.write("  AND not within segment-edge guard\n")
        fh.write("  AND no LLI bit-0 event in window\n")
        fh.write("  AND no LLI bit-1 event in window\n")
        fh.write("  AND no robust GF jump in window\n")
        fh.write("  AND no nonzero RINEX epoch flag in window\n")
        fh.write("  AND SP3 geometry available\n")
        fh.write(
            f"  AND elevation_deg >= {args.elevation_mask_deg:.1f}\n"
        )
        fh.write("\n")
        fh.write("Scientific status:\n")
        fh.write(
            "  SIGMA_PHI_GF_EQUIV_RAD is an experimental dual-frequency\n"
        )
        fh.write(
            "  geometry-free phase-fluctuation proxy, not reference-grade Phi60.\n"
        )
        fh.write(
            "  OLD_SINGLE_SIGMA_PHASE_A_RAD is retained only as a diagnostic\n"
        )
        fh.write(
            "  comparison and should not be interpreted as the final phase metric.\n"
        )

    total_seconds = time.perf_counter() - t0

    print()
    print(f"Metric calculation: {metric_seconds:.3f} s")
    print(f"Output minute rows: {len(out):,}")
    print(f"QC-valid GF rows: {len(valid):,}")
    print(f"QC-valid GPS rows: {len(valid_g):,}")
    print(f"QC-valid GLONASS rows: {len(valid_r):,}")

    if len(valid):
        print(
            "QC sigma_phi_gf_equiv_rad median: "
            f"{valid['sigma_phi_gf_equiv_rad'].median():.6f}"
        )
        print(
            "QC sigma_phi_gf_equiv_rad p95: "
            f"{valid['sigma_phi_gf_equiv_rad'].quantile(.95):.6f}"
        )
        print(
            "QC sigma_phi_gf_equiv_rad p99: "
            f"{valid['sigma_phi_gf_equiv_rad'].quantile(.99):.6f}"
        )
        imax = valid["sigma_phi_gf_equiv_rad"].idxmax()
        r = valid.loc[imax]
        print(
            "QC maximum: "
            f"{r['sigma_phi_gf_equiv_rad']:.6f} rad "
            f"at {r['window_start']} {r['sv']} "
            f"elev={r['elevation_deg']:.2f} deg"
        )

    print(f"CSV: {out_csv}")
    if out_parquet is not None:
        print(f"Parquet: {out_parquet}")
    print(f"Series summary: {summary_csv}")
    print(f"Manifest: {manifest}")
    print(f"Total wall: {total_seconds:.3f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
