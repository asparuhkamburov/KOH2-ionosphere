#!/usr/bin/env python3
from __future__ import annotations

"""
KOH2 continuous 10 Hz GNSS phase-fluctuation / C/N0 proxy processor
===================================================================

Purpose
-------
This is the boundary-aware companion to
`process_koh2_high_rate_scintillation_fast.py`.

Unlike the compatibility/reference implementation, this program concatenates
all hourly RINEX 3 observation files first, builds continuous satellite/signal
arcs across file boundaries, and only then applies detrending/filtering.

This specifically addresses artificial filter transients caused by restarting
a zero-phase Butterworth filter independently at every hourly RINEX boundary.

Scientific outputs
------------------
SIGMA_PHI_RAD
    60-s standard deviation of carrier phase after:
      * concatenating hourly files into continuous satellite/signal arcs;
      * splitting at true data gaps;
      * conservative gross cycle-slip splitting;
      * linear detrending;
      * sixth-order Butterworth high-pass filtering at 0.1 Hz.

S4_CNO_PROXY
    Uncorrected amplitude-scintillation proxy derived from RINEX Sxx/C/N0:
      * Sxx dB-Hz -> linear power-like quantity;
      * low-pass trend at 0.1 Hz;
      * normalized intensity standard deviation over 60 s.

IMPORTANT
---------
S4_CNO_PROXY is NOT a calibrated/reference ISMR S4 index.

SIGMA_PHI_RAD is NOT yet claimed as a reference-grade Phi60 product. A
conventional geodetic receiver can retain clock/oscillator, multipath,
cycle-slip, tracking-loop, and firmware effects.

The program therefore writes QC metadata for each minute:
    qc_edge_seconds
        Minimum time, in seconds, from the minute midpoint to the nearest
        start/end of the continuous filtered arc.
    qc_near_arc_edge
        True when that distance is below --edge-guard-seconds.
    qc_phase_segment_samples
        Number of samples in the phase segment that was filtered.
    qc_cn0_segment_samples
        Number of samples in the C/N0 segment that was filtered.

Recommended publication use:
    Exclude qc_near_arc_edge == True from quantitative sigma-phi analysis
    until the guard interval has been sensitivity-tested.

Dependencies
------------
    numpy
    pandas
    scipy

Input
-----
RINEX 3 observation files. The parser keeps only Lxx and Sxx observables.

Example
-------
python process_koh2_high_rate_scintillation_continuous.py ^
  --input-dir "D:\\...\\2025\\01\\01\\RINEX\\10Z_01H" ^
  --output-dir "D:\\...\\SCINT_CONTINUOUS" ^
  --file-token "_01H_10Z_"

No machine-specific runtime paths are embedded in the script.
"""

import argparse
import gzip
import math
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd
from scipy import signal


DEFAULT_WINDOW_SECONDS = 60.0
DEFAULT_HP_CUTOFF_HZ = 0.1
DEFAULT_FILTER_ORDER = 6
DEFAULT_MIN_COMPLETENESS = 0.80
DEFAULT_GAP_FACTOR = 2.5
DEFAULT_MIN_SAMPLE_RATE_HZ = 5.0
DEFAULT_SLIP_ABS_FLOOR_CYCLES = 1.0
DEFAULT_SLIP_MAD_FACTOR = 12.0
DEFAULT_MIN_CN0_DBHZ = 10.0
DEFAULT_MAX_CN0_DBHZ = 70.0

# A zero-phase IIR filter uses samples on both sides of each point.
# The validated sensitivity run uses a 120 s guard around true arc edges.
DEFAULT_EDGE_GUARD_SECONDS = 120.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Continuous/boundary-aware 10 Hz RINEX scintillation/proxy processor."
        )
    )

    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--input-file", type=Path)
    src.add_argument("--input-dir", type=Path)

    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--recursive", action="store_true")
    p.add_argument("--file-token", default="_01H_10Z_")

    p.add_argument("--window-seconds", type=float, default=DEFAULT_WINDOW_SECONDS)
    p.add_argument("--hp-cutoff-hz", type=float, default=DEFAULT_HP_CUTOFF_HZ)
    p.add_argument("--filter-order", type=int, default=DEFAULT_FILTER_ORDER)
    p.add_argument("--min-completeness", type=float, default=DEFAULT_MIN_COMPLETENESS)
    p.add_argument("--gap-factor", type=float, default=DEFAULT_GAP_FACTOR)
    p.add_argument("--min-sample-rate-hz", type=float, default=DEFAULT_MIN_SAMPLE_RATE_HZ)
    p.add_argument("--slip-abs-floor-cycles", type=float, default=DEFAULT_SLIP_ABS_FLOOR_CYCLES)
    p.add_argument("--slip-mad-factor", type=float, default=DEFAULT_SLIP_MAD_FACTOR)
    p.add_argument("--min-cn0-dbhz", type=float, default=DEFAULT_MIN_CN0_DBHZ)
    p.add_argument("--max-cn0-dbhz", type=float, default=DEFAULT_MAX_CN0_DBHZ)
    p.add_argument("--edge-guard-seconds", type=float, default=DEFAULT_EDGE_GUARD_SECONDS)

    p.add_argument("--parquet", action="store_true")
    p.add_argument("--overwrite", action="store_true")

    return p.parse_args()


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="ignore", newline="")
    return path.open("r", encoding="ascii", errors="ignore", newline="")


def parse_rinex3_header(stream: TextIO) -> tuple[dict[str, list[str]], str]:
    obs_types: dict[str, list[str]] = {}
    expected: dict[str, int] = {}
    current: str | None = None
    ss_unit = "UNKNOWN"
    version: float | None = None

    while True:
        line = stream.readline()
        if not line:
            raise EOFError("END OF HEADER not found.")

        label = line[60:80].strip() if len(line) >= 60 else ""

        if label == "RINEX VERSION / TYPE":
            try:
                version = float(line[:9])
            except ValueError:
                version = None

        elif label == "SYS / # / OBS TYPES":
            if line and line[0] != " ":
                current = line[0]
                try:
                    expected[current] = int(line[3:6])
                except ValueError:
                    expected[current] = 0
                obs_types[current] = []

            if current is not None:
                obs_types[current].extend(line[7:60].split())

        elif label == "SIGNAL STRENGTH UNIT":
            ss_unit = line[:20].strip() or "UNKNOWN"

        elif label == "END OF HEADER":
            break

    if version is not None and version < 3:
        raise ValueError(f"RINEX 3.x required; got {version}")

    for system, n in expected.items():
        if n > 0:
            obs_types[system] = obs_types[system][:n]

    if not obs_types:
        raise ValueError("No SYS / # / OBS TYPES records found.")

    return obs_types, ss_unit


def parse_epoch_ns(line: str) -> tuple[int, int, int]:
    f = line[1:].split()
    if len(f) < 8:
        raise ValueError(f"Malformed epoch record: {line.rstrip()}")

    year, month, day, hour, minute = map(int, f[:5])
    second = float(f[5])
    flag = int(f[6])
    n_sat = int(f[7])

    sec_i = int(math.floor(second))
    usec = int(round((second - sec_i) * 1_000_000))
    if usec >= 1_000_000:
        sec_i += 1
        usec -= 1_000_000

    dt = datetime(
        year, month, day, hour, minute, sec_i, usec,
        tzinfo=timezone.utc,
    )
    return int(round(dt.timestamp() * 1e9)), flag, n_sat


def parse_selected_observables(
    path: Path,
) -> tuple[str, int, int, dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]]:
    raw: dict[tuple[str, str], list[list[float | int]]] = defaultdict(lambda: [[], []])
    epoch_count = 0
    sat_record_count = 0

    with open_text(path) as stream:
        obs_types, ss_unit = parse_rinex3_header(stream)

        selected = {
            system: [
                (i, name)
                for i, name in enumerate(names)
                if name.startswith(("L", "S"))
            ]
            for system, names in obs_types.items()
        }

        while True:
            line = stream.readline()
            if not line:
                break
            if not line.startswith(">"):
                continue

            epoch_ns, flag, n_sat = parse_epoch_ns(line)
            epoch_count += 1

            if flag not in (0, 1):
                for _ in range(n_sat):
                    stream.readline()
                continue

            for _ in range(n_sat):
                rec = stream.readline()
                if not rec:
                    break

                sv = rec[:3].strip()
                if len(sv) < 2:
                    continue

                system = sv[0]
                if system not in obs_types:
                    continue

                sat_record_count += 1

                for i, obs_name in selected[system]:
                    start = 3 + 16 * i
                    if start >= len(rec):
                        continue
                    text = rec[start:start + 14].strip()
                    if not text:
                        continue
                    try:
                        value = float(text)
                    except ValueError:
                        continue

                    key = (sv, obs_name)
                    raw[key][0].append(epoch_ns)
                    raw[key][1].append(value)

    series = {
        key: (
            np.asarray(tv[0], dtype=np.int64),
            np.asarray(tv[1], dtype=float),
        )
        for key, tv in raw.items()
    }

    return ss_unit, epoch_count, sat_record_count, series


def collect_files(args: argparse.Namespace) -> list[Path]:
    if args.input_file is not None:
        p = args.input_file.expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(p)
        return [p]

    root = args.input_dir.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)

    it = root.rglob("*") if args.recursive else root.glob("*")
    token = args.file_token.upper()

    files = [
        p for p in it
        if p.is_file()
        and token in p.name.upper()
        and p.name.lower().endswith((".rnx", ".rnx.gz", ".obs", ".obs.gz"))
    ]
    return sorted(files)


def infer_sample_rate_hz(times_ns: np.ndarray) -> float:
    if len(times_ns) < 3:
        return math.nan
    dt = np.diff(times_ns.astype(np.int64)) / 1e9
    dt = dt[np.isfinite(dt) & (dt > 0)]
    if dt.size == 0:
        return math.nan
    med = float(np.median(dt))
    return 1.0 / med if med > 0 else math.nan


def split_at_time_gaps(
    times_ns: np.ndarray,
    nominal_dt: float,
    gap_factor: float,
) -> list[tuple[int, int]]:
    n = len(times_ns)
    if n == 0:
        return []
    if n == 1:
        return [(0, 1)]

    dt = np.diff(times_ns.astype(np.int64)) / 1e9
    cut_after = np.where(
        (~np.isfinite(dt))
        | (dt <= 0)
        | (dt > gap_factor * nominal_dt)
    )[0]

    starts = np.r_[0, cut_after + 1]
    stops = np.r_[cut_after + 1, n]
    return list(zip(starts.tolist(), stops.tolist()))


def split_at_gross_phase_slips(
    phase_cycles: np.ndarray,
    abs_floor: float,
    mad_factor: float,
) -> list[tuple[int, int]]:
    n = len(phase_cycles)
    if n < 5:
        return [(0, n)]

    dd = np.diff(phase_cycles, n=2)
    finite = np.isfinite(dd)
    if finite.sum() < 5:
        return [(0, n)]

    x = dd[finite]
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    robust_sigma = 1.4826 * mad
    threshold = max(abs_floor, mad_factor * robust_sigma)

    idx = np.where(
        np.isfinite(dd)
        & (np.abs(dd - med) > threshold)
    )[0]

    if idx.size == 0:
        return [(0, n)]

    cuts = np.unique(np.clip(idx + 1, 1, n - 1))
    starts = np.r_[0, cuts]
    stops = np.r_[cuts, n]
    return [
        (int(a), int(b))
        for a, b in zip(starts, stops)
        if b > a
    ]


def butter_sos(fs_hz: float, cutoff_hz: float, order: int, btype: str):
    return signal.butter(
        order,
        cutoff_hz,
        btype=btype,
        fs=fs_hz,
        output="sos",
    )


def phase_filtered(
    phase_cycles: np.ndarray,
    fs_hz: float,
    cutoff_hz: float,
    order: int,
) -> np.ndarray | None:
    if len(phase_cycles) < max(60, int(10 * fs_hz)):
        return None
    if not np.all(np.isfinite(phase_cycles)):
        return None

    centered = phase_cycles - phase_cycles[0]
    centered = signal.detrend(centered, type="linear")
    radians = centered * (2 * np.pi)

    try:
        return signal.sosfiltfilt(
            butter_sos(fs_hz, cutoff_hz, order, "highpass"),
            radians,
        )
    except ValueError:
        return None


def normalized_cn0(
    cn0: np.ndarray,
    fs_hz: float,
    cutoff_hz: float,
    order: int,
) -> np.ndarray | None:
    if len(cn0) < max(60, int(10 * fs_hz)):
        return None
    if not np.all(np.isfinite(cn0)):
        return None

    intensity = np.power(10.0, cn0 / 10.0)
    try:
        trend = signal.sosfiltfilt(
            butter_sos(fs_hz, cutoff_hz, order, "lowpass"),
            intensity,
        )
    except ValueError:
        return None

    good = np.isfinite(trend) & (trend > 0) & np.isfinite(intensity)
    if good.sum() < 10:
        return None

    out = np.full(intensity.shape, np.nan, dtype=float)
    out[good] = intensity[good] / trend[good]
    return out


def min_samples(args: argparse.Namespace, fs_hz: float) -> int:
    expected = max(1, int(round(args.window_seconds * fs_hz)))
    return int(math.floor(expected * args.min_completeness))


def arc_edge_distance_seconds(
    minute_mid_ns: int,
    arc_start_ns: int,
    arc_end_ns: int,
) -> float:
    return min(
        (minute_mid_ns - arc_start_ns) / 1e9,
        (arc_end_ns - minute_mid_ns) / 1e9,
    )


def phase_rows(
    station: str,
    sv: str,
    suffix: str,
    times_ns: np.ndarray,
    phase_cycles: np.ndarray,
    fs_hz: float,
    args: argparse.Namespace,
) -> pd.DataFrame:
    rows = []
    minimum = min_samples(args, fs_hz)
    nominal_dt = 1.0 / fs_hz

    for ga, gb in split_at_time_gaps(times_ns, nominal_dt, args.gap_factor):
        t_run = times_ns[ga:gb]
        p_run = phase_cycles[ga:gb]
        if len(p_run) < minimum:
            continue

        for sa, sb in split_at_gross_phase_slips(
            p_run,
            args.slip_abs_floor_cycles,
            args.slip_mad_factor,
        ):
            t = t_run[sa:sb]
            p = p_run[sa:sb]
            if len(p) < minimum:
                continue

            filtered = phase_filtered(
                p, fs_hz, args.hp_cutoff_hz, args.filter_order
            )
            if filtered is None:
                continue

            idx = pd.DatetimeIndex(t.astype("datetime64[ns]"))
            minute = idx.floor("min")
            tmp = pd.DataFrame({"minute": minute, "phi": filtered})

            arc_start_ns = int(t[0])
            arc_end_ns = int(t[-1])
            seg_samples = len(t)

            for m, g in tmp.groupby("minute", sort=True):
                vals = g["phi"].to_numpy(float)
                vals = vals[np.isfinite(vals)]
                if vals.size < minimum:
                    continue

                mid = m + pd.Timedelta(seconds=args.window_seconds / 2.0)
                mid_ns = int(mid.to_datetime64().astype("datetime64[ns]").astype(np.int64))
                edge_s = arc_edge_distance_seconds(
                    mid_ns, arc_start_ns, arc_end_ns
                )

                rows.append({
                    "window_start": m,
                    "window_mid": mid,
                    "station": station,
                    "sv": sv,
                    "system": sv[:1],
                    "signal": suffix,
                    "phase_obs": "L" + suffix,
                    "sample_rate_hz": fs_hz,
                    "n_phase": int(vals.size),
                    "sigma_phi_rad": float(np.std(vals, ddof=0)),
                    "qc_phase_segment_samples": int(seg_samples),
                    "qc_edge_seconds_phase": float(edge_s),
                    "qc_near_arc_edge_phase": bool(
                        edge_s <= args.edge_guard_seconds
                    ),
                })

    if not rows:
        return pd.DataFrame()

    return (
        pd.DataFrame(rows)
        .sort_values("n_phase", ascending=False)
        .drop_duplicates(
            subset=["window_start", "sv", "signal"],
            keep="first",
        )
        .sort_values(["window_start", "sv", "signal"])
        .reset_index(drop=True)
    )


def s4_rows(
    station: str,
    sv: str,
    suffix: str,
    times_ns: np.ndarray,
    cn0: np.ndarray,
    fs_hz: float,
    args: argparse.Namespace,
) -> pd.DataFrame:
    rows = []
    minimum = min_samples(args, fs_hz)
    nominal_dt = 1.0 / fs_hz

    for a, b in split_at_time_gaps(times_ns, nominal_dt, args.gap_factor):
        t = times_ns[a:b]
        s = cn0[a:b]
        if len(s) < minimum:
            continue

        norm = normalized_cn0(
            s, fs_hz, args.hp_cutoff_hz, args.filter_order
        )
        if norm is None:
            continue

        idx = pd.DatetimeIndex(t.astype("datetime64[ns]"))
        minute = idx.floor("min")
        tmp = pd.DataFrame(
            {"minute": minute, "norm": norm, "cn0": s}
        )

        arc_start_ns = int(t[0])
        arc_end_ns = int(t[-1])
        seg_samples = len(t)

        for m, g in tmp.groupby("minute", sort=True):
            n = g["norm"].to_numpy(float)
            c = g["cn0"].to_numpy(float)

            good = (
                np.isfinite(n)
                & np.isfinite(c)
                & (c >= args.min_cn0_dbhz)
                & (c <= args.max_cn0_dbhz)
            )
            n = n[good]
            c = c[good]
            if n.size < minimum:
                continue

            mean_i = float(np.mean(n))
            if not np.isfinite(mean_i) or mean_i <= 0:
                continue

            variance = max(
                0.0,
                float(np.mean(n * n)) - mean_i * mean_i,
            )
            s4 = math.sqrt(variance) / mean_i

            mid = m + pd.Timedelta(seconds=args.window_seconds / 2.0)
            mid_ns = int(mid.to_datetime64().astype("datetime64[ns]").astype(np.int64))
            edge_s = arc_edge_distance_seconds(
                mid_ns, arc_start_ns, arc_end_ns
            )

            rows.append({
                "window_start": m,
                "window_mid": mid,
                "station": station,
                "sv": sv,
                "system": sv[:1],
                "signal": suffix,
                "strength_obs": "S" + suffix,
                "sample_rate_hz": fs_hz,
                "n_cn0": int(n.size),
                "mean_cn0_dbhz": float(np.mean(c)),
                "std_cn0_dbhz": float(np.std(c, ddof=0)),
                "s4_cno_proxy": float(s4),
                "qc_cn0_segment_samples": int(seg_samples),
                "qc_edge_seconds_cn0": float(edge_s),
                "qc_near_arc_edge_cn0": bool(
                    edge_s <= args.edge_guard_seconds
                ),
            })

    if not rows:
        return pd.DataFrame()

    return (
        pd.DataFrame(rows)
        .drop_duplicates(
            subset=["window_start", "sv", "signal"],
            keep="first",
        )
        .sort_values(["window_start", "sv", "signal"])
        .reset_index(drop=True)
    )


def main() -> int:
    args = parse_args()

    if not (0 < args.min_completeness <= 1):
        raise ValueError("--min-completeness must be in (0, 1].")
    if args.min_cn0_dbhz >= args.max_cn0_dbhz:
        raise ValueError("Invalid C/N0 range.")
    if args.edge_guard_seconds < 0:
        raise ValueError("--edge-guard-seconds must be >= 0.")

    files = collect_files(args)
    if not files:
        print("No matching files found.")
        return 2

    outdir = args.output_dir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    station = files[0].name[:4].upper() if len(files[0].name) >= 4 else "UNKN"

    print("=" * 78)
    print("KOH2 CONTINUOUS HIGH-RATE SCINTILLATION/PROXY PROCESSOR")
    print("=" * 78)
    print("Files:", len(files))
    print("Window:", args.window_seconds, "s")
    print("Cutoff:", args.hp_cutoff_hz, "Hz")
    print("Filter order:", args.filter_order)
    print("Edge guard:", args.edge_guard_seconds, "s")

    t0 = time.perf_counter()

    merged: dict[tuple[str, str], list[tuple[np.ndarray, np.ndarray]]] = defaultdict(list)
    units = set()
    total_epochs = 0
    total_sat_records = 0

    for i, path in enumerate(files, 1):
        print(f"[{i}/{len(files)}] Parsing {path.name}")
        ss_unit, epochs, sat_records, series = parse_selected_observables(path)
        units.add(ss_unit)
        total_epochs += epochs
        total_sat_records += sat_records

        for key, pair in series.items():
            merged[key].append(pair)

    parse_s = time.perf_counter() - t0

    # Concatenate and sort each observable across all files.
    continuous: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    duplicate_epoch_count = 0

    for key, chunks in merged.items():
        times = np.concatenate([c[0] for c in chunks])
        values = np.concatenate([c[1] for c in chunks])

        order = np.argsort(times, kind="mergesort")
        times = times[order]
        values = values[order]

        if len(times) > 1:
            keep = np.r_[True, np.diff(times) != 0]
            duplicate_epoch_count += int((~keep).sum())
            times = times[keep]
            values = values[keep]

        continuous[key] = (times, values)

    print(f"Streaming parse + merge: {parse_s:.3f} s")
    print("Epoch records:", total_epochs)
    print("Satellite records:", total_sat_records)
    print("Continuous L/S series:", len(continuous))
    print("Duplicate epochs removed:", duplicate_epoch_count)

    t1 = time.perf_counter()
    satellites = sorted({sv for sv, _ in continuous})
    suffixes = sorted({
        obs[1:]
        for _, obs in continuous
        if obs.startswith(("L", "S"))
    })

    pieces = []

    for sv in satellites:
        for suffix in suffixes:
            phase_df = pd.DataFrame()
            s4_df = pd.DataFrame()

            pkey = (sv, "L" + suffix)
            skey = (sv, "S" + suffix)

            if pkey in continuous:
                t, y = continuous[pkey]
                fs = infer_sample_rate_hz(t)
                if np.isfinite(fs) and fs >= args.min_sample_rate_hz:
                    phase_df = phase_rows(
                        station, sv, suffix, t, y, fs, args
                    )

            if skey in continuous:
                t, y = continuous[skey]
                valid = (
                    np.isfinite(y)
                    & (y >= args.min_cn0_dbhz)
                    & (y <= args.max_cn0_dbhz)
                )
                tv = t[valid]
                yv = y[valid]
                if len(tv) >= 3:
                    fs = infer_sample_rate_hz(tv)
                    if np.isfinite(fs) and fs >= args.min_sample_rate_hz:
                        s4_df = s4_rows(
                            station, sv, suffix, tv, yv, fs, args
                        )

            keys = [
                "window_start", "window_mid", "station",
                "sv", "system", "signal",
            ]

            if not phase_df.empty and not s4_df.empty:
                m = phase_df.merge(
                    s4_df, on=keys, how="outer", suffixes=("", "_s4")
                )
                if "sample_rate_hz_s4" in m.columns:
                    m["sample_rate_hz"] = m["sample_rate_hz"].fillna(
                        m["sample_rate_hz_s4"]
                    )
                    m.drop(columns=["sample_rate_hz_s4"], inplace=True)
                pieces.append(m)
            elif not phase_df.empty:
                pieces.append(phase_df)
            elif not s4_df.empty:
                pieces.append(s4_df)

    result = (
        pd.concat(pieces, ignore_index=True, sort=False)
        if pieces else pd.DataFrame()
    )

    if not result.empty:
        result["signal_strength_unit_header"] = (
            ",".join(sorted(units)) if units else "UNKNOWN"
        )

        phase_flag = (
            result.get(
                "qc_near_arc_edge_phase",
                pd.Series(
                    False,
                    index=result.index,
                    dtype="boolean",
                ),
            )
            .astype("boolean")
            .fillna(False)
            .astype(bool)
        )

        cn0_flag = (
            result.get(
                "qc_near_arc_edge_cn0",
                pd.Series(
                    False,
                    index=result.index,
                    dtype="boolean",
                ),
            )
            .astype("boolean")
            .fillna(False)
            .astype(bool)
        )

        result["qc_near_arc_edge"] = phase_flag | cn0_flag

        result = (
            result.sort_values(["window_start", "sv", "signal"])
            .drop_duplicates(
                subset=["window_start", "sv", "signal"],
                keep="first",
            )
            .reset_index(drop=True)
        )

    metric_s = time.perf_counter() - t1

    if len(files) == 1:
        stem = files[0].stem
        if stem.lower().endswith(".rnx"):
            stem = stem[:-4]
        base = stem + "_SCINT_1MIN_CONTINUOUS"
    else:
        # Use UTC date span visible in outputs, rather than guessing from filenames.
        if not result.empty:
            d0 = pd.Timestamp(result["window_start"].min()).strftime("%Y%m%d")
            d1 = pd.Timestamp(result["window_start"].max()).strftime("%Y%m%d")
            base = f"{station}_{d0}_{d1}_SCINT_1MIN_CONTINUOUS"
        else:
            base = f"{station}_SCINT_1MIN_CONTINUOUS"

    csv_path = outdir / f"{base}.csv"
    result.to_csv(csv_path, index=False)

    if args.parquet:
        try:
            result.to_parquet(outdir / f"{base}.parquet", index=False)
        except Exception as exc:
            print("Parquet not written:", type(exc).__name__, exc)

    manifest = outdir / f"{base}_MANIFEST.txt"
    sigma_count = int(result["sigma_phi_rad"].notna().sum()) if "sigma_phi_rad" in result else 0
    s4_count = int(result["s4_cno_proxy"].notna().sum()) if "s4_cno_proxy" in result else 0
    edge_count = int(result["qc_near_arc_edge"].sum()) if "qc_near_arc_edge" in result else 0

    manifest.write_text(
        "\n".join([
            f"Station: {station}",
            f"Files processed: {len(files)}",
            f"Output rows: {len(result)}",
            f"Valid sigma_phi_rad rows: {sigma_count}",
            f"Valid s4_cno_proxy rows: {s4_count}",
            f"Rows flagged qc_near_arc_edge: {edge_count}",
            f"Window length: {args.window_seconds:.1f} s",
            f"High-pass cutoff for sigma_phi: {args.hp_cutoff_hz:.3f} Hz",
            f"Butterworth order: {args.filter_order}",
            f"Minimum window completeness: {args.min_completeness:.2f}",
            f"Edge guard: {args.edge_guard_seconds:.1f} s",
            "",
            "IMPORTANT:",
            "S4_CNO_PROXY is an uncorrected proxy derived from RINEX Sxx/CN0.",
            "It is not asserted to be a calibrated ISMR S4 index.",
            "SIGMA_PHI_RAD is not yet asserted to be a reference-grade Phi60 index.",
            "For quantitative analysis, consider excluding qc_near_arc_edge == True.",
            "",
            "Files:",
            *[f"  {p}" for p in files],
        ]) + "\n",
        encoding="utf-8",
    )

    print(f"Metric calculation: {metric_s:.3f} s")
    print("Output rows:", len(result))
    print("Valid sigma_phi_rad:", sigma_count)
    print("Valid s4_cno_proxy:", s4_count)
    print("QC near-edge rows:", edge_count)
    print("CSV:", csv_path)
    print("Manifest:", manifest)
    print(f"Total wall: {time.perf_counter() - t0:.3f} s")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
