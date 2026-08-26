#!/usr/bin/env python3
"""
Compute an experimental 60-s geometry-free phase-fluctuation proxy from the
paired output of inspect_koh2_phase_lli_gf.py.

The input must contain:
    epoch
    geometry_free_m
    L1C_lli / L2C_lli (or selected equivalents)
    optionally iono_phase_a_equiv_rad_from_gf

Method
------
1. Split at time gaps, LLI bit-0/bit-1 events, and robust GF jumps.
2. Linearly detrend each continuous geometry-free segment.
3. Apply a 6th-order zero-phase Butterworth high-pass at 0.1 Hz.
4. Convert the high-passed GF metres to equivalent ionospheric carrier-A
   phase radians using the already prepared input scaling when available.
5. Calculate population standard deviation (ddof=0) in 60-s windows.
6. Flag windows within 120 s of a segment edge.

This is an EXPERIMENTAL geometry-free phase-fluctuation proxy, not a
reference-grade ISMR Phi60 product.
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

WINDOW_S = 60.0
HP_HZ = 0.1
ORDER = 6
MIN_COMPLETENESS = 0.80
GAP_FACTOR = 2.5
EDGE_GUARD_S = 120.0
JUMP_MAD_FACTOR = 12.0
JUMP_ABS_FLOOR_M = 0.02


def parser():
    p = argparse.ArgumentParser()
    p.add_argument("--input-csv", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--phase-a", default="L1C")
    p.add_argument("--phase-b", default="L2C")
    p.add_argument("--window-seconds", type=float, default=WINDOW_S)
    p.add_argument("--hp-cutoff-hz", type=float, default=HP_HZ)
    p.add_argument("--filter-order", type=int, default=ORDER)
    p.add_argument("--edge-guard-seconds", type=float, default=EDGE_GUARD_S)
    return p


def robust_jump_flags(gf: np.ndarray) -> tuple[np.ndarray, float]:
    d = np.diff(gf)
    good = np.isfinite(d)
    flags = np.zeros(len(gf), dtype=bool)
    if good.sum() < 10:
        return flags, np.nan

    x = d[good]
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    robust = JUMP_MAD_FACTOR * 1.4826 * mad
    threshold = max(JUMP_ABS_FLOOR_M, robust)
    jump = np.abs(d - med) > threshold
    flags[1:] = np.where(np.isfinite(d), jump, True)
    return flags, threshold


def main():
    a = parser().parse_args()
    inp = Path(a.input_csv)
    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inp)
    df["epoch"] = pd.to_datetime(df["epoch"], errors="coerce")
    df["geometry_free_m"] = pd.to_numeric(df["geometry_free_m"], errors="coerce")
    df = df.dropna(subset=["epoch", "geometry_free_m"]).sort_values("epoch").reset_index(drop=True)

    if len(df) < 20:
        raise RuntimeError("Too few geometry-free samples.")

    dt = df["epoch"].diff().dt.total_seconds().to_numpy()
    positive_dt = dt[np.isfinite(dt) & (dt > 0)]
    sample_dt = float(np.median(positive_dt))
    fs = 1.0 / sample_dt
    gap_threshold = GAP_FACTOR * sample_dt

    if fs < 5.0:
        raise RuntimeError(f"Sample rate too low for this workflow: {fs:.3f} Hz")
    if not (0 < a.hp_cutoff_hz < fs/2):
        raise ValueError("Invalid high-pass cutoff for sample rate.")

    lli0 = np.zeros(len(df), dtype=bool)
    lli1 = np.zeros(len(df), dtype=bool)
    for obs in (a.phase_a, a.phase_b):
        col = f"{obs}_lli"
        if col in df.columns:
            x = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int).to_numpy()
            lli0 |= (x & 1) != 0
            lli1 |= (x & 2) != 0

    jump_flags, jump_threshold_m = robust_jump_flags(df["geometry_free_m"].to_numpy())

    new_segment = np.zeros(len(df), dtype=bool)
    new_segment[0] = True
    new_segment[1:] |= (~np.isfinite(dt[1:])) | (dt[1:] > gap_threshold)
    new_segment |= lli0 | lli1 | jump_flags

    segment_id = np.cumsum(new_segment) - 1
    df["segment_id"] = segment_id
    df["gf_jump_flag"] = jump_flags
    df["lli_loss_lock_any"] = lli0
    df["lli_half_cycle_any"] = lli1

    # If the inspector already supplied the equivalent carrier-A phase, use
    # its exact scale. Derive the metres->radians scale robustly from ratios.
    phase_col = "iono_phase_a_equiv_rad_from_gf"
    if phase_col not in df.columns:
        raise ValueError(
            "Input lacks iono_phase_a_equiv_rad_from_gf. "
            "Use the paired CSV from inspect_koh2_phase_lli_gf.py."
        )

    phase_equiv = pd.to_numeric(df[phase_col], errors="coerce").to_numpy()
    gf = df["geometry_free_m"].to_numpy()
    valid_scale = np.isfinite(phase_equiv) & np.isfinite(gf) & (np.abs(gf) > 1e-6)
    if valid_scale.sum() < 10:
        raise RuntimeError("Cannot determine GF metres-to-equivalent-radians scale.")
    scale_rad_per_m = float(np.median(phase_equiv[valid_scale] / gf[valid_scale]))

    sos = butter(
        a.filter_order,
        a.hp_cutoff_hz,
        btype="highpass",
        fs=fs,
        output="sos",
    )

    hp_gf = np.full(len(df), np.nan)
    edge_s = np.full(len(df), np.nan)
    seg_n = np.zeros(len(df), dtype=int)

    for seg, idx in df.groupby("segment_id").groups.items():
        ii = np.asarray(list(idx), dtype=int)
        if len(ii) < max(30, 3 * (2 * a.filter_order + 1)):
            continue

        t = (df.loc[ii, "epoch"] - df.loc[ii[0], "epoch"]).dt.total_seconds().to_numpy()
        y = df.loc[ii, "geometry_free_m"].to_numpy(dtype=float)

        if not np.all(np.isfinite(y)):
            continue

        # Remove constant + linear trend before the high-pass.
        coef = np.polyfit(t, y, 1)
        yd = y - np.polyval(coef, t)

        try:
            yf = sosfiltfilt(sos, yd)
        except ValueError:
            continue

        hp_gf[ii] = yf
        seg_n[ii] = len(ii)

        dist_start = t
        dist_end = t[-1] - t
        edge_s[ii] = np.minimum(dist_start, dist_end)

    df["gf_highpass_m"] = hp_gf
    df["iono_phase_a_highpass_rad"] = hp_gf * scale_rad_per_m
    df["qc_edge_seconds_gf"] = edge_s
    df["qc_near_gf_segment_edge"] = edge_s <= a.edge_guard_seconds
    df["gf_segment_samples"] = seg_n

    origin = df["epoch"].iloc[0].floor("min")
    sec_from_origin = (df["epoch"] - origin).dt.total_seconds()
    win_id = np.floor(sec_from_origin / a.window_seconds).astype(int)
    df["window_start"] = origin + pd.to_timedelta(win_id * a.window_seconds, unit="s")

    expected_n = int(round(fs * a.window_seconds))
    min_n = int(np.ceil(expected_n * MIN_COMPLETENESS))

    rows = []
    for ws, g in df.groupby("window_start", sort=True):
        x = g["iono_phase_a_highpass_rad"].dropna().to_numpy()
        n = len(x)
        valid = n >= min_n
        rows.append({
            "window_start": ws,
            "window_mid": ws + pd.to_timedelta(a.window_seconds/2, unit="s"),
            "sample_rate_hz": fs,
            "n_samples": n,
            "expected_samples": expected_n,
            "completeness": n/expected_n if expected_n else np.nan,
            "sigma_phi_gf_equiv_rad": float(np.std(x, ddof=0)) if valid else np.nan,
            "max_abs_hp_gf_m": float(np.max(np.abs(g["gf_highpass_m"].dropna())))
                if g["gf_highpass_m"].notna().any() else np.nan,
            "min_edge_seconds_gf": float(g["qc_edge_seconds_gf"].min())
                if g["qc_edge_seconds_gf"].notna().any() else np.nan,
            "qc_near_gf_segment_edge": bool(
                g["qc_near_gf_segment_edge"].fillna(True).any()
            ),
            "lli_loss_lock_count": int(g["lli_loss_lock_any"].sum()),
            "lli_half_cycle_count": int(g["lli_half_cycle_any"].sum()),
            "gf_jump_count": int(g["gf_jump_flag"].sum()),
        })

    out = pd.DataFrame(rows)
    out["qc_gf_analysis"] = (
        out["sigma_phi_gf_equiv_rad"].notna()
        & (~out["qc_near_gf_segment_edge"])
        & (out["lli_loss_lock_count"] == 0)
        & (out["lli_half_cycle_count"] == 0)
        & (out["gf_jump_count"] == 0)
    )

    out_csv = outdir / "GF_PHASE_PROXY_1MIN.csv"
    sample_csv = outdir / "GF_PHASE_PROXY_SAMPLES.csv"
    report = outdir / "GF_PHASE_PROXY_REPORT.txt"

    out.to_csv(out_csv, index=False)
    df.to_csv(sample_csv, index=False)

    good = out[out["qc_gf_analysis"]]
    lines = [
        "KOH2 GEOMETRY-FREE PHASE PROXY DIAGNOSTIC",
        "=" * 52,
        f"Input: {inp}",
        f"Sample rate: {fs:.6f} Hz",
        f"GF->carrier-A phase scale: {scale_rad_per_m:.9f} rad/m",
        f"Robust GF jump threshold: {jump_threshold_m:.6f} m",
        f"Raw GF jump flags: {int(jump_flags.sum())}",
        f"LLI loss-lock samples: {int(lli0.sum())}",
        f"LLI half-cycle samples: {int(lli1.sum())}",
        f"Windows: {len(out)}",
        f"QC-valid windows: {len(good)}",
    ]
    if len(good):
        imax = good["sigma_phi_gf_equiv_rad"].idxmax()
        r = good.loc[imax]
        lines += [
            f"QC median sigma_phi_gf_equiv_rad: {good['sigma_phi_gf_equiv_rad'].median():.6f}",
            f"QC p95 sigma_phi_gf_equiv_rad: {good['sigma_phi_gf_equiv_rad'].quantile(.95):.6f}",
            f"QC maximum sigma_phi_gf_equiv_rad: {r['sigma_phi_gf_equiv_rad']:.6f}",
            f"QC maximum window: {r['window_start']}",
        ]

    # Explicitly report the 14:16 window when present.
    target = pd.Timestamp("2025-01-01 14:16:00")
    hit = out[out["window_start"].eq(target)]
    if len(hit):
        r = hit.iloc[0]
        lines += [
            "",
            "14:16 UTC TARGET",
            f"sigma_phi_gf_equiv_rad: {r['sigma_phi_gf_equiv_rad']:.6f}",
            f"qc_gf_analysis: {bool(r['qc_gf_analysis'])}",
            f"min_edge_seconds_gf: {r['min_edge_seconds_gf']:.1f}",
            f"GF jump count: {int(r['gf_jump_count'])}",
            f"LLI loss-lock count: {int(r['lli_loss_lock_count'])}",
            f"LLI half-cycle count: {int(r['lli_half_cycle_count'])}",
        ]

    lines += [
        "",
        "Scientific status:",
        "Experimental dual-frequency geometry-free phase-fluctuation proxy.",
        "Not asserted to be reference-grade ISMR Phi60.",
    ]
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print("Minute CSV:", out_csv)
    print("Sample CSV:", sample_csv)
    print("Report:", report)


if __name__ == "__main__":
    main()
