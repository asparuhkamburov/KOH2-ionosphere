#!/usr/bin/env python3
"""
Targeted RINEX 3 carrier-phase LLI and geometry-free diagnostic.

Purpose
-------
Inspect a satellite over a short high-rate interval before changing the main
scintillation workflow. The script:
  * parses phase observations (Lxx) including RINEX LLI and SSI digits;
  * decodes LLI bits 0/1/2;
  * pairs two selected carrier phases at common epochs;
  * converts carrier cycles to metres;
  * forms a geometry-free phase combination;
  * reports raw inter-epoch jumps and LLI events.

RINEX 3 LLI meanings used:
  bit 0: loss of lock between previous and current observation;
  bit 1: half-cycle ambiguity/slip possible at current epoch;
  bit 2: BOC/MBOC tracking indicator.

This is a diagnostic script, not a final scintillation-index processor.
"""

from __future__ import annotations
import argparse
import gzip
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd

C_MPS = 299792458.0


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="ignore", newline="")
    return path.open("r", encoding="ascii", errors="ignore", newline="")


def parse_header(stream: TextIO):
    obs_types = {}
    expected = {}
    current = None
    glo_channels = {}

    while True:
        line = stream.readline()
        if not line:
            raise EOFError("END OF HEADER not found")
        label = line[60:80].strip() if len(line) >= 60 else ""

        if label == "SYS / # / OBS TYPES":
            if line and line[0] != " ":
                current = line[0]
                try:
                    expected[current] = int(line[3:6])
                except ValueError:
                    expected[current] = 0
                obs_types[current] = []
            if current is not None:
                obs_types[current].extend(line[7:60].split())

        elif label == "GLONASS SLOT / FRQ #":
            # Parse tokens such as R01  1 R02 -4 ...
            tokens = line[:60].split()
            # first token can be total number of slots
            i = 1 if tokens and tokens[0].lstrip("+-").isdigit() else 0
            while i + 1 < len(tokens):
                sv = tokens[i]
                if sv.startswith("R"):
                    try:
                        glo_channels[sv] = int(tokens[i+1])
                    except ValueError:
                        pass
                    i += 2
                else:
                    i += 1

        elif label == "END OF HEADER":
            break

    for sys, n in expected.items():
        if n > 0:
            obs_types[sys] = obs_types[sys][:n]
    return obs_types, glo_channels


def parse_epoch(line: str):
    f = line[1:].split()
    year, month, day, hour, minute = map(int, f[:5])
    sec = float(f[5])
    flag = int(f[6])
    nsat = int(f[7])
    sec_i = int(math.floor(sec))
    usec = int(round((sec-sec_i)*1e6))
    if usec >= 1_000_000:
        sec_i += 1
        usec -= 1_000_000
    t = pd.Timestamp(datetime(year,month,day,hour,minute,sec_i,usec,
                              tzinfo=timezone.utc)).tz_localize(None)
    return t, flag, nsat


def lli_decode(ch: str):
    if not ch or not ch.strip():
        lli = 0
    else:
        try:
            lli = int(ch)
        except ValueError:
            lli = 0
    return lli, bool(lli & 1), bool(lli & 2), bool(lli & 4)


def carrier_frequency_hz(sv: str, obs: str, glo_channels: dict[str,int]):
    band = obs[1]
    sys = sv[0]

    if sys == "G":
        return {
            "1": 1575.42e6,
            "2": 1227.60e6,
            "5": 1176.45e6,
        }.get(band)

    if sys == "R":
        k = glo_channels.get(sv)
        if k is None:
            return None
        if band == "1":
            return (1602.0 + 0.5625*k)*1e6
        if band == "2":
            return (1246.0 + 0.4375*k)*1e6
        if band == "3":
            return 1202.025e6

    if sys == "E":
        return {
            "1": 1575.42e6,
            "5": 1176.45e6,
            "7": 1207.14e6,
            "8": 1191.795e6,
            "6": 1278.75e6,
        }.get(band)

    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-file", required=True)
    ap.add_argument("--sv", required=True, help="e.g. R02")
    ap.add_argument("--phase-a", default=None, help="e.g. L1C")
    ap.add_argument("--phase-b", default=None, help="e.g. L2C")
    ap.add_argument("--start", default=None, help="ISO time, e.g. 2025-01-01T14:10:00")
    ap.add_argument("--end", default=None, help="ISO time, e.g. 2025-01-01T14:20:00")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    path = Path(args.input_file)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    start = pd.Timestamp(args.start) if args.start else None
    end = pd.Timestamp(args.end) if args.end else None

    rows = []
    with open_text(path) as f:
        obs_types, glo_channels = parse_header(f)

        sys = args.sv[0]
        if sys not in obs_types:
            raise ValueError(f"System {sys} not present in RINEX header")

        phase_names = [x for x in obs_types[sys] if x.startswith("L")]
        if not phase_names:
            raise ValueError(f"No phase observables for system {sys}")

        phase_a = args.phase_a or phase_names[0]
        phase_b = args.phase_b or (phase_names[1] if len(phase_names) > 1 else None)

        selected = [(i, name) for i, name in enumerate(obs_types[sys])
                    if name in {phase_a, phase_b}]

        print("RINEX phase observables:", " ".join(phase_names))
        print("Selected:", phase_a, phase_b)
        if sys == "R":
            print("GLONASS channel:", args.sv, glo_channels.get(args.sv))

        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith(">"):
                continue

            epoch, flag, nsat = parse_epoch(line)

            # Observation records still have to be consumed even outside range.
            for _ in range(nsat):
                rec = f.readline()
                if not rec:
                    break

                if rec[:3].strip() != args.sv:
                    continue

                if start is not None and epoch < start:
                    continue
                if end is not None and epoch > end:
                    continue

                for i, obs in selected:
                    pos = 3 + 16*i
                    field = rec[pos:pos+16]
                    if len(field) < 14:
                        continue
                    txt = field[:14].strip()
                    if not txt:
                        continue
                    try:
                        val = float(txt)
                    except ValueError:
                        continue

                    lli, b0, b1, b2 = lli_decode(field[14:15])
                    ssi_txt = field[15:16].strip()
                    ssi = int(ssi_txt) if ssi_txt.isdigit() else np.nan

                    rows.append({
                        "epoch": epoch,
                        "epoch_flag": flag,
                        "sv": args.sv,
                        "phase_obs": obs,
                        "phase_cycles": val,
                        "lli": lli,
                        "lli_loss_lock_bit0": b0,
                        "lli_half_cycle_bit1": b1,
                        "lli_boc_tracking_bit2": b2,
                        "ssi": ssi,
                    })

    long = pd.DataFrame(rows)
    if long.empty:
        raise RuntimeError("No selected phase observations found in requested interval")

    long = long.sort_values(["epoch","phase_obs"]).reset_index(drop=True)
    long_path = outdir / f"{args.sv}_PHASE_LLI_LONG.csv"
    long.to_csv(long_path, index=False)

    # Wide paired diagnostic
    value = long.pivot_table(index="epoch", columns="phase_obs",
                             values="phase_cycles", aggfunc="first")
    lli = long.pivot_table(index="epoch", columns="phase_obs",
                           values="lli", aggfunc="first")
    ef = long.groupby("epoch")["epoch_flag"].max()

    paired = pd.DataFrame(index=value.index)
    paired["epoch_flag"] = ef

    if phase_a in value.columns:
        paired[f"{phase_a}_cycles"] = value[phase_a]
        paired[f"{phase_a}_lli"] = lli.get(phase_a, 0)
    if phase_b and phase_b in value.columns:
        paired[f"{phase_b}_cycles"] = value[phase_b]
        paired[f"{phase_b}_lli"] = lli.get(phase_b, 0)

    fa = carrier_frequency_hz(args.sv, phase_a, glo_channels)
    fb = carrier_frequency_hz(args.sv, phase_b, glo_channels) if phase_b else None

    if fa is not None and phase_a in value.columns:
        la = C_MPS/fa
        paired[f"{phase_a}_m"] = value[phase_a]*la
        paired[f"d{phase_a}_m"] = paired[f"{phase_a}_m"].diff()

    if fb is not None and phase_b in value.columns:
        lb = C_MPS/fb
        paired[f"{phase_b}_m"] = value[phase_b]*lb
        paired[f"d{phase_b}_m"] = paired[f"{phase_b}_m"].diff()

    if (
        fa is not None and fb is not None
        and phase_a in value.columns and phase_b in value.columns
    ):
        paired["geometry_free_m"] = paired[f"{phase_a}_m"] - paired[f"{phase_b}_m"]
        paired["delta_geometry_free_m"] = paired["geometry_free_m"].diff()

        # Scale GF ionospheric range to equivalent phase fluctuation on carrier A.
        # sigma ignores the overall sign.
        gf_to_a_m = (fb*fb)/(fa*fa - fb*fb)
        paired["iono_phase_a_equiv_rad_from_gf"] = (
            paired["geometry_free_m"] * gf_to_a_m * (2*np.pi*fa/C_MPS)
        )
        paired["delta_iono_phase_a_equiv_rad_from_gf"] = (
            paired["iono_phase_a_equiv_rad_from_gf"].diff()
        )

    paired["lli_bit0_any"] = False
    paired["lli_bit1_any"] = False
    for obs in [phase_a, phase_b]:
        col = f"{obs}_lli"
        if col in paired.columns:
            x = pd.to_numeric(paired[col], errors="coerce").fillna(0).astype(int)
            paired["lli_bit0_any"] |= (x & 1).astype(bool)
            paired["lli_bit1_any"] |= (x & 2).astype(bool)

    paired = paired.reset_index()
    pair_path = outdir / f"{args.sv}_PHASE_LLI_GF_PAIRED.csv"
    paired.to_csv(pair_path, index=False)

    n_lli0 = int(long["lli_loss_lock_bit0"].sum())
    n_lli1 = int(long["lli_half_cycle_bit1"].sum())
    n_lli2 = int(long["lli_boc_tracking_bit2"].sum())

    report = [
        f"RINEX PHASE LLI / GEOMETRY-FREE DIAGNOSTIC — {args.sv}",
        "="*64,
        f"Input: {path}",
        f"Interval: {paired['epoch'].min()} -> {paired['epoch'].max()}",
        f"Selected carriers: {phase_a}, {phase_b}",
        f"Phase samples: {len(long)}",
        f"LLI bit0 (loss lock) samples: {n_lli0}",
        f"LLI bit1 (half-cycle) samples: {n_lli1}",
        f"LLI bit2 (BOC/MBOC) samples: {n_lli2}",
    ]
    if sys == "R":
        report.append(f"GLONASS frequency channel: {glo_channels.get(args.sv)}")
    if fa and fb:
        report += [
            f"Carrier A frequency: {fa/1e6:.6f} MHz",
            f"Carrier B frequency: {fb/1e6:.6f} MHz",
            f"Nondispersive cycles scaling B/A: {fb/fa:.9f}",
            f"Ionospheric phase scaling B/A: {fa/fb:.9f}",
        ]
    if "delta_geometry_free_m" in paired:
        x = paired["delta_geometry_free_m"].dropna()
        if len(x):
            med = float(np.median(x))
            mad = float(np.median(np.abs(x-med)))
            report += [
                f"delta GF median: {med:.6e} m",
                f"delta GF MAD: {mad:.6e} m",
                f"delta GF max abs: {float(np.max(np.abs(x))):.6e} m",
            ]

    report_path = outdir / f"{args.sv}_PHASE_LLI_GF_REPORT.txt"
    report_path.write_text("\n".join(report)+"\n", encoding="utf-8")

    print("\n".join(report))
    print("Long CSV:", long_path)
    print("Paired GF CSV:", pair_path)
    print("Report:", report_path)


if __name__ == "__main__":
    main()
