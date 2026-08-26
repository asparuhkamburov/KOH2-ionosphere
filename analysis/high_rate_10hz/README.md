# High-rate 10 Hz GNSS phase-fluctuation analysis — EXPERIMENTAL

This module contains the publication-oriented high-rate KOH2 processing workflow
for 10 Hz RINEX 3 observations.

**Status: EXPERIMENTAL.**

The module is included for transparency, reproducibility, and further validation.
It is **not** part of the repository's validated `PASS` claims and does **not**
produce a reference-grade ISMR `Phi60` or calibrated `S4` product.

## Scientific rationale

An initial single-frequency carrier-phase metric (`SIGMA_PHI_RAD`) was implemented
with 60 s windows, linear detrending, and a sixth-order zero-phase Butterworth
high-pass filter at 0.1 Hz. Continuous processing and arc-edge QC removed large
hour-boundary filter artifacts, but cross-frequency diagnostics showed that the
remaining single-frequency phase fluctuations scaled much more closely with
carrier frequency than with the inverse-frequency scaling expected for first-order
ionospheric phase.

Accordingly, the primary experimental phase metric in this module is now based on
a **dual-frequency geometry-free carrier-phase combination** formed in metres
before detrending and high-pass filtering.

The primary output is:

```text
SIGMA_PHI_GF_EQUIV_RAD
```

This is an **experimental geometry-free phase-fluctuation proxy**, expressed as
an equivalent carrier-A ionospheric phase fluctuation in radians. It must not be
described as validated `Phi60`.

The amplitude-side quantity:

```text
S4_CNO_PROXY
```

is derived from RINEX `Sxx` / C/N0-like observables and is **uncalibrated**. It
must not be described as standard or calibrated ISMR `S4`.

## Repository contents

```text
analysis/high_rate_10hz/
├── process_koh2_high_rate_scintillation_continuous_qc120.py
├── add_sp3_geometry_to_scintillation.py
├── process_koh2_gf_scintillation_continuous.py
├── diagnostics/
│   ├── inspect_koh2_phase_lli_gf.py
│   └── analyze_koh2_gf_phase_proxy.py
├── validation/
├── README.md
├── STATUS.md
├── VALIDATION.md
└── requirements-high-rate.txt
```

## Validated representative configuration

Representative day:

```text
2025-01-01 / DOY 001
```

High-rate RINEX sampling:

```text
10 Hz
```

Primary processing defaults:

```text
window length                60 s
high-pass cutoff             0.1 Hz
Butterworth filter order     6
minimum completeness         0.80
continuous arc-edge guard    120 s
GF jump MAD factor           12
GF jump absolute floor       0.020 m
primary elevation mask       30 deg
sensitivity masks            35 and 40 deg
```

Carrier pairs used in the representative KOH2 data:

```text
GPS       L1C - L2W/L2X, selected automatically per satellite
GLONASS   L1C - L2C, with FDMA channel frequencies from RINEX header metadata
```

## Processing sequence

The representative workflow is:

```text
hourly 10 Hz RINEX 3
        |
        v
continuous single-frequency diagnostic
        |
        +--> S4_CNO_PROXY
        |
        +--> diagnostic SIGMA_PHI_RAD
        |
        v
SP3 geometry / elevation
        |
        v
dual-frequency geometry-free processor
        |
        v
LLI + gap + GF-jump segmentation
        |
        v
linear detrending
        |
        v
0.1 Hz zero-phase Butterworth high-pass
        |
        v
60 s SIGMA_PHI_GF_EQUIV_RAD
        |
        v
SP3 + elevation + arc-edge analysis QC
```

The single-frequency stage remains in the public module because it documents the
diagnostic path that revealed the common-mode problem and supplies the validated
geometry table used by the current workflow. Its `SIGMA_PHI_RAD` column should
not be interpreted as the final physical phase-scintillation metric.

## Step 1 — continuous single-frequency diagnostic

```bat
python process_koh2_high_rate_scintillation_continuous_qc120.py ^
  --input-dir "<RINEX_10HZ_DIRECTORY>" ^
  --output-dir "<CONTINUOUS_OUTPUT_DIRECTORY>" ^
  --file-token "_01H_10Z_"
```

The important QC columns include separate phase and C/N0 arc-edge diagnostics.
The 120 s guard was selected after an empirical edge-sensitivity test on the
representative day.

## Step 2 — add SP3 geometry

```bat
python add_sp3_geometry_to_scintillation.py ^
  --input-csv "<CONTINUOUS_1MIN_CSV>" ^
  --sp3-file "<DAY1.SP3>" ^
  --sp3-file "<DAY2.SP3>" ^
  --output-dir "<GEOMETRY_OUTPUT_DIRECTORY>" ^
  --elevation-mask-deg 30
```

The validated KOH2 station ECEF coordinates are supplied as script defaults but
can be overridden from the command line.

The geometry implementation interpolates SP3 ECEF coordinates inside available
orbit coverage and does not extrapolate beyond the SP3 epoch range.

**Time-system caveat:** input window labels and SP3 epoch labels are matched
directly. The user must ensure that the observation and SP3 epoch labels are on
a consistent time scale for the intended analysis.

## Step 3 — continuous geometry-free processing

```bat
python process_koh2_gf_scintillation_continuous.py ^
  --input-dir "<RINEX_10HZ_DIRECTORY>" ^
  --geometry-csv "<CONTINUOUS_SP3_GEOMETRY_QC.csv>" ^
  --output-dir "<GF_OUTPUT_DIRECTORY>" ^
  --file-token "_01H_10Z_"
```

The final analysis mask `qc_gf_analysis` requires:

```text
finite SIGMA_PHI_GF_EQUIV_RAD
AND outside the 120 s segment-edge guard
AND no LLI bit-0 event in the minute
AND no LLI bit-1 event in the minute
AND no robust GF jump in the minute
AND no nonzero RINEX epoch flag in the minute
AND SP3 geometry available
AND elevation >= selected mask
```

The RINEX LLI interpretation follows the RINEX 3 convention: bit 0 indicates a
possible cycle slip/loss of lock; bit 1 indicates possible half-cycle ambiguity;
bit 2 is retained by the diagnostic parser but is not used as the main GF
segmentation criterion.

## Diagnostic scripts

`diagnostics/inspect_koh2_phase_lli_gf.py` reads raw RINEX observation fields,
including LLI/SSI digits, and forms a short-interval geometry-free L1/L2 series.

`diagnostics/analyze_koh2_gf_phase_proxy.py` independently processes that paired
series for a targeted interval. It was used as a cross-check against the full-day
continuous processor.

## Representative 24-hour result

For 2025-01-01:

```text
RINEX files processed             24
epoch records                864,000
satellite records         15,306,248
paired GF samples         14,257,040
satellite GF series               52
output minute rows             23,867
QC-valid GF rows               12,605
QC-valid GPS rows               7,342
QC-valid GLONASS rows           5,263
```

For the primary `elevation >=30 deg` mask:

```text
median SIGMA_PHI_GF_EQUIV_RAD    0.066473 rad
p95                              0.096742 rad
p99                              0.118494 rad
maximum                          0.191216 rad
```

The maximum occurs very close to the 30 deg elevation cutoff, so elevation
sensitivity is essential. The validation package therefore also records 35 and
40 deg sensitivity results.

## Event-level context

A strict `elevation >=40 deg` sensitivity analysis identified an upper-tail
multi-link cluster around 16:50–17:10 UTC. The median did not show a comparable
bulk increase.

A pyOASIS ROTI/DTEC comparison provides **partial satellite-specific support**,
especially for GLONASS R04, where the GF enhancement overlaps elevated ROTI and
DTEC. Other links, including R21, do not show the same agreement. The GPS
pyOASIS products available in the event interval do not include the strongest
GF GPS candidates.

Therefore the event is described conservatively as a geometry-free
phase-fluctuation enhancement associated with a disturbed ionospheric interval,
not as validated scintillation.

## What this module does not establish

This module does not establish:

- calibration against a dedicated scintillation monitor;
- equivalence to receiver-produced `Phi60`;
- equivalence to calibrated ISMR `S4`;
- a universal physical threshold for `SIGMA_PHI_GF_EQUIV_RAD`;
- causal attribution of individual events to geomagnetic forcing;
- complete multi-day or multi-year validation of the high-rate method.

See `VALIDATION.md` and `STATUS.md` for the exact evidence and release wording.

## Dependencies

See `requirements-high-rate.txt`.

`pyarrow` is optional and needed only when Parquet output is requested.

## Data policy

Raw 10 Hz RINEX, SP3 products, large generated CSV/Parquet files, and external
data caches are not included in this source module. Small validation summaries
are retained under `validation/`.
