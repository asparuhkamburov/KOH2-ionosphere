# KOH2 PyTECGg processing

Publication-oriented wrapper for the PyTECGg workflow used to derive calibrated
slant TEC (sTEC), vertical TEC (vTEC), and station Vertical Equivalent TEC
(VEq) from KOH2 RINEX observations.

## Scientific processing defaults

The publication script preserves the operational processing settings:

```text
Input sampling token : _30S_
IPP height           : 350000 m
Elevation mask       : 30 deg
Selection mode       : quality
Minimum arc length   : 30
Threshold abs        : 5.0
Threshold std        : 5.0
Threshold jump       : 10.0
Polynomial degree    : 3
Batch size           : 30 epochs
```

## Daily processing sequence

```text
30 s RINEX observation segment(s)
    -> parse and concatenate
    -> remove duplicate epoch/SV/observable records
    -> BKG broadcast navigation
    -> PyTECGg GNSS context + ephemerides
    -> geometry-free phase/code + Melbourne-Wubbena combinations
    -> arc extraction / levelling
    -> satellite coordinates
    -> IPP / azimuth / elevation, 30 deg mask
    -> TEC calibration
    -> calibrated sTEC / vTEC
    -> station VEq
    -> CSV / Parquet / daily PNG products
```

## Expected directory layout

```text
DATA_ROOT/
└── 2025/
    └── 01/
        └── 01/
            ├── RINEX/
            │   └── *_30S_*.rnx
            ├── PyTECGg_INPUT/
            │   └── NAV/
            └── PyTECGg_OUTPUT/
```

If no navigation file is already present, PyTECGg downloads the BKG BRDC
navigation product.

## Usage

Whole year:

```bat
python process_koh2_pytecgg.py ^
  --year 2025 ^
  --data-root "D:\path\to\data"
```

Single validation date:

```bat
python process_koh2_pytecgg.py ^
  --year 2025 ^
  --date 2025-01-01 ^
  --data-root "D:\path\to\data" ^
  --force
```

`--data-root` may point either to the parent of the year folders or directly
to the requested year directory.

## Outputs

For each processed day the script writes products such as:

```text
KOH2_001_2025_PyTECGg_GEOMETRY.parquet
KOH2_001_2025_PyTECGg_TEC.parquet
KOH2_001_2025_PyTECGg_TEC.csv
KOH2_001_2025_PyTECGg_VEQ.parquet
KOH2_001_2025_PyTECGg_VEQ.csv
KOH2_001_2025_PyTECGg_MANIFEST.txt
```

Daily PNG plots are also generated.

## Validation status

The publication-oriented wrapper was validated on KOH2 data for 2025-01-01
(DOY 001) using the same prepared 30-second RINEX observation file and the same
existing BKG broadcast navigation product as the operational workflow.

The resulting TEC and VEq Parquet files were bit-for-bit identical to the
corresponding operational outputs. Exact SHA-256 checksums and test scope are
documented in `VALIDATION.md`.

This is a representative spot check and does not imply that every day in the
full archive was independently checksum-validated.
