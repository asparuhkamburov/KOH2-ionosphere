# KOH2 Madrigal validation

This directory contains two complementary Madrigal analyses.

## 1. `validate_tec_madrigal_2019_2026.py`

Primary observation-based comparison of:

- PyTECGg VEq at KOH2;
- PyTECGg VTEC at PyTECGg IPPs;
- pyOASIS VTEC at pyOASIS IPPs;

against MIT Haystack / CEDAR Madrigal MAPGPS gridded GNSS VTEC.

The scientific matching method is unchanged from the operational V4 FAST
workflow:

- nearest Madrigal 5-minute epoch;
- absolute time difference <= 180 s;
- nearest available Madrigal grid cell;
- horizontal separation <= 80 km;
- elevation cutoff 30 degrees for IPP comparisons;
- residual = KOH2 method TEC - Madrigal VTEC.

The daily Madrigal receiver-site list is also inspected when available. If
KOH2 contributed to a day's MAPGPS product, that comparison must not be
described as fully independent validation.

## 2. `compare_igs_madrigal_direct_2019_2026.py`

Diagnostic reference-product comparison. It removes PyTECGg and pyOASIS and
compares IGS Final GIM directly with cached Madrigal VTEC at the exact Madrigal
bin epochs and locations inside a fixed KOH2 regional box.

Residual convention:

```text
IGS VTEC - Madrigal VTEC
```

## Dependency

Both scripts import:

```text
validate_tec_igs_2019_2026.py
```

from the same directory.

## Validation status

Both publication wrappers passed representative numerical-equivalence
validation for KOH2 on 2025-01-01 (DOY 001).

- KOH2 vs Madrigal: 3/3 daily-statistics rows exactly equal to the operational
  V4 FAST workflow.
- Direct IGS vs Madrigal: 1/1 daily-statistics row exactly equal to the
  operational direct-reference workflow.

See `VALIDATION.md` for details.
