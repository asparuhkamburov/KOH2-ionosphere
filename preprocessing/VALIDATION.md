# Validation record

## KOH2 2025-01-01 (DOY 001)

The publication-oriented GFZRNX preprocessing wrapper was tested on the
24 hourly KOH2 RINEX files for 2025-01-01.

Input sampling token:

```text
10Z
```

The tested processing sequence was:

```text
24 hourly RINEX files
    -> GFZRNX concatenation + -smp 1
    -> KOH200ATA_R_20250010000_01D_01S_MO.rnx
    -> GFZRNX -smp 30
    -> KOH200ATA_R_20250010000_01D_30S_MO.rnx
```

The whole-file SHA-256 values differed from the historical operational files
because the RINEX headers contain processing-time metadata.

To isolate the GNSS observation content, SHA-256 was therefore computed only
for the RINEX body after `END OF HEADER`.

### Daily 1 s product

Publication wrapper:

```text
d6f56ca15a4507e290d76f43c722be3dea179c9a12e78f929e2e84ef200ad9ef
```

Historical operational product:

```text
d6f56ca15a4507e290d76f43c722be3dea179c9a12e78f929e2e84ef200ad9ef
```

Result: **bit-for-bit identical observation body**.

### Daily 30 s product

Publication wrapper:

```text
9b1cb59276bb980be2d5b16c33e528ed82050e561fc08332592af7cb383be8b4
```

Historical operational product:

```text
9b1cb59276bb980be2d5b16c33e528ed82050e561fc08332592af7cb383be8b4
```

Result: **bit-for-bit identical observation body**.

## Scope of this validation

This test validates the hourly-to-daily concatenation and 10 Hz -> 1 Hz ->
30 s sampling path for the tested 2025 day.

It does **not** by itself validate:

- the Trimble T02 -> RINEX conversion wrapper;
- the legacy `*.YYo` header-harmonization/RINEX-3 conversion stage;
- all years or all possible input layouts.

Additional representative tests can be recorded before the final monograph
release.


## KOH2 2019-12-28 (DOY 362)

A second validation test was performed on a representative **1 Hz** KOH2 day
from the 2019-2024 observation period.

Input:

```text
24 hourly RINEX files
sampling token: 01S
```

Processing sequence:

```text
24 hourly 1 s RINEX files
    -> GFZRNX concatenation + -smp 1
    -> KOH200ATA_R_20193620000_01D_01S_MO.rnx
    -> GFZRNX -smp 30
    -> KOH200ATA_R_20193620000_01D_30S_MO.rnx
```

As for the 2025 validation, SHA-256 was computed only for the RINEX observation
body after `END OF HEADER`.

### Daily 1 s product

Publication wrapper and historical operational product:

```text
8a4eb9b8c4fc999abf163205d547bd2282e62aadd55fa73aa5513288752b9510
```

Result: **bit-for-bit identical observation body**.

### Daily 30 s product

Publication wrapper and historical operational product:

```text
96db4492274d8ceed36d525ba7bda851dd5e92e645567e2e7c64949a6e4d789a
```

Result: **bit-for-bit identical observation body**.

## Current validation coverage

The publication-oriented GFZRNX wrapper has now reproduced the historical
operational observation bodies for both principal KOH2 sampling regimes:

- 2019-2024: representative 1 Hz day validated (2019-12-28);
- 2025-2026: representative 10 Hz day validated through 1 s and 30 s products
  (2025-01-01).

The remaining preprocessing path that has not yet been independently spot-
checked is the legacy `*.YYo` -> CRUX header harmonization / RINEX 3
standardization stage.
