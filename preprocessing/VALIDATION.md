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
