# PyTECGg validation

## Representative validation test

The publication-oriented `process_koh2_pytecgg.py` wrapper was validated
against the operational PyTECGg workflow using:

```text
Station : KOH2
Date    : 2025-01-01
DOY     : 001
```

The test used the same prepared daily 30-second RINEX observation file:

```text
KOH200ATA_R_20250010000_01D_30S_MO.rnx
```

and the same existing BKG broadcast navigation file:

```text
BRDC00IGS_R_20250010000_01D_MN.rnx.gz
```

The wrapper processed one daily 30-second RINEX file and reported:

```text
Parsed observation records : 365,757
Duplicate records removed  : 0
Records entering PyTECGg   : 365,757
GNSS systems               : G, R
Prepared ephemerides       : 57
Linear-combination rows    : 51,038
Arc-levelled rows          : 51,038
Rows after elevation mask  : 27,622
TEC rows                   : 27,622
VEq rows                   : 27,622
```

A PyTECGg warning was issued for missing GLONASS ephemeris information for
R09 for one or more epochs. The processing completed successfully.

## TEC output equivalence

Publication wrapper output and operational reference:

```text
KOH2_001_2025_PyTECGg_TEC.parquet
```

produced the same SHA-256 checksum:

```text
1a359836ee3321940224eae4b4ac53b4ea1695edfa353e1c2c54005f12cd9845
```

Result: **bit-for-bit identical TEC Parquet product**.

## VEq output equivalence

Publication wrapper output and operational reference:

```text
KOH2_001_2025_PyTECGg_VEQ.parquet
```

produced the same SHA-256 checksum:

```text
78ad2464de6cc919b150f475b3d64792310795373a3b4bbccc7f72487304dd54
```

Result: **bit-for-bit identical VEq Parquet product**.

## Validation conclusion

For the representative KOH2 day 2025-01-01, the refactored
publication-oriented wrapper reproduces the operational PyTECGg TEC and VEq
Parquet outputs bit-for-bit while preserving the original scientific
processing parameters and sequence.

This validation is a representative spot check. It does not mean that every
day in the complete KOH2 archive was independently checksum-validated.
