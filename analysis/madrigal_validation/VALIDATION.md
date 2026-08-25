# Validation

## Representative validation date

Both publication-oriented Madrigal wrappers were checked against their
operational counterparts using:

```text
Station : KOH2
Date    : 2025-01-01
DOY     : 001
```

The existing production PyTECGg, pyOASIS, Madrigal cache products and IGS Final
IONEX cache were reused.

## 1. KOH2 processors vs Madrigal

Script:

```text
validate_tec_madrigal_2019_2026.py
```

The publication wrapper reproduced all three daily comparisons:

```text
PyTECGg_VEq_vs_Madrigal_station
PyTECGg_VTEC_vs_Madrigal_IPP
pyOASIS_VTEC_vs_Madrigal_IPP
```

Representative-day matched samples:

```text
PyTECGg VEq station : 1,130 / 2,881
PyTECGg VTEC IPP    : 8,786 / 21,007
pyOASIS VTEC IPP    : 11,081 / 27,150
```

Exact comparison with the operational daily-statistics table after filtering
to year=2025 and doy=001:

```text
new rows: 3
old rows: 3
same columns: True
EXACT EQUALITY: True
```

The Madrigal daily receiver-site check reported:

```text
KOH2 exact input: False
nearest input receiver: sgp3
nearest input distance: 398.6 km
```

## 2. Direct IGS Final GIM vs Madrigal

Script:

```text
compare_igs_madrigal_direct_2019_2026.py
```

Residual convention:

```text
IGS VTEC - Madrigal VTEC
```

Representative-day result:

```text
n    = 13,053
bias = +4.634950 TECU
MAE  = 5.705847 TECU
RMSE = 6.753109 TECU
r    = 0.801072
```

Exact comparison with the operational daily-statistics table after filtering
to year=2025 and doy=001:

```text
new rows: 1
old rows: 1
same columns: True
EXACT EQUALITY: True
```

## Validation conclusion

For KOH2 on 2025-01-01, both publication wrappers reproduce their operational
daily statistics exactly.

These are representative wrapper-equivalence spot checks. They do not imply
that every available observation day from 2019-2026 was independently
checksum-validated.
