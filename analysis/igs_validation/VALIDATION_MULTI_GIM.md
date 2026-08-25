# Multi-GIM validation

## Representative wrapper-equivalence test

The publication-oriented `validate_tec_multi_gim_2019_2026.py` wrapper was
validated against the existing operational multi-GIM workflow using:

```text
Station: KOH2
Date: 2025-01-01
DOY: 001
Reference GIMs: CODE, ESA, JPL, UPC
```

The test used the same production PyTECGg and pyOASIS outputs and the same
cached final IONEX products used by the operational workflow.

The publication wrapper completed all expected comparisons:

```text
3 comparison types × 4 GIM products = 12 daily-statistics rows
```

For every GIM, the matched sample counts were:

```text
PyTECGg VEq at station : 2,881
PyTECGg VTEC at IPP    : 21,007
pyOASIS VTEC at IPP    : 27,150
```

## Exact table equivalence

The publication result:

```text
KOH2_2019_2026_multiGIM_daily_statistics.csv
```

was compared with the operational multi-year table after filtering the
operational table to:

```text
year = 2025
doy  = 001
```

Both tables were sorted by:

```text
reference_gim
comparison
```

The comparison returned:

```text
new rows: 12
old rows: 12
same columns: True
EXACT EQUALITY: True
```

Therefore, all identifying fields and all numerical statistics for the
representative day are exactly equal between the publication wrapper and the
operational workflow.

## Validation conclusion

For KOH2 on 2025-01-01, the publication-oriented multi-GIM wrapper reproduces
the operational daily statistics exactly for CODE, ESA, JPL and UPC and for all
three comparison types.

This is a representative wrapper-equivalence spot check. It does not imply
that every observation day from 2019-2026 was independently revalidated.

CODE, ESA, JPL and UPC final GIMs are products of IGS ionospheric analysis
centers. The comparison should therefore be described as a multi-GIM or
inter-product robustness comparison, not as four statistically independent
validations.
