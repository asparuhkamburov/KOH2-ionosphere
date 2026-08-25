# Validation

## Representative test

```text
Station : KOH2
Date    : 2025-01-01
DOY     : 001
```

Input availability on the representative day:

```text
PyIRI hourly epochs : 24
PyTECGg matched      : 24
IGS matched          : 24
Madrigal matched     : 11
All four available   : 11
```

The publication wrapper was compared directly with the operational workflow.

## Common-hour values table

```text
new rows: 24
old rows: 24
same columns: True
EXACT EQUALITY: True
```

## Daily pairwise statistics table

```text
new rows: 6
old rows: 6
same columns: True
EXACT EQUALITY: True
```

The six daily pairwise comparisons on the representative day were:

```text
IGS_minus_Madrigal
IGS_minus_PyIRI
Madrigal_minus_PyIRI
PyTECGg_VEq_minus_IGS
PyTECGg_VEq_minus_Madrigal
PyTECGg_VEq_minus_PyIRI
```

## Validation conclusion

For KOH2 on 2025-01-01, the publication-oriented common-hour wrapper reproduces
the operational hourly values and pairwise daily statistics exactly.

This is a representative wrapper-equivalence spot check. It does not imply
that every available day from 2019-2026 was independently revalidated.
