# PyIRI common-hour comparison

This workflow places four station-level TEC series on the exact hourly PyIRI
epochs:

1. PyIRI climatological VTEC;
2. PyTECGg VEq at KOH2;
3. IGS Final GIM interpolated at KOH2;
4. nearest available Madrigal VTEC cell around KOH2.

## Scientific interpretation

PyIRI is used as an empirical climatological background. Because F10.7 is an
input to PyIRI, this comparison must not be described as an independent
solar-activity validation.

## Matching

- PyTECGg: nearest VEq epoch within 60 s of the hourly PyIRI epoch.
- IGS Final GIM: interpolated to the exact hourly PyIRI epoch at KOH2.
- Madrigal: nearest epoch within 180 s and nearest available spatial cell
  within 80 km.

## Validation status

Representative numerical equivalence was confirmed for KOH2 on 2025-01-01
(DOY 001).

The publication wrapper reproduced both operational output tables exactly:

```text
common_hour_values.csv
  new rows: 24
  old rows: 24
  same columns: True
  EXACT EQUALITY: True

common_hour_daily_statistics.csv
  new rows: 6
  old rows: 6
  same columns: True
  EXACT EQUALITY: True
```

See `VALIDATION.md` for the detailed validation record.
