# Validation status

## Result

Syntax validation: **PASS**.

Complete supplied-dataset numerical equivalence: **PASS**.

The publication-oriented script was tested against the operational implementation
using the complete supplied 2019-2026 common-hour CSV rather than a one-day spot
check. The processing run used 7,824 input rows.

## Recorded strict-common-hour coverage

```text
Common observation days: 325
Strict common hours:      6,983
Months with data:         19
```

## Numerical comparison

The publication and operational versions were compared using exact pandas
`DataFrame.equals()` equality for the following CSV products:

```text
KOH2_monthly_strict_common_hour_statistics.csv
  publication rows: 96
  operational rows: 96
  same columns: True
  EXACT EQUALITY: True

KOH2_daily_strict_common_hour_means.csv
  publication rows: 325
  operational rows: 325
  same columns: True
  EXACT EQUALITY: True

KOH2_monthly_available_data_statistics.csv
  publication rows: 96
  operational rows: 96
  same columns: True
  EXACT EQUALITY: True
```

No numerical differences were detected in the three compared CSV products.
Accordingly, the publication-oriented implementation reproduces the operational
CSV-level numerical products exactly for the complete supplied common-hour
dataset.

## Validation scope

This result demonstrates implementation equivalence for the supplied harmonized
common-hour dataset. It should not be interpreted as an independent scientific
validation of the underlying TEC products.

Figure files are not used for bit-for-bit validation because rendering metadata
can differ across Matplotlib versions and platforms even when the numerical data
are identical.

## Syntax check

The publication script is checked with:

```bat
python -m py_compile plot_koh2_monthly_real_vtec_f107.py
```

The recorded syntax check completed successfully.

## Runtime note

Pandas may emit warnings when timezone-aware timestamps are converted to monthly
periods. These warnings did not stop execution and did not affect the validated
CSV outputs.
