# KOH2 Solar Cycle 25 monthly VTEC + F10.7

This directory contains the publication-oriented script:

```text
plot_koh2_monthly_real_vtec_f107.py
```

It uses the validated common-hour table produced by the PyIRI/common-hour
workflow.

## Recommended scientific product

The strict-common-hour result is the recommended direct inter-product figure.

The aggregation is:

```text
strict common hourly epochs
-> one daily mean per VTEC product
-> equal-day monthly mean
```

Therefore each available observation day receives equal weight in the monthly
mean.

F10.7 is first reduced to one value per day and then averaged monthly.

Missing months are retained as gaps and are not synthetically interpolated.
Coverage is reported separately through the number of strict common hours and
common observation days per month.

## Compared series

- PyTECGg VEq
- IGS Final GIM VTEC at KOH2
- Madrigal GNSS VTEC near KOH2
- PyIRI climatological VTEC
- NASA OMNI daily F10.7

PyIRI is used as empirical climatological context rather than as an independent
validation product.

## Main outputs

```text
KOH2_daily_strict_common_hour_means.csv
KOH2_monthly_strict_common_hour_statistics.csv
KOH2_monthly_mean_VTEC_F107_strict_common_hours.png
KOH2_monthly_mean_VTEC_F107_strict_common_hours.pdf
KOH2_monthly_mean_VTEC_F107_strict_common_hours.svg
```

An available-data alternative is also produced:

```text
KOH2_monthly_available_data_statistics.csv
KOH2_monthly_mean_VTEC_F107_available_data.png
```

This alternative is descriptive only because the products may use different
valid epochs. Use the strict-common-hour products for direct inter-product
comparison.

## Usage

```bat
python plot_koh2_monthly_real_vtec_f107.py ^
 --input-file "path\to\KOH2_2019_2026_common_hour_values.csv" ^
 --output-dir "path\to\MONTHLY_FIGURES"
```

If `--output-dir` is omitted, `MONTHLY_FIGURES` is created beside the input
CSV. The default date window is 2019-01-01 through 2026-12-31.

Optional switches:

```text
--start-date YYYY-MM-DD
--end-date YYYY-MM-DD
--no-sem-bands
--no-high-activity-shading
```

## Validation

The publication-oriented script passed syntax validation and complete supplied
dataset numerical equivalence testing against the operational implementation.
See `VALIDATION.md` for the recorded test and results.
