# KOH2 Solar Cycle 25 analysis

This directory contains the publication-oriented Solar Cycle 25 / geomagnetic-analysis scripts for KOH2 (2019-2026).

## Main workflow

1. `analyze_koh2_solar_cycle25_2019_2026.py`
   - builds the daily solar/geomagnetic master table;
   - merges validated KOH2 common-hour VTEC, optional IGS/Madrigal validation statistics, OMNI2, SYM-H and SILSO data;
   - calculates F10.7 + seasonal background models, correlations, background residuals, storm/quiet summaries and equal-day yearly summaries.
2. `analyze_koh2_solar_cycle25_sensitivity.py`
   - compares ordinary daily OLS with equal-year weighted WLS;
   - performs pooled and within-year storm/quiet sensitivity analyses with multiple-testing control.
3. `analyze_koh2_solar_cycle25_lagged_bootstrap.py`
   - evaluates lag 0/+1/+2 calendar-day geomagnetic associations;
   - uses year-cluster and 7-day calendar moving-block bootstrap confidence intervals.
4. `plot_koh2_monthly_real_vtec_f107.py`
   - creates the validated strict-common-hour monthly VTEC + F10.7 figure and supporting CSV tables.
5. `figures/`
   - contains four portable publication figure suites for validation benchmarks, solar/reference divergence, geomagnetic response, and multi-GIM/coverage summaries.

`compare_csv_outputs.py` is a validation helper for exact old-vs-new CSV comparisons.

## Scientific interpretation

- PyIRI is an empirical climatological background, not an independent solar-activity validation, because F10.7 is an input to the model.
- CODE/ESA/JPL/UPC comparisons are multi-GIM/inter-product robustness checks, not four independent validations.
- Solar/geomagnetic results should be described as *associated with*, *related to*, or *corresponding to* activity unless a separate causal analysis is performed.
- No empirical TEC bias correction is applied.
- Calendar-day lags are true calendar lags, not previous available KOH2 observation rows.

## Reproducibility

The publication scripts do not contain machine-specific data paths. Inputs and output locations are supplied using CLI arguments. For strict equivalence testing of the main analysis, reuse the same existing `_INDEX_CACHE` as the operational run so the input index files are identical.

See `CMD_TEST_COMMANDS.txt` for the Windows CMD test sequence, including the downstream figure-suite checks.
