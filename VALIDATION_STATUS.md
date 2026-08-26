# Consolidated validation status

This file summarizes implementation-equivalence checks for the publication-oriented KOH2 processing and analysis repository. Scientific validation and software-output equivalence are distinct concepts; the scope of each result is stated below.

## Preprocessing

Status: **PASS**.

The Trimble/RINEX/GFZRNX preprocessing workflow was checked on representative 2019 and 2025 cases, including T02 conversion and GFZRNX processing equivalence. See `preprocessing/VALIDATION.md` for details.

## pyOASIS

Status: **PASS — representative day**.

Representative date: 2025-01-01 (DOY 001).

Final TEC product:

```text
KOH2_001_2025_L1L2.TEC
```

Operational and publication outputs were bit-for-bit identical. SHA-256:

```text
4280d8a6451adc158cb0b37f350f5b67ef45482587d57e6f6b7096dcfac5634f
```

This establishes implementation equivalence for the tested day, not a checksum of every archive day.

## PyTECGg

Status: **PASS — representative day**.

Representative date: 2025-01-01.

Processing counts recorded for the validation run:

```text
observations                 365,757
ephemerides                       57
LC/arc rows                    51,038
rows after elevation cutoff   27,622
TEC/VEq rows                  27,622
```

Final checksums:

```text
TEC: 1a359836ee3321940224eae4b4ac53b4ea1695edfa353e1c2c54005f12cd9845
VEQ: 78ad2464de6cc919b150f475b3d64792310795373a3b4bbccc7f72487304dd54
```

The missing GLONASS R09 ephemeris warning observed in the representative run was nonfatal.

## IGS Final GIM validation

Status: **PASS — representative day**.

Representative date: 2025-01-01.

Matched values:

```text
PyTECGg VEq at station       2,881
PyTECGg VTEC at IPPs        21,007
pyOASIS VTEC at IPPs        27,150
```

The three representative `DAILY_MATCHED` Parquet products were bit-for-bit identical between operational and publication implementations.

Method scope includes elevation >= 30 degrees, bilinear IONEX spatial interpolation, linear temporal interpolation, and residual convention `method - IGS`.

## Multi-GIM robustness

Status: **PASS — representative day**.

References: CODE, ESA, JPL, UPC.

Representative 2025-01-01 comparison produced the expected 12 rows with exact equality between operational and publication outputs.

This is an inter-product/multi-GIM robustness test, not four statistically independent validations.

## Madrigal comparison

Status: **PASS — representative day**.

Representative date: 2025-01-01.

Matched values:

```text
PyTECGg VEq       1,130
PyTECGg VTEC      8,786
pyOASIS          11,081
```

Operational and publication daily-statistics tables were exactly equal.

For that date, KOH2 was not an exact Madrigal input receiver in the checked site list; the nearest listed receiver was `sgp3` at approximately 398.6 km.

Direct IGS–Madrigal diagnostic on the same date:

```text
n      13,053
bias   +4.634950 TECU
MAE     5.705847 TECU
RMSE    6.753109 TECU
r       0.801072
```

The direct diagnostic output was exactly equal between implementations.

## PyIRI climatological context

Status: **PASS — representative day**.

Configuration includes PyIRI 0.1.7, URSI coefficients, hourly UT sampling, and vertical integration from 90 to 2000 km in 5 km steps.

Representative 2025-01-01:

```text
F10.7        211.9 sfu
mean VTEC    25.367487 TECU
median VTEC  25.556149 TECU
```

Hourly CSV SHA-256:

```text
02bc39204e04ed7769aeaf449bef93ec020ad6611cc1cc5b0d68b72b75cab1fa
```

PyIRI is used as empirical climatological context; F10.7 is model forcing, so it is not an independent solar-activity validation.

## PyIRI common-hour layer

Status: **PASS — representative day**.

Representative 2025-01-01 availability:

```text
PyIRI       24
PyTECGg     24
IGS         24
Madrigal    11
ALL4        11
```

Matching rules include PyTECGg <= 60 s, Madrigal <= 180 s and <= 80 km, and IGS interpolation at the exact hourly KOH2 epochs.

The 24-row values table and 6-row daily-statistics table were exactly equal between operational and publication implementations.

## Solar Cycle 25 monthly strict-common-hour product

Status: **PASS — complete supplied common-hour dataset**.

Validation run:

```text
input rows               7,824
common observation days    325
strict common hours       6,983
months with strict data      19
```

Exact pandas DataFrame equality was obtained for:

```text
KOH2_monthly_strict_common_hour_statistics.csv   96 rows
KOH2_daily_strict_common_hour_means.csv         325 rows
KOH2_monthly_available_data_statistics.csv       96 rows
```

## Solar Cycle 25 main analysis

Status: **PASS — complete supplied analysis dataset**.

The publication implementation reproduced all nine compared operational CSV products with exact pandas DataFrame equality:

```text
KOH2_2019_2026_solar_geomagnetic_daily_master.csv                 326 rows
KOH2_2019_2026_yearly_equal_day_summary.csv                         8 rows
KOH2_2019_2026_correlation_all_days.csv                             60 rows
KOH2_2019_2026_correlation_quiet_days.csv                           20 rows
KOH2_2019_2026_background_regression_F107_season.csv                10 rows
KOH2_2019_2026_background_regression_coefficients.csv               60 rows
KOH2_2019_2026_background_residual_geomagnetic_correlations.csv     40 rows
KOH2_2019_2026_validation_error_geomagnetic_correlations.csv        72 rows
KOH2_2019_2026_storm_vs_quiet_summary.csv                           66 rows
```

Overall exact equality: **TRUE**.

The 326-day daily master and the 325-day strict four-product subset are not contradictory: one daily record lacks a valid PyTECGg value and is excluded from PyTECGg-dependent analyses.

## Solar Cycle 25 sensitivity analysis

Status: **PASS**.

The equal-year WLS / OLS sensitivity products and storm-versus-quiet pooled, within-year, and equal-year summary CSVs were confirmed exactly equal to their operational counterparts.

## Lagged geomagnetic bootstrap analysis

Status: **PASS**.

Default settings retained from the operational workflow:

```text
lags                    0, +1, +2 calendar days
moving-block length     7 calendar days
bootstrap replicates    1000
random seed             20260822
```

A result is marked `bootstrap_robust=True` only when both the year-cluster and 7-day moving-block 95% confidence intervals exclude zero with the same sign.

All compared lagged-bootstrap CSV outputs were confirmed exactly equal to operational results.

## Solar Cycle figure suite 01–04

Status: **PASS for numerical figure-support tables**.

The publication-oriented figure scripts were syntax checked and their generated CSV products were compared with operational outputs. The user validation run reported all comparisons as exact equality `TRUE`.

Rendered PNG/PDF/SVG files are intentionally not used for bit-for-bit validation because rendering metadata, fonts, or Matplotlib versions can change binary files without changing plotted numerical values.

## High-rate 10 Hz analysis

Status: **EXPERIMENTAL — reproducible/internal checks completed; not validated Phi60/S4**.

Representative date: 2025-01-01 (DOY 001).

Completed checks include:

```text
one-hour streaming implementation equivalence
two-hour continuous boundary-artifact test
120 s arc-edge sensitivity
SP3 geometry/elevation screening
RINEX LLI inspection
cross-frequency common-mode diagnostic
targeted R02 geometry-free cross-check
full 24-hour continuous geometry-free run
30/35/40/45 deg elevation sensitivity
event-level pyOASIS ROTI/DTEC comparison
```

One-hour streaming equivalence reproduced 2,345 reference rows with exact
`(window_start, sv, signal)` keys. Maximum numerical differences were
approximately `1.27e-8 rad` for the legacy `sigma_phi_rad` and `1.43e-14` for
`S4_CNO_PROXY`.

The old hourly implementation produced hundreds-of-radians phase spikes at
hour boundaries. Continuous processing removed the 01:00 UTC artifact
(old maximum `325.427 rad`; continuous maximum `1.326 rad`), demonstrating
that the original spike was a filter-restart artifact.

Cross-frequency scaling showed that the remaining single-frequency
`SIGMA_PHI_RAD` was strongly affected by nondispersive/common-mode content.
The primary experimental phase metric was therefore changed to a
dual-frequency geometry-free quantity:

```text
SIGMA_PHI_GF_EQUIV_RAD
```

For the representative 24-hour run:

```text
files processed                           24
paired GF samples                 14,257,040
satellite GF series                       52
output minute rows                    23,867
QC-valid GF rows                      12,605
QC-valid GPS rows                      7,342
QC-valid GLONASS rows                  5,263
median GF sigma                      0.066473 rad
p95                                  0.096742 rad
p99                                  0.118494 rad
```

A targeted R02 diagnostic and the full-day processor agreed to approximately
`1.1e-8 rad` for the 2025-01-01 14:16 UTC representative minute.

Event-level comparison around 16:50–17:10 UTC provides partial
satellite-specific support from pyOASIS ROTI/DTEC, especially for R04, but
does not independently validate the full multi-satellite event as
scintillation.

`S4_CNO_PROXY` remains an uncalibrated C/N0-derived proxy.

The high-rate module is included in the repository for transparency and
reproducibility but is **outside the repository's PASS validation claims**.
See `analysis/high_rate_10hz/README.md`, `STATUS.md`, and `VALIDATION.md`.
