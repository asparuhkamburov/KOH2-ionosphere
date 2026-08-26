# Validation record — high-rate 10 Hz module

## Scope

This document records software-equivalence, QC, sensitivity, and event-level
cross-checks for the experimental KOH2 10 Hz module.

It does **not** claim validation of a standard ISMR `Phi60` or calibrated `S4`
product.

Representative date: **2025-01-01 (DOY 001)**.

## 1. Streaming parser equivalence

One 10 Hz hourly RINEX file was processed with the publication-oriented
streaming parser and compared with the earlier GeoRinex-based reference table.

```text
reference rows                         2,345
streaming rows                         2,345
(window_start, sv, signal) keys        exact
max |delta sigma_phi_rad|              1.2689895e-08 rad
max |delta S4_CNO_PROXY|               1.4297591e-14
```

Status: **PASS for the tested one-hour implementation-equivalence case**.

This does not validate the physical interpretation of either metric.

## 2. Hour-boundary artifact test

The old workflow filtered hourly files independently. At the 01:00 UTC
boundary, the old and continuous implementations gave:

```text
old hourly sigma_phi maximum            325.427228869 rad
continuous sigma_phi maximum              1.326357826 rad
old hourly median                       180.632564355 rad
continuous median                         0.127061361 rad
continuous QC edge flags at 01:00         0 / 43
```

The reduction demonstrates that the hundreds-of-radians 01:00 spike was an
hourly filter-restart artifact rather than a robust geophysical feature.

Status: **PASS as an internal boundary/QC test**.

## 3. Arc-edge sensitivity

On the representative full day, the strongest single-frequency phase values
were concentrated close to true filtered-arc edges.

After excluding the first/last 120 s of each phase segment:

```text
median sigma_phi_rad                     0.109704 rad
p95                                      0.161772 rad
p99                                      0.201967 rad
maximum                                  0.569141 rad
```

The 120 s guard was therefore retained for the experimental workflow.

## 4. Cross-frequency diagnostic

After continuous processing, SP3 geometry, elevation >=30 deg, and edge QC,
the observed phase-sigma ratios were much closer to nondispersive frequency
scaling than to the inverse-frequency scaling expected for first-order
ionospheric carrier phase.

Examples:

```text
GPS L5X/L1C observed median ratio         0.7463
nondispersive f5/f1                       0.7468
first-order ionospheric f1/f5             1.3391

R02 L2C/L1C observed median ratio         0.7868
nondispersive f2/f1                       0.7778
first-order ionospheric f1/f2             1.2857
```

Interpretation: the original single-frequency `SIGMA_PHI_RAD` contains
substantial nondispersive/common-mode content and should not be promoted to a
physical `Phi60` result.

This diagnostic motivated the geometry-free implementation.

## 5. R02 LLI / geometry-free targeted diagnostic

For R02, 14:05–14:25 UTC:

```text
phase samples                            24,002
LLI bit-0 samples                            0
LLI bit-1 samples                            0
LLI bit-2 samples                            0
maximum |delta geometry-free|          0.008647 m
robust GF jump threshold              0.024767 m
raw GF jump flags                            0
```

The old single-frequency R02/L1C peak near 14:16 UTC did not survive the
geometry-free transformation:

```text
old single-frequency sigma             about 0.3165 rad
GF equivalent sigma at 14:16             0.079334 rad
GF interval median                       0.085421 rad
```

Thus the old peak was not a robust ionospheric phase-scintillation candidate.

## 6. Targeted vs full-day GF implementation cross-check

For R02 at 2025-01-01 14:16 UTC:

```text
targeted diagnostic                      0.0793336498826335 rad
full-day continuous processor            0.0793336607155998 rad
absolute difference                      1.0833e-08 rad
```

Status: **PASS for representative implementation consistency**.

## 7. Full 24-hour continuous GF run

Representative 2025-01-01 result:

```text
files processed                                 24
epoch records                              864,000
satellite records                       15,306,248
paired geometry-free samples           14,257,040
satellite GF series                             52
output minute rows                           23,867
QC-valid GF rows                             12,605
QC-valid GPS rows                             7,342
QC-valid GLONASS rows                         5,263
```

Primary elevation >=30 deg distribution:

```text
median                                    0.066473 rad
p95                                       0.096742 rad
p99                                       0.118494 rad
maximum                                   0.191216 rad
```

The maximum occurs at approximately 30.37 deg elevation and is therefore
sensitive to the elevation mask.

## 8. Elevation sensitivity

```text
mask     N       median      p95       p99       maximum
>=30   12605     0.06647    0.09674   0.11849   0.19122
>=35   10493     0.06457    0.08997   0.10692   0.18978
>=40    8713     0.06281    0.08561   0.10069   0.12090
>=45    6929     0.06097    0.08480   0.09954   0.11019
```

The residual GF amplitude is elevation-dependent. The publication workflow
therefore retains 30 deg as the primary mask and uses 35/40 deg as explicit
sensitivity cases rather than silently replacing the primary threshold.

## 9. Event-level upper-tail cluster

At elevation >=40 deg, an upper-tail multi-link cluster occurs around
16:50–17:10 UTC.

Using a per-satellite daily p95 threshold:

```text
16:30–16:50 UTC   minutes with >=2 high links    0 / 20
16:50–17:10 UTC   minutes with >=2 high links    8 / 20
                  minutes with >=3 high links    2 / 20
17:10–17:30 UTC   minutes with >=2 high links    0 / 20
```

The event-window median remains close to the pre/post medians, so the effect is
an upper-tail/subset-link enhancement rather than a bulk shift.

Representative simultaneous links include G01/R21 near 16:56 UTC, G01/G10 near
16:57 UTC, and R21/R04 near 17:04 UTC.

## 10. pyOASIS ROTI/DTEC event cross-check

The comparison is **partial and satellite-specific**.

R04 provides the strongest supporting case at elevation >=40 deg:

```text
R04 event ROTI maximum                       0.86942
R04 daily median ROTI                        0.06023
R04 daily p95 ROTI                           0.38620

R04 event |DTEC| maximum                     4.93158 TECU
R04 daily median |DTEC|                      1.69705 TECU
R04 daily p95 |DTEC|                         3.91065 TECU
```

The R04 GF enhancement near 17:03–17:05 overlaps the elevated ROTI/DTEC
interval.

R21 does not show the same ROTI/DTEC confirmation during the strongest GF
interval, and the available pyOASIS GPS products do not contain the strongest
GF GPS candidates G01 and G10.

Therefore the pyOASIS comparison supports a disturbed ionospheric interval but
does not independently validate the full multi-satellite GF cluster as
scintillation.

## 11. Full-day matched association with pyOASIS

Nearest-time matched comparisons at elevation >=40 deg show only modest
association overall and stronger correspondence for GLONASS than GPS.

Within-satellite rank correlations:

```text
               ALL       GPS       GLONASS
GF vs ROTI     0.197     0.004      0.250
GF vs |DTEC|   0.283    -0.036      0.371
```

These results reinforce the need for constellation-specific caution.

## Final status

**EXPERIMENTAL — reproducible and internally tested, but not externally
validated as a standard scintillation product.**

The recommended interpretation is:

> Experimental dual-frequency geometry-free phase-fluctuation proxy, with
> partial event-level support from independent pyOASIS ionospheric disturbance
> indicators.

The module remains outside the repository's `PASS` validation claims.
