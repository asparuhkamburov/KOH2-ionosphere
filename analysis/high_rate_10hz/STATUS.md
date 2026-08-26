# High-rate 10 Hz module status

## Status: EXPERIMENTAL

This module is intentionally included in the public repository, but it is not
included in the repository's validated `PASS` claim.

The current primary phase quantity is:

```text
SIGMA_PHI_GF_EQUIV_RAD
```

It is an experimental dual-frequency geometry-free phase-fluctuation proxy.

The current amplitude-side quantity is:

```text
S4_CNO_PROXY
```

It is an uncalibrated C/N0-derived proxy.

Neither quantity should be reported as a reference-grade ISMR scintillation
index.

## Why the status is EXPERIMENTAL rather than PENDING

The implementation is no longer untested. The following checks have been
completed:

- one-hour streaming-parser equivalence against the earlier GeoRinex-based
  implementation;
- two-hour boundary test demonstrating removal of the artificial hourly
  filter-restart artifact;
- empirical 120 s edge-guard sensitivity assessment;
- SP3 geometry and elevation screening;
- RINEX LLI inspection;
- cross-frequency scaling diagnostic demonstrating contamination of the
  original single-frequency phase metric by nondispersive/common-mode content;
- full 24-hour continuous geometry-free processing for 2025-01-01;
- targeted R02 geometry-free implementation cross-check;
- 30/35/40/45 deg elevation sensitivity;
- event-level comparison with pyOASIS ROTI and DTEC.

These checks justify inclusion as an experimental reproducible method, but not
promotion to a validated physical scintillation product.

## Approved publication wording

Preferred:

> Experimental dual-frequency geometry-free phase-fluctuation proxy derived
> from 10 Hz geodetic RINEX observations.

For the representative event:

> A geometry-free phase-fluctuation enhancement associated with a disturbed
> ionospheric interval was observed, with partial satellite-specific support
> from pyOASIS ROTI/DTEC.

Avoid:

- "validated Phi60";
- "standard scintillation index";
- "ISMR-equivalent Phi60";
- "calibrated S4";
- causal statements unless supported by a separate causal analysis.

## Promotion criteria

Promotion beyond EXPERIMENTAL should require at least one strong external
event-level validation route, preferably dedicated scintillation-monitor data
or another independently calibrated high-rate GNSS scintillation product, plus
repeatability across additional disturbed and quiet days.
