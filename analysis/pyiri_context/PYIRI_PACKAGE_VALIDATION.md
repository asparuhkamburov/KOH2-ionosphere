# PyIRI context package validation

This package contains two publication-oriented PyIRI workflows:

1. `run_pyiri_koh2_2019_2026.py`
   - representative hourly output is bit-for-bit identical to the operational
     output for KOH2 on 2025-01-01;
2. `common_hour/compare_pyiri_common_hour_koh2_2019_2026.py`
   - representative hourly values (24/24 rows) and daily pairwise statistics
     (6/6 rows) are exactly equal to the operational workflow.

PyIRI is used as empirical climatological context, not as an independent
validation source.
