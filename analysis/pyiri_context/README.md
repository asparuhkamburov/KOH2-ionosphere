# KOH2 PyIRI climatological context

`run_pyiri_koh2_2019_2026.py` generates empirical/climatological IRI VTEC at
the fixed KOH2 coordinates for the available GNSS observation days.

Scientific settings retained from the operational workflow:

- PyIRI implementation of IRI climatology;
- URSI foF2 coefficient option;
- NASA OMNI daily observed F10.7 forcing;
- hourly UT epochs, 00:00 through 23:00;
- vertical electron-density integration from 90 to 2000 km;
- 5 km vertical integration step.

## Scientific interpretation

PyIRI is used as an **empirical climatological background/context**. It should
not be described as an independent validation of the KOH2 GNSS-derived TEC.

Daily F10.7 is an explicit model input, so PyIRI must also not be treated as an
independent test of the relationship between TEC and solar activity.

## Representative wrapper-equivalence test

```bat
python run_pyiri_koh2_2019_2026.py ^
  --data-root "F:\KOH2data" ^
  --output-root "D:\path\to\PYIRI_OUTPUT" ^
  --indices-root "<PYIRI_INDICES_ROOT>" ^
  --year 2025 ^
  --date 2025-01-01
```

For validation, reusing the existing operational `_INDICES` directory ensures
that the publication wrapper uses exactly the same cached NASA OMNI F10.7
input file as the operational run.


## Validation status

The publication wrapper passed representative bit-for-bit equivalence
validation for KOH2 on 2025-01-01 (DOY 001).

The 24-hour PyIRI hourly output file was identical to the operational output.
See `VALIDATION.md` for the exact SHA-256 checksum and test scope.
