# Validation

## Representative wrapper-equivalence test

The publication-oriented `run_pyiri_koh2_2019_2026.py` wrapper was checked
against the existing operational workflow using:

```text
Station : KOH2
Date    : 2025-01-01
DOY     : 001
PyIRI   : 0.1.7
Option  : URSI
UT grid : hourly, 24 samples
TEC integration : 90-2000 km, 5 km step
F10.7   : 211.9 sfu
Source  : NASA OMNI daily F10.7
```

The existing cached solar-index file was reused so that the model input was
identical to the operational run.

Representative result:

```text
mean VTEC   = 25.367487 TECU
median VTEC = 25.556149 TECU
range       = 21.138 .. 29.593 TECU
```

## Bit-for-bit output equivalence

File:

```text
KOH2_2025_001_PyIRI_URSI_hourly.csv
```

SHA-256 for both publication and operational outputs:

```text
02bc39204e04ed7769aeaf449bef93ec020ad6611cc1cc5b0d68b72b75cab1fa
```

Therefore the representative hourly PyIRI output is bit-for-bit identical.

## Interpretation

PyIRI is used here as an empirical climatological ionospheric background driven
by solar activity. It is not treated as an independent validation data source.

This is a representative wrapper-equivalence spot check and does not imply
that every modeled day from 2019-2026 was independently checksum-validated.
