# KOH2 validation against IGS Final GIM

`validate_tec_igs_2019_2026.py` is the publication-oriented wrapper for the
validated KOH2 IGS Final GIM comparison workflow.

The scientific comparison logic is retained from the operational V4 script:

- PyTECGg VEq vs IGS Final GIM at the KOH2 station position;
- PyTECGg VTEC vs IGS Final GIM at the PyTECGg IPP;
- pyOASIS VTEC vs IGS Final GIM at the pyOASIS IPP;
- 30 degree elevation cutoff for IPP comparisons;
- bilinear latitude/longitude interpolation of IONEX;
- linear interpolation in time;
- residual convention: GNSS-derived VTEC minus IGS GIM VTEC.

The publication refactor changes runtime configuration only: data/output paths,
year selection and an optional representative-date filter are supplied by CLI.

## Full study interval

```bat
python validate_tec_igs_2019_2026.py ^
  --data-root "E:\KOH2data"
```

## One representative validation date

```bat
python validate_tec_igs_2019_2026.py ^
  --data-root "E:\KOH2data" ^
  --output-root "D:\path\to\TEST_IGS_GITHUB" ^
  --ionex-cache-root "E:\KOH2data\TEC_VALIDATION_IGS_2019_2026\_IGS_IONEX" ^
  --year 2025 ^
  --date 2025-01-01
```

`--output-root` should be set to a separate test location during wrapper
validation so the existing operational validation products are not overwritten.

Earthdata credentials are still read from the normal Windows `_netrc`/`.netrc`
configuration when an IONEX file must be downloaded.


## Validation status

The publication wrapper has been validated on the representative KOH2 day
2025-01-01 (DOY 001). All three daily matched Parquet comparison products were
bit-for-bit identical to the corresponding operational V4 outputs.

See `VALIDATION.md` for exact SHA-256 checksums and test scope.
