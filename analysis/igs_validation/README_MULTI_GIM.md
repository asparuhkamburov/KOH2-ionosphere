# KOH2 multi-GIM comparison

`validate_tec_multi_gim_2019_2026.py` extends the validated IGS Final combined
GIM workflow to individual final GIM products from:

- CODE
- ESA
- JPL
- UPC

This is a **multi-GIM / inter-product comparison**, not four statistically
independent validations.

The script deliberately imports the publication core:

```text
validate_tec_igs_2019_2026.py
```

from the same directory and reuses its KOH2 coordinates, PyTECGg reader,
pyOASIS reader, IONEX parser, spatial/temporal interpolation, elevation cutoff,
residual definition and statistics.

## Representative test

```bat
python validate_tec_multi_gim_2019_2026.py ^
  --data-root "F:\KOH2data" ^
  --output-root "D:\path\to\TEST_MULTI_GIM\OUTPUT" ^
  --gim-cache-root "<MULTIGIM_IONEX_CACHE_ROOT>" ^
  --year 2025 ^
  --date 2025-01-01
```

`validate_tec_igs_2019_2026.py` must be in the same directory.

The scientific defaults are unchanged from the operational script. Runtime
paths and optional representative-date selection are supplied by CLI.


## Validation status

Representative-day wrapper equivalence was confirmed for KOH2 on 2025-01-01.
The publication and operational workflows produced exactly the same 12
daily-statistics rows across CODE, ESA, JPL and UPC after matching the same
three comparison types.

See `VALIDATION_MULTI_GIM.md` for the test scope and exact comparison result.
