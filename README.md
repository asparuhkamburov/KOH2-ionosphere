# KOH2 GNSS Ionospheric TEC Processing and Validation

Publication-oriented processing and analysis scripts supporting the KOH2 GNSS ionospheric TEC study at Livingston Island, Antarctica.

The repository documents the reproducible software workflow from RINEX preparation through pyOASIS and PyTECGg TEC estimation, external GNSS TEC comparisons, PyIRI climatological context, and Solar Cycle 25 / geomagnetic-association analyses for 2019–2026.

## Repository structure

```text
KOH2-ionosphere/
├── preprocessing/
├── pyoasis/
├── pytecgg/
├── analysis/
│   ├── igs_validation/
│   ├── madrigal_validation/
│   ├── pyiri_context/
│   │   └── common_hour/
│   ├── solar_cycle/
│   │   └── figures/
│   └── high_rate_10hz/
├── README.md
├── VALIDATION_STATUS.md
├── RELEASE_CHECKLIST.md
├── CITATION.cff
├── requirements.txt
└── .gitignore
```

Each validated module contains its own `README.md` and/or `VALIDATION.md` with method-specific details and validation scope.

## Scientific workflow

```text
Trimble / RINEX observations
        |
        v
preprocessing
        |
        +-------------------+
        |                   |
        v                   v
     pyOASIS             PyTECGg
        |                   |
        +---------+---------+
                  |
                  v
         external TEC comparison
         IGS / multi-GIM / Madrigal
                  |
                  v
          PyIRI common-hour context
                  |
                  v
       Solar Cycle 25 / geomagnetic
       analysis, sensitivity, bootstrap
```

## Validation status

The principal publication workflows are validated against their operational implementations. Validation scope differs by module and is stated explicitly; representative-day validation must not be interpreted as whole-archive checksum validation.

| Module | Status | Validation scope |
| --- | --- | --- |
| Preprocessing | PASS | Representative 2019/2025 cases and GFZRNX/T02 equivalence checks |
| pyOASIS | PASS | 2025-01-01 final TEC bit-for-bit equivalence |
| PyTECGg | PASS | 2025-01-01 TEC and VEq bit-for-bit equivalence |
| IGS Final GIM | PASS | Representative 2025-01-01 wrapper/output equivalence |
| Multi-GIM (CODE/ESA/JPL/UPC) | PASS | Representative 2025-01-01, 12 comparison rows |
| Madrigal | PASS | Representative 2025-01-01 numerical equivalence |
| PyIRI climatological context | PASS | Representative 2025-01-01 hourly output equivalence |
| PyIRI strict common-hour layer | PASS | Representative 2025-01-01 exact table equality |
| Solar Cycle 25 main analysis | PASS | Complete supplied dataset; 9/9 compared CSV products exactly equal |
| Solar sensitivity analysis | PASS | Operational vs publication CSV equality |
| Lagged bootstrap analysis | PASS | Operational vs publication CSV equality |
| Monthly/figure suite | PASS | Numerical CSV products exactly equal; rendered images excluded from bit-for-bit checks |
| High-rate 10 Hz scintillation | PENDING | Implementation available separately; numerical equivalence not established |

See `VALIDATION_STATUS.md` for the consolidated record.

## Interpretation safeguards

- **PyIRI is climatological context, not an independent solar-activity validation.** Daily F10.7 is an input to the model.
- **CODE/ESA/JPL/UPC are a multi-GIM/inter-product robustness comparison**, not four independent validations.
- **Madrigal is an observation-derived external GNSS processing chain.** Its degree of independence depends on the contributing receiver network.
- The direct **IGS–Madrigal** comparison is a reference-product diagnostic.
- Solar and geomagnetic findings should be described as **associations/relationships**, not causal effects without a separate causal design.
- No empirical TEC bias correction is applied in the publication workflow.
- `S4_CNO_PROXY` from the high-rate workflow is an **uncalibrated proxy**, not a standard ISMR S4 index.

## Data and external software

Raw GNSS observations, downloaded geophysical products, caches, credentials, and generated result archives are intentionally not distributed in this source repository.

External tools/services used by parts of the workflow include GFZRNX, pyOASIS, IGS/CDDIS products, Madrigal, NASA/SPDF OMNI, WDC-SILSO, and PyIRI. Refer to the module documentation for required inputs and access requirements.

## Python environment

The scripts were developed and tested in the project Python environment. A broad dependency list is provided in `requirements.txt`; before archival publication, freeze the exact tested environment if strict environment reproduction is required.

Syntax checking for an individual script can be performed with:

```bat
python -m py_compile path\to\script.py
```

## Reproducibility and figures

Numerical validation is performed on CSV/Parquet/TEC products as appropriate. PNG/PDF/SVG files are not used for bit-for-bit equivalence because rendering metadata can vary with Matplotlib, fonts, and platform while underlying numerical data remain identical.

## Release status

This package is a **release candidate** for GitHub/Zenodo archiving. Before public release, complete the author metadata in `CITATION.cff`, choose and add a license, freeze/version the environment as desired, add the repository URL, create the release tag, and insert the Zenodo DOI after deposition. See `RELEASE_CHECKLIST.md`.
