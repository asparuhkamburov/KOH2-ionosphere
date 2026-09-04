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
├── SOFTWARE_CITATIONS.md
├── CITATION.cff
├── LICENSE
├── LICENSE_NOTICE.md
├── THIRD_PARTY_LICENSES.md
├── requirements.txt
├── requirements-pyoasis.txt
├── requirements-pytecgg.txt
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

The separate `analysis/high_rate_10hz/` branch processes 10 Hz RINEX observations as an **EXPERIMENTAL** geometry-free phase-fluctuation workflow. It is not included in the validated TEC `PASS` claims.

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
| High-rate 10 Hz phase-fluctuation analysis | **EXPERIMENTAL** | 2025-01-01 implementation/QC/sensitivity/event-level checks; not validated Phi60/S4 |

See `VALIDATION_STATUS.md` for the consolidated record.

## Interpretation safeguards

- **PyIRI is climatological context, not an independent solar-activity validation.** Daily F10.7 is an input to the model.
- **CODE/ESA/JPL/UPC are a multi-GIM/inter-product robustness comparison**, not four independent validations.
- **Madrigal is an observation-derived external GNSS processing chain.** Its degree of independence depends on the contributing receiver network.
- The direct **IGS–Madrigal** comparison is a reference-product diagnostic.
- Solar and geomagnetic findings should be described as **associations/relationships**, not causal effects without a separate causal design.
- No empirical TEC bias correction is applied in the publication workflow.
- The high-rate `SIGMA_PHI_GF_EQUIV_RAD` product is an **experimental dual-frequency geometry-free phase-fluctuation proxy**, not a validated ISMR `Phi60` index.
- `S4_CNO_PROXY` is an **uncalibrated C/N0-derived proxy**, not calibrated or standard ISMR `S4`.
- The earlier single-frequency high-rate `SIGMA_PHI_RAD` is retained only as a diagnostic because cross-frequency testing showed substantial nondispersive/common-mode content.

## Data and external software

Raw GNSS observations, downloaded geophysical products, caches, credentials, and generated result archives are intentionally not distributed in this source repository.

Parts of the workflow rely on third-party scientific software, including GFZRNX, pyOASIS/OASIS, and PyTECGg. These packages and executables are not redistributed by this repository unless explicitly stated. Users should obtain them from their official distribution channels and comply with their respective licence terms.

External scientific software used to produce results should be cited separately from this repository. The recommended citations and project links are collected in `SOFTWARE_CITATIONS.md`.

External data products and services used by parts of the workflow include IGS/CDDIS products, Madrigal, NASA/SPDF OMNI, and WDC-SILSO. Refer to the module documentation for required inputs and access requirements.

## External software citations

When publishing results produced with this workflow, please cite the relevant third-party software in addition to citing this repository.

Key software citations include:

- **GFZRNX** — Nischan, T. (2016). *GFZRNX - RINEX GNSS Data Conversion and Manipulation Toolbox*. GFZ Data Services. https://doi.org/10.5880/GFZ.1.1.2016.002
- **PyTECGg** — Ventriglia, V., Guerra, M., Okoh, D., Vermicelli, P., Ciraolo, L., & Cesaroni, C. (2026). *PyTECGg: total electron content calibration with GNSS data*. SoftwareX, 34, 102737. https://doi.org/10.1016/j.softx.2026.102737
- **pyOASIS / OASIS** — use the current citation instructions supplied by the OASIS project. At the time this repository documentation was prepared, the project listed the manuscript *Introducing OASIS: An Open-Access System for Ionospheric Studies* as submitted, together with a software-repository citation.

See `SOFTWARE_CITATIONS.md` for full details.

## Python environment

The scripts were developed and tested in project Python environments. Workflow-specific direct dependencies for the two principal TEC-processing wrappers are provided in `requirements-pyoasis.txt` and `requirements-pytecgg.txt`. The root `requirements.txt` is retained as a backward-compatible convenience entry point for the pyOASIS wrapper. Analysis modules may require additional packages documented in their module READMEs. Before archival publication, freeze the exact tested environment(s) if strict environment reproduction is required.

Syntax checking for an individual script can be performed with:

```bat
python -m py_compile path\to\script.py
```

## Reproducibility and figures

Numerical validation is performed on CSV/Parquet/TEC products as appropriate. PNG/PDF/SVG files are not used for bit-for-bit equivalence because rendering metadata can vary with Matplotlib, fonts, and platform while underlying numerical data remain identical.

## AI-assisted development

Parts of the source code and documentation in this repository were developed with assistance from OpenAI ChatGPT.

The scientific workflow, processing parameters, validation design, interpretation, testing, and final review were performed by the repository author. AI-generated or AI-assisted code was reviewed, modified where necessary, and validated against the operational processing results before inclusion in the publication-oriented repository.

The author takes responsibility for the released software and its scientific use.

## License and third-party software

Except where otherwise noted, original source code authored specifically for the KOH2-ionosphere workflow is released under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

Copyright (C) 2026 Asparuh Kamburov.

The repository uses separately distributed third-party scientific software that retains its own licensing terms and is not relicensed by KOH2-ionosphere. In particular:

- **PyTECGg** is an external Python dependency licensed upstream under **GPL-3.0-or-later**.
- **OASIS / pyOASIS** is an external Python dependency licensed upstream under **CC BY-NC 4.0**. The upstream non-commercial-use condition remains applicable to pyOASIS itself.
- **GFZRNX** is an external executable distributed by GFZ under its own license terms; it is not included in this repository.

The KOH2 `pyoasis` and `pytecgg` scripts are workflow/orchestration wrappers. The current repository does not vendor the upstream source code of OASIS/pyOASIS or PyTECGg.

See [`LICENSE`](LICENSE) for the GPLv3 text, [`LICENSE_NOTICE.md`](LICENSE_NOTICE.md) for the licensing boundary, [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) for the third-party license summary, and [`SOFTWARE_CITATIONS.md`](SOFTWARE_CITATIONS.md) for scientific attribution and citation guidance.

## Citing this repository

The `CITATION.cff` file contains citation metadata for this repository itself.

Third-party scientific software used by the workflow must be cited separately; see `SOFTWARE_CITATIONS.md`.

## Release status

This repository is being prepared for GitHub/Zenodo archival release. The high-rate 10 Hz module is included explicitly as EXPERIMENTAL and remains outside the validated `PASS` claims. Author metadata, third-party software citations and licensing boundaries, AI-assisted development disclosure, and the GPL-3.0-or-later licence are documented.

Before archival publication, add the final GitHub repository URL, create the release tag/version, add the release date/version metadata to `CITATION.cff`, and insert the Zenodo DOI after deposition. See `RELEASE_CHECKLIST.md`.
