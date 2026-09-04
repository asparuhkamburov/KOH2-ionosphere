# Release checklist

## Completed for this release candidate

- [x] Publication-oriented processing and analysis modules assembled in one repository tree.
- [x] Runtime machine-specific absolute paths removed from publication scripts.
- [x] Command-line interfaces used for portable publication workflows where required.
- [x] Python syntax validation completed for all included `.py` files.
- [x] Module-level validation documentation retained.
- [x] Consolidated validation scope documented in `VALIDATION_STATUS.md`.
- [x] Solar Cycle 25 analytical core validated against operational outputs.
- [x] Solar figure suite 01–04 numerical support tables validated against operational outputs.
- [x] Scientific wording reviewed to distinguish validation, inter-product comparison, climatological context, and association from causation.
- [x] Raw GNSS data, downloaded caches, credentials, and generated result archives excluded from the source repository.
- [x] Author metadata and ORCID added to `CITATION.cff`.
- [x] GPL-3.0-or-later selected and `LICENSE` added.
- [x] Licensing boundary documented in `LICENSE_NOTICE.md`.
- [x] Third-party license summary documented in `THIRD_PARTY_LICENSES.md`.
- [x] pyOASIS and PyTECGg wrappers marked with SPDX/provenance headers.
- [x] Workflow-specific `requirements-pyoasis.txt` and `requirements-pytecgg.txt` added.
- [x] External software citations documented in `SOFTWARE_CITATIONS.md`.
- [x] AI-assisted development disclosure added to the root `README.md`.
- [x] High-rate 10 Hz module included with explicit **EXPERIMENTAL** status and excluded from repository `PASS` claims.
- [x] High-rate README/STATUS/VALIDATION documents added with explicit warnings that `SIGMA_PHI_GF_EQUIV_RAD` is not validated `Phi60` and `S4_CNO_PROXY` is not calibrated `S4`.
- [x] Representative 2025-01-01 high-rate implementation, boundary, elevation-sensitivity, LLI, geometry-free, and ROTI/DTEC event checks documented.

## Required before public GitHub / Zenodo release

- [ ] Add the final GitHub repository URL to release metadata/documentation.
- [ ] Decide whether to freeze exact Python package versions and, if desired, generate a requirements/environment freeze from the validated environment.
- [ ] Create the final release tag/version (for example `v1.0.0`).
- [ ] Add the final `version` and `date-released` fields to `CITATION.cff`.
- [ ] Create/archive the GitHub release in Zenodo.
- [ ] Add the Zenodo DOI to `CITATION.cff` and the root `README.md`.
- [ ] Perform a final check that repository links, citations, licence references, and module paths resolve correctly.

## Third-party software and data

Before manuscript submission, verify that the bibliography contains the citations required for the external software and data products actually used in the reported analyses.

Current software citation guidance is collected in `SOFTWARE_CITATIONS.md`.
