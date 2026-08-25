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
- [x] Raw GNSS data, downloaded caches, credentials, and generated result archives excluded from the source package.

## Required before public GitHub / Zenodo release

- [ ] Replace the author placeholders in `CITATION.cff` with the intended release authors.
- [ ] Choose and add an explicit software license (`LICENSE`).
- [ ] Add the final GitHub repository URL to release metadata/documentation.
- [ ] Decide whether to freeze exact Python package versions and, if so, generate an environment lock/requirements freeze from the validated environment.
- [ ] Create the final release tag/version.
- [ ] Deposit the tagged release in Zenodo and add the DOI to `CITATION.cff` / README.
- [ ] Review third-party software/data citation requirements for GFZRNX, pyOASIS, IGS/CDDIS, Madrigal, NASA/SPDF OMNI, WDC-SILSO, and PyIRI.

## High-rate 10 Hz block

The high-rate scintillation block remains outside the validated-release claim. Either:

- complete and document its numerical-equivalence validation before promoting it to validated status; or
- retain the explicit `PENDING` status in `analysis/high_rate_10hz/STATUS.md` for the public release.
