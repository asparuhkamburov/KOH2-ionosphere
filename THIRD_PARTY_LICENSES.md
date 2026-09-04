# Third-Party Software and License Summary

This document identifies principal third-party software used by the **KOH2-ionosphere** workflow. It is an informational summary. The authoritative terms are the licenses and conditions published by the respective upstream projects.

Third-party software is not automatically covered by the KOH2-ionosphere GPL-3.0-or-later declaration merely because a project script imports or invokes it.

## 1. OASIS / pyOASIS

**Role in KOH2-ionosphere:** GNSS ionospheric processing used by the pyOASIS workflow, including the standard pyOASIS processing functions called by the KOH2 wrapper.

**Distribution model in this repository:** external dependency; upstream source code is not vendored or redistributed.

**Upstream project:** https://github.com/giorgiopicanco/OASIS

**Upstream license:** Creative Commons Attribution-NonCommercial 4.0 International (**CC BY-NC 4.0**).

The KOH2 wrapper imports and calls the separately installed pyOASIS package. The OASIS/pyOASIS license remains applicable to that third-party software. In particular, users must observe the upstream attribution requirements and non-commercial-use restriction.

The repository's GPL-3.0-or-later license does not relicense OASIS/pyOASIS and does not remove its upstream restrictions.

---

## 2. PyTECGg

**Role in KOH2-ionosphere:** GNSS Total Electron Content processing, including geometry, arc extraction, calibrated TEC, and vertical-equivalent TEC products used by the KOH2 workflow.

**Distribution model in this repository:** external dependency; upstream source code is not vendored or redistributed.

**Upstream project:** https://github.com/viventriglia/PyTECGg

**Upstream license:** GNU General Public License v3.0 or later (**GPL-3.0-or-later**).

The KOH2 PyTECGg wrapper imports and calls the separately installed PyTECGg package. PyTECGg retains its upstream copyright and licensing terms.

---

## 3. GFZRNX

**Role in KOH2-ionosphere:** RINEX conversion, manipulation, concatenation, sampling, metadata harmonisation, and staging for downstream processing.

**Distribution model in this repository:** external executable; not redistributed.

**Official service page:** https://gnss.gfz.de/services/gfzrnx

**Licensing:** GFZ distributes GFZRNX under its own proprietary license terms. GFZ currently provides a scientific-use license for eligible scientific users and separate commercial licensing for commercial use. Registration and compliance with the applicable GFZRNX terms are the user's responsibility.

KOH2-ionosphere does not grant any rights to GFZRNX.

---

## 4. Python scientific ecosystem

The repository also uses third-party Python packages. Depending on the workflow, these include packages such as:

- NumPy
- pandas
- SciPy
- Matplotlib
- Polars
- Requests
- h5py
- PyIRI
- MadrigalWeb
- unlzw3
- urllib3
- and dependencies installed transitively by pyOASIS and PyTECGg.

These are independent third-party projects and retain their respective upstream licenses. This file does not reproduce every dependency's license text; users who redistribute a complete environment, container, wheel bundle, or binary distribution should review the license metadata of the exact package versions being redistributed.

Workflow-specific installation files are provided where appropriate, but they are dependency specifications rather than relicensing instruments.

---

## 5. External data products, models, and services

Parts of the workflow use external GNSS and geophysical resources, including products or services from IGS/CDDIS, BKG, Madrigal, NASA/SPDF OMNI, WDC-SILSO, and other providers documented in the module READMEs.

Data licenses, access policies, acknowledgement requirements, and citation requirements are separate from the KOH2-ionosphere software license and remain applicable.

---

## 6. KOH2-ionosphere original code

Except where otherwise indicated, original source code authored specifically for **KOH2-ionosphere** is released under:

**GNU General Public License v3.0 or later (GPL-3.0-or-later)**

See [`LICENSE`](LICENSE) and [`LICENSE_NOTICE.md`](LICENSE_NOTICE.md).

Copyright (C) 2026 Asparuh Kamburov

---

## Disclaimer

This summary is provided to improve transparency about software provenance and third-party dependencies. It is not legal advice. If this summary conflicts with an authoritative upstream license or set of terms, the authoritative upstream terms govern the corresponding third-party component.
