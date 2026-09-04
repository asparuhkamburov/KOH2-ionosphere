# KOH2-ionosphere License Notice

Copyright (C) 2026 Asparuh Kamburov

Except where otherwise noted, original source code developed specifically for the **KOH2-ionosphere** project is licensed under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

The complete GPLv3 license text is provided in the repository's [`LICENSE`](LICENSE) file.

## Scope of the project license

The GPL-3.0-or-later declaration applies to original KOH2-ionosphere source code for which the project author holds the applicable rights. It does not replace, modify, or supersede the licenses or terms of independently developed third-party software used by the workflow.

The current repository does **not** vendor the source code of pyOASIS/OASIS or PyTECGg. The KOH2 processing wrappers use separately installed third-party packages through their Python APIs. GFZRNX is invoked as an external executable and is not distributed with this repository.

If third-party source code is added to the repository in a future version, that material must retain the licensing and attribution requirements applicable to the corresponding upstream project and must be identified separately.

## pyOASIS / OASIS

Parts of the KOH2-ionosphere workflow use **OASIS / pyOASIS** as an external Python dependency.

- Upstream project: https://github.com/giorgiopicanco/OASIS
- Upstream license: **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)**

KOH2-ionosphere does not relicense OASIS/pyOASIS. The `pyoasis/process_koh2_pyoasis.py` wrapper contains KOH2-specific workflow and orchestration code and calls the separately installed pyOASIS package through its public Python API; the pyOASIS source code is not included in that file or elsewhere in this repository.

Users of the pyOASIS processing path remain responsible for complying with the upstream OASIS/pyOASIS license, including its attribution and non-commercial-use conditions.

## PyTECGg

Parts of the KOH2-ionosphere workflow use **PyTECGg** as an external Python dependency.

- Upstream project: https://github.com/viventriglia/PyTECGg
- Upstream license: **GNU General Public License v3.0 or later (GPL-3.0-or-later)**

KOH2-ionosphere does not relicense PyTECGg. The `pytecgg/process_koh2_pytecgg.py` wrapper contains KOH2-specific workflow and orchestration code and calls the separately installed PyTECGg package through its Python API; the PyTECGg source code is not vendored in this repository.

## GFZRNX

Some preprocessing and pyOASIS workflow steps invoke **GFZRNX** as an external executable.

GFZRNX is not distributed with KOH2-ionosphere. GFZ provides separate licensing terms, including a scientific-use license and commercial licensing options. Users must obtain GFZRNX from its official distribution channel and comply with the license applicable to their use.

Official service page: https://gnss.gfz.de/services/gfzrnx

## Other dependencies, data products, and services

Other Python libraries, command-line utilities, GNSS products, data services, models, and archives used by KOH2-ionosphere remain subject to their respective upstream licenses, terms of use, access conditions, and attribution requirements.

See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) and [`SOFTWARE_CITATIONS.md`](SOFTWARE_CITATIONS.md) for additional information.

## Scientific attribution

Software licensing and scientific citation are separate requirements. Users should cite KOH2-ionosphere and the relevant third-party scientific software, algorithms, GNSS products, models, and datasets according to the recommendations of their respective authors and providers.

This notice is intended to document software provenance and licensing boundaries; it is not legal advice.
