# External Scientific Software and Recommended Citations

This repository contains publication-oriented scripts developed for the KOH2 GNSS ionospheric TEC study. Several stages of the workflow rely on third-party scientific software.

These external software packages are not authored by the maintainer of this repository and should be cited separately when results produced with them are used in a publication.

Third-party software is not redistributed by this repository unless explicitly stated. Users should obtain each package or executable from its official source and comply with the corresponding licence terms.

## GFZRNX

**Role in the KOH2 workflow**

GFZRNX is used for RINEX manipulation tasks such as format conversion, file concatenation/splicing, sampling, header editing/harmonisation, metadata inspection, and preparation of observation files for subsequent processing.

**Recommended citation**

Nischan, T. (2016). *GFZRNX - RINEX GNSS Data Conversion and Manipulation Toolbox*. GFZ Data Services.  
https://doi.org/10.5880/GFZ.1.1.2016.002

**Official project page**

https://gnss.gfz.de/services/gfzrnx

GFZ explicitly requests this citation when GFZRNX is used.

---

## PyTECGg

**Role in the KOH2 workflow**

PyTECGg is used for GNSS Total Electron Content processing, including satellite-level TEC estimation and station-level equivalent vertical TEC products used in the KOH2 comparison and validation workflow.

**Recommended publication citation**

Ventriglia, V., Guerra, M., Okoh, D., Vermicelli, P., Ciraolo, L., & Cesaroni, C. (2026). *PyTECGg: total electron content calibration with GNSS data*. SoftwareX, 34, 102737.  
https://doi.org/10.1016/j.softx.2026.102737

**Package/project**

https://pypi.org/project/pytecgg/

The software publication should be cited when PyTECGg-derived results are reported.

---

## pyOASIS / OASIS

**Role in the KOH2 workflow**

pyOASIS (Open-Access System for Ionospheric Studies) is used as an external GNSS ionospheric-processing package in the KOH2 workflow. Its processing chain is used to derive TEC and related ionospheric indices from GNSS observations.

**Current citation guidance from the OASIS project**

At the time this file was prepared, the OASIS project listed the following manuscript as **submitted**:

Picanço, G. A. S., Fagundes, P. R., Prol, F. S., Denardini, C. M., Mendoza, L. P. O., Pillat, V. G., Rodrigues, I., Christovam, A. L., Meza, A. M., Natali, M. P., Romero-Hernández, E., Aguirre-Gutierrez, R., Agyei-Yeboah, E., & Muella, M. T. A. H. (2025). *Introducing OASIS: An Open-Access System for Ionospheric Studies*. GPS Solutions. Manuscript submitted for publication.

The project also supplies a software-repository citation:

Picanço, G. A. S. (2025). *OASIS: Open-Access System for Ionospheric Studies* [Software]. GitHub.  
https://github.com/giorgiopicanco/OASIS

**Official package/project pages**

https://pypi.org/project/pyOASIS/  
https://github.com/giorgiopicanco/OASIS

Because the project currently labels the manuscript as submitted, users should check the official OASIS project page before final manuscript submission and replace the provisional reference if a final journal citation or DOI has become available.

---

## How these citations relate to CITATION.cff

`CITATION.cff` describes how to cite the **KOH2 repository itself**.

It does not replace citations to the external scientific software used by the workflow.

For a publication based on this repository, the recommended practice is therefore to cite:

1. the KOH2 repository/software release;
2. GFZRNX where GFZRNX was used;
3. PyTECGg where PyTECGg-derived results were used;
4. pyOASIS/OASIS where pyOASIS-derived results were used;
5. the external scientific data products and services relevant to the specific analysis.

## Suggested wording for a Methods section

A concise acknowledgement in a manuscript can be written along the following lines:

> GNSS observation files were prepared and harmonised using GFZRNX. Ionospheric TEC was independently processed using pyOASIS/OASIS and PyTECGg, followed by the validation and analysis procedures implemented in the KOH2 publication workflow. The corresponding software packages are cited separately in the reference list.

The exact wording should be adapted to the processing steps actually used in the publication.
