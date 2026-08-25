# KOH2 GNSS Ionospheric Processing

Reproducible processing scripts used for the analysis of GNSS-derived
ionospheric TEC at station **KOH2**, Livingston Island, Antarctica.

This repository is intended to accompany a scientific monograph and related
publications. The repository should contain the authors' processing and
analysis code, but **not** redistribute third-party executables, third-party
datasets, Earthdata credentials, or the full GNSS observation archive.

## Processing chain

The observational workflow is:

```text
Trimble Alloy
    |
    +-- T02 hourly raw files
    |
    +-- Trimble Convert To RINEX
    |
    +-- RINEX 3.04
    |
    +-- GFZRNX + project scripts
    |      - header harmonization
    |      - file naming
    |      - concatenation to daily files
    |      - subsampling where required
    |      - RINEX-version conversion where required
    |
    +-- PyTECGg / pyOASIS
    |
    +-- TEC validation and Solar Cycle 25 analysis
```

## `process_koh2_pyoasis.py`

This is the publication-oriented version of the annual pyOASIS automation
script. It reproduces the operational processing sequence used by the project:

1. traverse the `YEAR/MONTH/DAY` archive;
2. select a suitable daily/partial-day RINEX observation file;
3. obtain the corresponding GFZ MGEX rapid SP3 orbit from CDDIS;
4. prepare a short-name RINEX 2.x staging file with GFZRNX;
5. run:
   - `SP3intp`
   - `RNXclean`
   - `RNXlevelling`
   - `ROTIcalc`
   - `DTECcalc`
   - `SIDXcalc`
   - `TECcalc`;
6. skip a day when the expected final TEC output is already present, unless
   `--force` is supplied.

### Example on Windows

```bat
python -m py_compile process_koh2_pyoasis.py

python process_koh2_pyoasis.py ^
  --year 2023 ^
  --data-root "E:\KOH2data\pyOASIS" ^
  --gfzrnx "D:\GNSS\gfzrnx_2.2.0_win10_64.exe"
```

`--data-root` may also point directly to the annual folder:

```bat
python process_koh2_pyoasis.py ^
  --year 2023 ^
  --data-root "E:\KOH2data\pyOASIS\2023" ^
  --gfzrnx "D:\GNSS\gfzrnx_2.2.0_win10_64.exe"
```

### Environment-variable alternative

```bat
set KOH2_DATA_ROOT=E:\KOH2data\pyOASIS
set GFZRNX_PATH=D:\GNSS\gfzrnx_2.2.0_win10_64.exe

python process_koh2_pyoasis.py --year 2023
```

## External dependencies

The workflow relies on external software/data services that are **not**
distributed in this repository:

- Trimble Convert To RINEX
- GFZRNX
- pyOASIS
- NASA CDDIS / Earthdata access for precise GNSS products

Earthdata authentication should be configured outside the source code, for
example through the user's `.netrc` / `_netrc` file.

## Data organization

Expected annual structure:

```text
<DATA_ROOT>\
  2023\
    01\
      01\
        RINEX\
        PRODUCTS\
        pyOASIS_INPUT\
        pyOASIS_OUTPUT\
      02\
        ...
    02\
      ...
```

Only the `RINEX` observation input is required initially. Product/input/output
folders are created as needed.


## Validation of the refactored workflow

The publication-oriented script was functionally checked against an existing
operational pyOASIS result for **KOH2, 2025-01-01 (DOY 001)**.

The refactored workflow completed the full processing chain and generated:

```text
KOH2_001_2025_L1L2.TEC
```

The final TEC file was compared with the previously generated operational
result using SHA-256. Both files produced the identical checksum:

```text
4280d8a6451adc158cb0b37f350f5b67ef45482587d57e6f6b7096dcfac5634f
```

This confirms bit-for-bit identity of the final TEC output for the validation
day.

The script also supports a restricted single-day test through:

```bat
python process_koh2_pyoasis.py ^
  --year 2025 ^
  --date 2025-01-01 ^
  --data-root "D:\path\to\data" ^
  --gfzrnx "D:\path\to\gfzrnx.exe"
```

## Reproducibility

Before creating the monograph release:

1. verify the code against the final processing environment;
2. record exact package versions;
3. create a tagged GitHub release, e.g. `v1.0-monograph`;
4. archive that release in Zenodo;
5. cite the resulting DOI in the monograph.

## Citation

A `CITATION.cff` template is included. Replace the author placeholders and,
after archiving the release, add the final DOI and repository URL.

## License

Choose and add a software license before public release. MIT or BSD-3-Clause
are common permissive choices, but the final choice should reflect the
authors' and institution's requirements.
