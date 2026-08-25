# KOH2 pyOASIS processing workflow

This directory contains the publication-oriented pyOASIS processing wrapper
used for KOH2 GNSS ionospheric processing.

## Script

```text
process_koh2_pyoasis.py
```

The wrapper automates the daily preparation and pyOASIS processing chain for a
KOH2 year/month/day archive.

## Processing workflow

For each valid `YEAR/MONTH/DAY` directory, the script:

1. searches the day's `RINEX` directory for a suitable observation file;
2. selects the pyOASIS input using the following priority:
   - full-day 30 s mixed RINEX;
   - another available 30 s mixed RINEX segment;
   - full-day 15 s mixed RINEX;
3. obtains the corresponding GFZ MGEX rapid precise orbit product:
   `GFZ0MGXRAP_YYYYDDD0000_01D_05M_ORB.SP3.gz`;
4. extracts the SP3 product when required;
5. converts/prepares a short-name RINEX 2.x staging observation file using
   GFZRNX;
6. runs the pyOASIS processing sequence:

```text
SP3intp
-> RNXclean
-> RNXlevelling
-> ROTIcalc
-> DTECcalc
-> SIDXcalc
-> TECcalc
```

7. checks for the expected final TEC product.

The script skips a day when the expected final TEC file already exists unless
`--force` is supplied.

## Expected archive structure

Typical daily layout:

```text
YEAR/
└── MM/
    └── DD/
        ├── RINEX/
        ├── PRODUCTS/
        │   ├── SP3/
        │   ├── CLK/
        │   └── NAV/
        ├── pyOASIS_INPUT/
        │   └── RINEX/
        └── pyOASIS_OUTPUT/
            ├── ORBITS/
            ├── RINEX/
            └── INDICES/
                ├── ROTI/
                ├── DTEC/
                ├── SIDX/
                └── TEC/
```

The expected final TEC product is:

```text
pyOASIS_OUTPUT/INDICES/TEC/KOH2_DDD_YYYY_L1L2.TEC
```

## Command-line use

Example:

```bat
python process_koh2_pyoasis.py ^
 --year 2025 ^
 --data-root "F:\KOH2data" ^
 --gfzrnx "D:\MGU\Projects\AntarcticaGNSS\KOH2data\gfzrnx_2.2.0_win10_64.exe"
```

The data root may be either the parent directory containing annual folders or
the requested annual folder itself.

Important options include:

```text
--year
--data-root
--gfzrnx
--station
--cddis-base
--timeout
--delete-gz
--force
```

`KOH2_DATA_ROOT` and `GFZRNX_PATH` may also be supplied as environment
variables.

## External dependencies

The workflow requires:

- Python;
- pyOASIS;
- GFZRNX;
- `requests`;
- a working CDDIS / NASA Earthdata configuration when SP3 files must be
  downloaded.

GFZRNX and pyOASIS are external software and are not redistributed in this
repository.

## Scientific role

This script is a processing wrapper. It does not redefine the scientific
algorithms implemented by pyOASIS; it automates input selection, orbit
preparation, file staging and execution of the standard pyOASIS processing
functions.

## Validation

The publication workflow passed a representative bit-for-bit output
equivalence check for KOH2 on 2025-01-01.

See [`VALIDATION.md`](VALIDATION.md) for the validation scope and checksum.
