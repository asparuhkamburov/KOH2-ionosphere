# KOH2 preprocessing

Publication-oriented wrappers for the preprocessing steps used for the KOH2
GNSS observation data set.

## Processing chain

```text
Trimble Alloy
    |
    v
hourly T02 observations
    |
    | Trimble ConvertToRINEX
    v
hourly RINEX
    |
    | GFZRNX
    | - selected header/metadata harmonization
    | - RINEX 3 naming standardization
    | - hourly -> daily concatenation
    | - 1 s daily sampling
    | - 30 s daily subsampling
    v
daily RINEX 1 s / 30 s
```

The raw KOH2 observations were stored as Trimble T02 files, one file per hour.
For the project, the observations were recorded at 1 Hz during 2019-2024 and
10 Hz during 2025-2026.

## Files

### `convert_t02_to_rinex.py`

Recursively calls the locally installed Trimble ConvertToRINEX program for
each `.T02` file. The converter output is written to the source file's
directory, reproducing the historical project wrapper.

The exact RINEX version is controlled by the installed Trimble converter
configuration. In the KOH2 project it was configured to produce RINEX 3.04.

Example:

```bat
python convert_t02_to_rinex.py ^
  --root "D:\path\to\KOH2data\2025" ^
  --converter "C:\Program Files (x86)\Trimble\ConvertToRinex\convertToRINEX.exe"
```

A non-destructive command preview is available with `--dry-run`.

### `prepare_koh2_rinex.py`

Runs the GFZRNX steps used to prepare the RINEX observation files.

For legacy `*.YYo` files, the script can apply `KOH2.crux` and create
standardized RINEX 3 output using:

```text
-fout ::RX3::00,ATA
-crux KOH2.crux
-hded
-site KOH2
```

For hourly RINEX 3 files, it creates the daily 1 s product using:

```text
-finp KOH200ATA_R_*_01H_*_MO.rnx
-fout ::RX3::00,ATA
-smp 1
-site KOH2
```

The daily 1 s product is then subsampled to 30 s using:

```text
-smp 30
```

Example for a fresh test directory containing one intended hourly stream:

```bat
python prepare_koh2_rinex.py ^
  --root "D:\path\to\test-day" ^
  --gfzrnx "D:\path\to\gfzrnx_2.2.0_win10_64.exe" ^
  --crux "KOH2.crux"
```

To reproduce the historical overwrite behavior, add:

```text
--overwrite
```

For a command-only preview, add:

```text
--dry-run
```

If a directory contains more than one hourly sampling token (for example both
`01S` and `10Z` files), the publication script skips that directory rather
than mixing two sampling streams. Use a more specific `--hourly-pattern` or
point `--root` to the intended stream.

## KOH2 CRUX metadata

`KOH2.crux` harmonizes the following fields:

- MARKER NAME: `KOH200ATA`
- MARKER NUMBER: `66026M002`
- APPROX POSITION XYZ:
  `1453335.2992 -2554570.1548 -5641700.7402` m

The CRUX file does not set antenna metadata.

## External software

These wrappers require software that is not distributed in this repository:

- Trimble ConvertToRINEX
- GFZRNX

Users are responsible for obtaining and using those programs under their
respective terms.

## Safety and provenance

The publication scripts intentionally do not reproduce the historical cleanup
operations that deleted intermediate `*.YYo` or `*.YYmix` files. Raw data and
generated products remain outside the source-code repository.

## Validation status

The hourly-to-daily preprocessing path has been validated for both principal
KOH2 sampling regimes:

- 2019-12-28 (DOY 362): hourly 1 Hz -> daily 1 s -> daily 30 s;
- 2025-01-01 (DOY 001): hourly 10 Hz -> daily 1 s -> daily 30 s.

For both validation days, the RINEX observation bodies of the generated daily
products were bit-for-bit identical to the historical operational products.

See `VALIDATION.md` for the checksums and validation scope.

Before the final monograph release, the remaining legacy `*.YYo` + CRUX
header-harmonization/RINEX-3 standardization path should also be spot-checked.
