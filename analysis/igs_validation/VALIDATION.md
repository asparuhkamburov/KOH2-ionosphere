# Validation

## Representative wrapper-equivalence test

The publication-oriented `validate_tec_igs_2019_2026.py` wrapper was validated
against the existing operational V4 workflow using:

```text
Station   : KOH2
Date      : 2025-01-01
DOY       : 001
Reference : IGS Final GIM
```

The test used the same production PyTECGg and pyOASIS outputs and the same
cached IGS Final IONEX product:

```text
IGS0OPSFIN_20250010000_01D_02H_GIM.INX
```

The publication wrapper reported:

```text
PyTECGg days found : 1
pyOASIS days found : 1
Union of days      : 1

PyTECGg VEq station matches : 2,881
PyTECGg IPP matches         : 21,007
pyOASIS IPP matches         : 27,150
```

## Bit-for-bit output equivalence

### PyTECGg VEq vs IGS at station position

```text
KOH2_2025001_PyTECGg_VEq_vs_IGS_station.parquet
SHA-256:
75025e40bc6999fe2c3a07b861efc448a3417f2554712c456ae482ac43590a4e
```

Publication and operational files were bit-for-bit identical.

### PyTECGg VTEC vs IGS at PyTECGg IPP

```text
KOH2_2025001_PyTECGg_VTEC_vs_IGS_IPP.parquet
SHA-256:
a95bcfd1ee3b5a59baca8efca93d244b0f67ee2c43280cd9048c5938c2681f81
```

Publication and operational files were bit-for-bit identical.

### pyOASIS VTEC vs IGS at pyOASIS IPP

```text
KOH2_2025001_pyOASIS_VTEC_vs_IGS_IPP.parquet
SHA-256:
005f242c2ab49010260d2e240320b15d5bad861d3f336530f0e606b82e2e3ca9
```

Publication and operational files were bit-for-bit identical.

## Validation conclusion

For the representative KOH2 day 2025-01-01, the publication-oriented wrapper
reproduces all three daily matched scientific comparison products from the
operational V4 workflow bit-for-bit.

This is a representative wrapper-equivalence spot check. It does not imply
that every observation day from 2019-2026 was independently checksum-validated,
and it is not an independent validation of the IGS Final GIM itself.
