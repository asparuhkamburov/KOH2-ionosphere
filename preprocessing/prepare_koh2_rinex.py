#!/usr/bin/env python3
"""
Prepare KOH2 RINEX observation files with GFZRNX.

Publication-oriented consolidation of the GFZRNX operations used during the
KOH2 preprocessing workflow:

1. Legacy RINEX observation files (*.YYo) can be converted/standardized to
   RINEX 3 while selected station metadata are harmonized from KOH2.crux.
2. Hourly RINEX 3 observation files are concatenated and sampled to 1 s.
3. The resulting daily 1 s RINEX is subsampled to 30 s.

The numerical GFZRNX operations mirror the historical project commands:
    -fout ::RX3::00,ATA
    -site KOH2
    -crux KOH2.crux -hded          (legacy header harmonization)
    -smp 1                         (daily 1 s product)
    -smp 30                        (daily 30 s product)

The script intentionally performs no raw-data deletion or folder reorganization.
"""

from __future__ import annotations

import argparse
import fnmatch
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


LEGACY_RINEX_RE = re.compile(r"\.\d{2}o$", re.IGNORECASE)
SAMPLING_TOKEN_RE = re.compile(r"_01H_([^_]+)_", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    default_crux = Path(__file__).with_name("KOH2.crux")

    parser = argparse.ArgumentParser(
        description=(
            "Harmonize legacy KOH2 RINEX files and create daily 1 s / 30 s "
            "RINEX products with GFZRNX."
        )
    )
    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        help="Root directory to scan recursively.",
    )
    parser.add_argument(
        "--gfzrnx",
        required=True,
        type=Path,
        help="Path to the GFZRNX executable.",
    )
    parser.add_argument(
        "--crux",
        type=Path,
        default=default_crux,
        help=f"CRUX file for header harmonization (default: {default_crux.name}).",
    )
    parser.add_argument(
        "--station",
        default="KOH2",
        help="Four-character station code passed to GFZRNX (default: KOH2).",
    )
    parser.add_argument(
        "--rinex3-suffix",
        default="ATA",
        help="RINEX 3 naming suffix used in ::RX3:: output (default: ATA).",
    )
    parser.add_argument(
        "--hourly-pattern",
        default="KOH200ATA_R_*_01H_*_MO.rnx",
        help=(
            "Hourly RINEX 3 filename pattern used for daily concatenation. "
            "Default reproduces the historical KOH2 pattern."
        ),
    )
    parser.add_argument(
        "--skip-legacy-header",
        action="store_true",
        help="Do not process legacy *.YYo observation files.",
    )
    parser.add_argument(
        "--skip-daily-1s",
        action="store_true",
        help="Do not concatenate hourly files to a daily 1 s file.",
    )
    parser.add_argument(
        "--skip-daily-30s",
        action="store_true",
        help="Do not create the daily 30 s product.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Pass -f to GFZRNX, allowing generated outputs to be overwritten. "
            "The historical project batch scripts used -f."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show GFZRNX commands without executing them.",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop after the first GFZRNX failure.",
    )
    return parser.parse_args()


def gfzrnx_output_spec(suffix: str) -> str:
    return f"::RX3::00,{suffix}"


def run_gfzrnx(
    executable: Path,
    args: list[str],
    cwd: Path,
    dry_run: bool,
) -> bool:
    command = [str(executable), *args]

    print("  CMD:", subprocess.list2cmdline(command))

    if dry_run:
        return True

    try:
        result = subprocess.run(
            command,
            cwd=str(cwd),
            check=True,
            text=True,
            capture_output=True,
        )
        if result.stdout.strip():
            print(result.stdout.rstrip())
        if result.stderr.strip():
            print(result.stderr.rstrip(), file=sys.stderr)
        return True

    except subprocess.CalledProcessError as exc:
        print(f"  ERROR: GFZRNX exit code {exc.returncode}", file=sys.stderr)
        if exc.stdout:
            print(exc.stdout.rstrip())
        if exc.stderr:
            print(exc.stderr.rstrip(), file=sys.stderr)
        return False


def find_legacy_observation_files(root: Path) -> list[Path]:
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and LEGACY_RINEX_RE.search(path.name)
    )


def find_hourly_directories(root: Path, pattern: str) -> dict[Path, list[Path]]:
    grouped: dict[Path, list[Path]] = defaultdict(list)

    for path in root.rglob("*.rnx"):
        if path.is_file() and fnmatch.fnmatch(path.name.lower(), pattern.lower()):
            grouped[path.parent].append(path)

    return {
        folder: sorted(files)
        for folder, files in sorted(grouped.items(), key=lambda item: str(item[0]))
    }


def sampling_tokens(files: list[Path]) -> set[str]:
    tokens: set[str] = set()

    for path in files:
        match = SAMPLING_TOKEN_RE.search(path.name)
        if match:
            tokens.add(match.group(1).upper())

    return tokens


def daily_1s_files(folder: Path) -> list[Path]:
    return sorted(
        path
        for path in folder.glob("*_01D_01S_*.rnx")
        if path.is_file()
    )


def legacy_header_stage(args: argparse.Namespace) -> tuple[int, int]:
    if args.skip_legacy_header:
        print("Legacy header stage: SKIPPED")
        return 0, 0

    crux = args.crux.resolve()
    if not crux.is_file():
        print(f"ERROR: CRUX file does not exist: {crux}", file=sys.stderr)
        return 0, 1

    legacy_files = find_legacy_observation_files(args.root.resolve())

    print()
    print("=" * 80)
    print("STAGE 1 - LEGACY RINEX HEADER HARMONIZATION / RINEX 3 STANDARDIZATION")
    print("=" * 80)
    print(f"Legacy observation files: {len(legacy_files)}")

    ok_count = 0
    fail_count = 0

    for path in legacy_files:
        print()
        print("Processing:", path)

        cmd_args = [
            "-finp",
            path.name,
            "-fout",
            gfzrnx_output_spec(args.rinex3_suffix),
            "-crux",
            str(crux),
            "-hded",
            "-site",
            args.station,
        ]
        if args.overwrite:
            cmd_args.append("-f")

        ok = run_gfzrnx(
            args.gfzrnx,
            cmd_args,
            path.parent,
            args.dry_run,
        )

        if ok:
            ok_count += 1
        else:
            fail_count += 1
            if args.stop_on_error:
                break

    return ok_count, fail_count


def daily_stage(args: argparse.Namespace) -> tuple[int, int, int]:
    root = args.root.resolve()
    grouped = find_hourly_directories(root, args.hourly_pattern)

    print()
    print("=" * 80)
    print("STAGE 2 - HOURLY -> DAILY 1 s -> DAILY 30 s")
    print("=" * 80)
    print(f"Hourly-data directories: {len(grouped)}")

    processed_dirs = 0
    fail_count = 0
    ambiguous_dirs = 0

    for folder, hourly_files in grouped.items():
        print()
        print("-" * 80)
        print("Folder:", folder)
        print("Hourly files:", len(hourly_files))

        tokens = sampling_tokens(hourly_files)
        if tokens:
            print("Hourly sampling token(s):", ", ".join(sorted(tokens)))

        # Prevent accidental mixing of two different hourly sampling streams.
        if len(tokens) > 1:
            ambiguous_dirs += 1
            print(
                "  SKIPPED: multiple hourly sampling rates were found in the "
                "same directory. Point --root to the intended stream or use a "
                "more specific --hourly-pattern.",
                file=sys.stderr,
            )
            continue

        processed_dirs += 1

        if not args.skip_daily_1s:
            print("Creating daily 1 s RINEX...")

            cmd_args = [
                "-finp",
                args.hourly_pattern,
                "-fout",
                gfzrnx_output_spec(args.rinex3_suffix),
                "-smp",
                "1",
                "-site",
                args.station,
            ]
            if args.overwrite:
                cmd_args.append("-f")

            ok = run_gfzrnx(
                args.gfzrnx,
                cmd_args,
                folder,
                args.dry_run,
            )
            if not ok:
                fail_count += 1
                if args.stop_on_error:
                    break
                continue

        if args.skip_daily_30s:
            continue

        if args.dry_run:
            print(
                "DRY RUN: daily 30 s stage would use the daily "
                "*_01D_01S_*.rnx product created above."
            )
            continue

        one_second_files = daily_1s_files(folder)

        if not one_second_files:
            print(
                "  ERROR: no *_01D_01S_*.rnx file found after daily 1 s stage.",
                file=sys.stderr,
            )
            fail_count += 1
            if args.stop_on_error:
                break
            continue

        for daily_path in one_second_files:
            print("Creating daily 30 s RINEX from:", daily_path.name)

            cmd_args = [
                "-finp",
                daily_path.name,
                "-fout",
                gfzrnx_output_spec(args.rinex3_suffix),
                "-smp",
                "30",
                "-site",
                args.station,
            ]
            if args.overwrite:
                cmd_args.append("-f")

            ok = run_gfzrnx(
                args.gfzrnx,
                cmd_args,
                folder,
                args.dry_run,
            )

            if not ok:
                fail_count += 1
                if args.stop_on_error:
                    return processed_dirs, fail_count, ambiguous_dirs

    return processed_dirs, fail_count, ambiguous_dirs


def main() -> int:
    args = parse_args()

    args.root = args.root.resolve()
    args.gfzrnx = args.gfzrnx.resolve()
    args.crux = args.crux.resolve()

    if not args.root.is_dir():
        print(f"ERROR: root directory does not exist: {args.root}", file=sys.stderr)
        return 2

    if not args.gfzrnx.is_file():
        print(f"ERROR: GFZRNX executable does not exist: {args.gfzrnx}", file=sys.stderr)
        return 2

    print("=" * 80)
    print("KOH2 GFZRNX PREPROCESSING")
    print("=" * 80)
    print(f"Root            : {args.root}")
    print(f"GFZRNX          : {args.gfzrnx}")
    print(f"CRUX            : {args.crux}")
    print(f"Station         : {args.station}")
    print(f"Hourly pattern  : {args.hourly_pattern}")
    print(f"Overwrite (-f)  : {args.overwrite}")
    print(f"Dry run         : {args.dry_run}")

    legacy_ok, legacy_failed = legacy_header_stage(args)

    if legacy_failed and args.stop_on_error:
        return 1

    daily_dirs, daily_failed, ambiguous_dirs = daily_stage(args)

    print()
    print("=" * 80)
    print("PREPROCESSING SUMMARY")
    print("=" * 80)
    print(f"Legacy files processed successfully : {legacy_ok}")
    print(f"Legacy failures                     : {legacy_failed}")
    print(f"Hourly-data directories processed   : {daily_dirs}")
    print(f"Ambiguous hourly directories skipped: {ambiguous_dirs}")
    print(f"Daily-stage failures                : {daily_failed}")

    failures = legacy_failed + daily_failed
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
