#!/usr/bin/env python3
"""
Convert Trimble T02 observation files with Trimble ConvertToRINEX.

This is a publication-oriented wrapper around the conversion step used for the
KOH2 GNSS data set. It does not bundle or replace Trimble ConvertToRINEX.

The historical KOH2 workflow used one T02 file per hour and wrote the converter
output into the same directory as the source file.

Note
----
The exact RINEX version is controlled by the installed Trimble ConvertToRINEX
configuration. In the KOH2 project the converter was configured to produce
RINEX 3.04.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recursively convert Trimble T02 files using ConvertToRINEX."
    )
    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        help="Root directory containing T02 files.",
    )
    parser.add_argument(
        "--converter",
        required=True,
        type=Path,
        help="Path to Trimble convertToRINEX.exe.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show commands without executing the converter.",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop immediately if one conversion fails.",
    )
    return parser.parse_args()


def find_t02_files(root: Path) -> list[Path]:
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() == ".t02"
    )


def main() -> int:
    args = parse_args()

    root = args.root.resolve()
    converter = args.converter.resolve()

    if not root.is_dir():
        print(f"ERROR: root directory does not exist: {root}", file=sys.stderr)
        return 2

    if not converter.is_file():
        print(f"ERROR: converter executable does not exist: {converter}", file=sys.stderr)
        return 2

    t02_files = find_t02_files(root)

    print("=" * 80)
    print("KOH2 T02 -> RINEX CONVERSION")
    print("=" * 80)
    print(f"Root      : {root}")
    print(f"Converter : {converter}")
    print(f"T02 files : {len(t02_files)}")
    print(f"Dry run   : {args.dry_run}")
    print()

    if not t02_files:
        print("No T02 files found.")
        return 0

    converted = 0
    failed = 0

    for index, t02_path in enumerate(t02_files, start=1):
        output_dir = t02_path.parent

        command = [
            str(converter),
            str(t02_path),
            "-p",
            str(output_dir),
        ]

        print(f"[{index}/{len(t02_files)}] {t02_path}")
        print("  Output directory:", output_dir)

        if args.dry_run:
            print("  DRY RUN:", subprocess.list2cmdline(command))
            converted += 1
            continue

        try:
            result = subprocess.run(
                command,
                check=True,
                text=True,
                capture_output=True,
            )
            if result.stdout.strip():
                print(result.stdout.rstrip())
            if result.stderr.strip():
                print(result.stderr.rstrip(), file=sys.stderr)
            converted += 1

        except subprocess.CalledProcessError as exc:
            failed += 1
            print(f"  ERROR: converter exit code {exc.returncode}", file=sys.stderr)
            if exc.stdout:
                print(exc.stdout.rstrip())
            if exc.stderr:
                print(exc.stderr.rstrip(), file=sys.stderr)

            if args.stop_on_error:
                break

    print()
    print("=" * 80)
    print("CONVERSION SUMMARY")
    print("=" * 80)
    print(f"Found     : {len(t02_files)}")
    print(f"Processed : {converted}")
    print(f"Failed    : {failed}")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
