from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare operational and publication CSV outputs using exact pandas DataFrame equality.")
    p.add_argument("--old-dir", required=True, type=Path)
    p.add_argument("--new-dir", required=True, type=Path)
    p.add_argument("--files", nargs="+", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    all_ok = True
    for name in args.files:
        old_path = args.old_dir / name
        new_path = args.new_dir / name
        print("\n" + name)
        if not old_path.is_file() or not new_path.is_file():
            print("  MISSING:", "old" if not old_path.is_file() else "", "new" if not new_path.is_file() else "")
            all_ok = False
            continue
        old = pd.read_csv(old_path)
        new = pd.read_csv(new_path)
        same_columns = list(old.columns) == list(new.columns)
        exact = new.equals(old)
        print("  new rows:", len(new))
        print("  old rows:", len(old))
        print("  same columns:", same_columns)
        print("  EXACT EQUALITY:", exact)
        all_ok &= same_columns and exact
    print("\nOVERALL EXACT EQUALITY:", all_ok)
    raise SystemExit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
