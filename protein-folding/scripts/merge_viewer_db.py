#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.common import PUBLISHED_ROOT, VIEWER_DB_PATH
from utils.viewer_db import (
    discover_published_viewers,
    normalize_existing_row,
    read_rows,
    row_key,
    write_rows,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge published viewer CSV files.")
    parser.add_argument("--published-root", type=Path, default=PUBLISHED_ROOT)
    parser.add_argument("--output", type=Path, default=VIEWER_DB_PATH)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    merged: dict[tuple[str, str, str], dict[str, str]] = {}
    for path in discover_published_viewers(args.published_root.resolve()):
        for source in read_rows(path):
            row = normalize_existing_row(source)
            merged[row_key(row)] = row
    rows = [merged[key] for key in sorted(merged)]
    write_rows(args.output.resolve(), rows)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
