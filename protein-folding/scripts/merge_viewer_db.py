#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.common import PUBLISHED_ROOT, VIEWER_DB_PATH, VIEWER_ROOT
from utils.experiment_notes import export_experiment_notes
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


def resolve_bundle_artifact(bundle_dir: Path, source_path: str) -> Path:
    path = Path(source_path)
    if not source_path:
        return bundle_dir
    if path.exists():
        return path
    parts = path.parts
    for index in range(len(parts)):
        candidate = bundle_dir / Path(*parts[index:])
        if candidate.exists():
            return candidate
    return bundle_dir / path.name


def relative_to_viewer_root(path: Path, viewer_root: Path) -> str:
    return path.relative_to(viewer_root).as_posix()


def main() -> int:
    args = parse_args()
    viewer_root = args.output.resolve().parent
    merged: dict[tuple[str, str, str], dict[str, str]] = {}
    for path in discover_published_viewers(args.published_root.resolve()):
        bundle_dir = path.parent
        for source in read_rows(path):
            row = normalize_existing_row(source)
            row["folder"] = relative_to_viewer_root(bundle_dir, viewer_root)
            structure_path = resolve_bundle_artifact(bundle_dir, row.get("structure_path", ""))
            confidence_path = resolve_bundle_artifact(bundle_dir, row.get("confidence_json", ""))
            affinity_path = resolve_bundle_artifact(bundle_dir, row.get("affinity_json", ""))
            row["structure_path"] = relative_to_viewer_root(structure_path, viewer_root)
            row["confidence_json"] = relative_to_viewer_root(confidence_path, viewer_root)
            if row.get("affinity_json"):
                row["affinity_json"] = relative_to_viewer_root(affinity_path, viewer_root)
            merged[row_key(row)] = row
    rows = [merged[key] for key in sorted(merged)]
    write_rows(args.output.resolve(), rows)
    export_experiment_notes(rows, viewer_root, VIEWER_ROOT.parent.parent)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
