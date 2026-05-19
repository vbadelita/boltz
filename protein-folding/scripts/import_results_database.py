#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import sys
from collections import defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.common import PUBLISHED_ROOT, VIEWER_ROOT
from utils.job_spec import slugify
from utils.viewer_db import base_row, normalize_existing_row, read_rows, write_rows

WORKSPACE_ROOT = VIEWER_ROOT.parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Import rows from an existing shared results database into viewer/published."
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=WORKSPACE_ROOT / "results_viewer" / "database.csv",
        help="Existing shared database CSV to import from.",
    )
    parser.add_argument(
        "--published-root",
        type=Path,
        default=PUBLISHED_ROOT,
        help="Destination published root inside the standalone viewer.",
    )
    parser.add_argument(
        "--experiment",
        action="append",
        default=[],
        help="Only import matching experiment names. Repeatable.",
    )
    return parser.parse_args()


def group_key(row: dict[str, str]) -> tuple[str, str]:
    experiment_slug = row["experiment_slug"] or slugify(row["experiment"])
    folder = row.get("folder", "").strip()
    if folder:
        return experiment_slug, slugify(folder)
    structure_path = Path(row.get("structure_path", ""))
    return experiment_slug, slugify(structure_path.parts[0] if structure_path.parts else "bundle")


def resolve_source_path(path_value: str) -> Path | None:
    if not path_value:
        return None
    path = Path(path_value)
    if path.is_absolute():
        return path
    candidate = WORKSPACE_ROOT / path
    return candidate if candidate.exists() else None


def relative_artifact_path(row: dict[str, str], path_value: str) -> Path | None:
    if not path_value:
        return None
    path = Path(path_value)
    folder = row.get("folder", "").strip()
    if folder:
        folder_path = Path(folder)
        if path.parts[: len(folder_path.parts)] == folder_path.parts:
            return Path(*path.parts[len(folder_path.parts) :])
    return path


def import_rows(
    database_path: Path,
    published_root: Path,
    experiments: set[str],
) -> int:
    rows = [normalize_existing_row(row) for row in read_rows(database_path)]
    if experiments:
      rows = [row for row in rows if row.get("experiment", "") in experiments]

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[group_key(row)].append(row)

    imported = 0
    for (experiment_slug, bundle_slug), group_rows in sorted(grouped.items()):
        bundle_dir = published_root / experiment_slug / bundle_slug
        if bundle_dir.exists():
            shutil.rmtree(bundle_dir)
        bundle_dir.mkdir(parents=True, exist_ok=True)

        rewritten_rows = []
        copied_paths: set[Path] = set()
        for row in group_rows:
            new_row = base_row(**row)
            new_row["experiment_slug"] = experiment_slug
            new_row["target_slug"] = row["target_slug"] or slugify(
                row["target_id"] or row["target_label"]
            )
            new_row["folder"] = bundle_dir.relative_to(VIEWER_ROOT).as_posix()

            for field in ("structure_path", "confidence_json", "affinity_json"):
                source_path = resolve_source_path(row.get(field, ""))
                relative_path = relative_artifact_path(row, row.get(field, ""))
                if source_path is None or relative_path is None:
                    continue
                destination = bundle_dir / relative_path
                if source_path.exists() and destination not in copied_paths:
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source_path, destination)
                    copied_paths.add(destination)
                new_row[field] = destination.relative_to(VIEWER_ROOT).as_posix()
            rewritten_rows.append(new_row)

        write_rows(bundle_dir / "viewer.csv", rewritten_rows)
        imported += 1
    return imported


def main() -> int:
    args = parse_args()
    count = import_rows(
        args.database.resolve(),
        args.published_root.resolve(),
        set(args.experiment),
    )
    print(f"Imported {count} bundle(s) from {args.database}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
