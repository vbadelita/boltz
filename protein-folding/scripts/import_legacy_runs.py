#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.common import PUBLISHED_ROOT, REPO_ROOT, repo_relative, viewer_relative
from utils.job_spec import slugify
from utils.viewer_db import base_row, normalize_existing_row, read_rows, write_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Import legacy Boltz run-local viewer bundles into protein-folding/viewer/published."
    )
    parser.add_argument(
        "--source",
        action="append",
        default=[],
        help="Legacy run directory containing viewer.csv. Repeatable.",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=REPO_ROOT.parent / "runs",
        help="Directory to scan for legacy run folders.",
    )
    parser.add_argument(
        "--published-root",
        type=Path,
        default=PUBLISHED_ROOT,
        help="Destination published root inside the viewer web root.",
    )
    parser.add_argument(
        "--copy-inputs",
        action="store_true",
        help="Copy the legacy run inputs/ directory when present.",
    )
    return parser.parse_args()


def legacy_run_dirs(args: argparse.Namespace) -> list[Path]:
    if args.source:
        return [Path(item).resolve() for item in args.source]
    return sorted(path.parent for path in args.source_root.resolve().glob("*/viewer.csv"))


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def relative_to_run(path_value: str, run_dir: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path.relative_to(run_dir)
    parts = path.parts
    if parts and parts[0] == "runs":
        path = Path(*parts[2:]) if len(parts) > 2 else Path()
    return path


def resolve_source_artifact(path_value: str, run_dir: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    if path.parts and path.parts[0] == "runs":
        workspace_root = run_dir.parent.parent
        return workspace_root / path
    return run_dir / path


def import_run(run_dir: Path, published_root: Path, copy_inputs: bool) -> tuple[str, Path, int]:
    viewer_path = run_dir / "viewer.csv"
    rows = [normalize_existing_row(row) for row in read_rows(viewer_path)]
    if not rows:
        return "empty", run_dir, 0

    experiment_slug = rows[0]["experiment_slug"] or slugify(rows[0]["experiment"])
    bundle_dir = published_root / experiment_slug / run_dir.name
    if bundle_dir.exists():
        shutil.rmtree(bundle_dir)
    bundle_dir.mkdir(parents=True, exist_ok=True)

    rewritten_rows = []
    copied_paths: set[Path] = set()
    for row in rows:
        structure_source = resolve_source_artifact(row["structure_path"], run_dir)
        confidence_source = resolve_source_artifact(row["confidence_json"], run_dir)
        structure_rel = relative_to_run(row["structure_path"], run_dir)
        confidence_rel = relative_to_run(row["confidence_json"], run_dir)
        structure_dest = bundle_dir / structure_rel
        confidence_dest = bundle_dir / confidence_rel
        if structure_source.exists() and structure_dest not in copied_paths:
            copy_file(structure_source, structure_dest)
            copied_paths.add(structure_dest)
        if confidence_source.exists() and confidence_dest not in copied_paths:
            copy_file(confidence_source, confidence_dest)
            copied_paths.add(confidence_dest)

        new_row = base_row(**row)
        new_row["experiment_slug"] = experiment_slug
        new_row["target_slug"] = row["target_slug"] or slugify(row["target_id"] or row["target_label"])
        new_row["folder"] = viewer_relative(bundle_dir)
        new_row["structure_path"] = viewer_relative(structure_dest)
        new_row["confidence_json"] = viewer_relative(confidence_dest)
        rewritten_rows.append(new_row)

    if copy_inputs and (run_dir / "inputs").exists():
        shutil.copytree(run_dir / "inputs", bundle_dir / "inputs", dirs_exist_ok=True)

    write_rows(bundle_dir / "viewer.csv", rewritten_rows)
    return experiment_slug, bundle_dir, len(rewritten_rows)


def main() -> int:
    args = parse_args()
    imported = 0
    for run_dir in legacy_run_dirs(args):
        if not (run_dir / "viewer.csv").exists():
            continue
        experiment_slug, bundle_dir, row_count = import_run(
            run_dir.resolve(), args.published_root.resolve(), args.copy_inputs
        )
        print(
            f"Imported {run_dir.as_posix()} -> {bundle_dir.relative_to(REPO_ROOT)} "
            f"({row_count} rows, {experiment_slug})"
        )
        imported += 1
    print(f"Imported {imported} legacy run bundle(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
