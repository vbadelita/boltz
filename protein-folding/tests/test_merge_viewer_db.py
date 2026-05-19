from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path


def load_module(name: str, relative_path: str):
    repo_root = Path(__file__).resolve().parents[1]
    scripts_root = repo_root / "scripts"
    if str(scripts_root) not in sys.path:
        sys.path.insert(0, str(scripts_root))
    module_path = repo_root / relative_path
    spec = importlib.util.spec_from_file_location(name, module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


viewer_db = load_module("pf_viewer_db", "scripts/utils/viewer_db.py")


def write_viewer_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=viewer_db.DATABASE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def test_merge_viewer_db_deduplicates_by_experiment_target_and_sample(tmp_path: Path) -> None:
    viewer_root = tmp_path / "viewer"
    published_root = viewer_root / "published"
    output_path = viewer_root / "database.csv"
    row_a = viewer_db.base_row(
        experiment="Example",
        experiment_slug="example",
        folder="published/example/run-1",
        target_id="target_001",
        target_slug="target-001",
        target_label="Target 001",
        sample_rank="0",
        status="ok",
        structure_path="published/example/run-1/predictions/target_001_model_0.cif",
        confidence_json="published/example/run-1/predictions/confidence_target_001_model_0.json",
    )
    row_b = dict(row_a)
    row_b["complex_plddt"] = "0.91"

    write_viewer_csv(published_root / "example" / "run-1" / "viewer.csv", [row_a])
    write_viewer_csv(published_root / "example" / "run-2" / "viewer.csv", [row_b])

    subprocess.run(  # noqa: S603
        [
            sys.executable,
            str(Path(__file__).resolve().parents[1] / "scripts" / "merge_viewer_db.py"),
            "--published-root",
            str(published_root),
            "--output",
            str(output_path),
        ],
        check=True,
    )

    with output_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    assert rows[0]["complex_plddt"] == "0.91"
