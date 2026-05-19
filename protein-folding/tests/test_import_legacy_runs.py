from __future__ import annotations

import csv
import importlib.util
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


import_legacy = load_module("pf_import_legacy", "scripts/import_legacy_runs.py")
viewer_db = load_module("pf_viewer_db_import_legacy", "scripts/utils/viewer_db.py")


def test_import_run_copies_viewer_artifacts_and_rewrites_paths(tmp_path: Path, monkeypatch) -> None:
    repo_root = tmp_path
    run_dir = repo_root / "runs" / "legacy_run"
    pred_dir = run_dir / "target_a" / "boltz_results_target_a" / "predictions" / "target_a"
    pred_dir.mkdir(parents=True)
    structure_source = pred_dir / "target_a_model_0.cif"
    confidence_source = pred_dir / "confidence_target_a_model_0.json"
    structure_source.write_text("data_test\n")
    confidence_source.write_text("{}\n")

    viewer_path = run_dir / "viewer.csv"
    viewer_path.parent.mkdir(parents=True, exist_ok=True)
    row = viewer_db.base_row(
        experiment="Legacy Example",
        folder="runs",
        target_id="target_a",
        target_label="Target A",
        sample_rank="0",
        status="ok",
        structure_path=structure_source.relative_to(repo_root).as_posix(),
        confidence_json=confidence_source.relative_to(repo_root).as_posix(),
    )
    with viewer_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=viewer_db.DATABASE_COLUMNS)
        writer.writeheader()
        writer.writerow(row)

    monkeypatch.setattr(import_legacy, "REPO_ROOT", repo_root)
    monkeypatch.setattr(
        import_legacy,
        "repo_relative",
        lambda path: path.relative_to(repo_root).as_posix(),
    )
    monkeypatch.setattr(
        import_legacy,
        "viewer_relative",
        lambda path: path.relative_to(repo_root / "protein-folding" / "viewer").as_posix(),
    )
    published_root = repo_root / "protein-folding" / "viewer" / "published"
    experiment_slug, bundle_dir, row_count = import_legacy.import_run(
        run_dir, published_root, copy_inputs=False
    )

    assert experiment_slug == "legacy-example"
    assert row_count == 1
    imported_csv = bundle_dir / "viewer.csv"
    assert imported_csv.exists()
    imported_rows = viewer_db.read_rows(imported_csv)
    assert imported_rows[0]["folder"] == bundle_dir.relative_to(repo_root / "protein-folding" / "viewer").as_posix()
    assert imported_rows[0]["structure_path"].startswith("published/legacy-example/legacy_run/")
    assert (repo_root / "protein-folding" / "viewer" / imported_rows[0]["structure_path"]).exists()
