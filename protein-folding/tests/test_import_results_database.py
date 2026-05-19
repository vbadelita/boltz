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


import_db = load_module("pf_import_db", "scripts/import_results_database.py")
viewer_db = load_module("pf_viewer_db_import_db", "scripts/utils/viewer_db.py")


def test_import_results_database_copies_rows_into_viewer_bundle(tmp_path: Path, monkeypatch) -> None:
    workspace_root = tmp_path
    viewer_root = workspace_root / "protein-folding" / "viewer"
    published_root = viewer_root / "published"
    source_dir = workspace_root / "flu_trimer_ligand_experiments" / "runs" / "run_a" / "target_a"
    source_dir.mkdir(parents=True)
    structure = source_dir / "target_a_model_0.cif"
    confidence = source_dir / "confidence_target_a_model_0.json"
    structure.write_text("data_test\n")
    confidence.write_text("{}\n")

    db_path = workspace_root / "results_viewer" / "database.csv"
    db_path.parent.mkdir(parents=True)
    row = viewer_db.base_row(
        experiment="Flu HA trimer ligand screen",
        folder="flu_trimer_ligand_experiments",
        target_id="target_a",
        target_slug="target-a",
        target_label="Target A",
        sample_rank="0",
        status="ok",
        structure_path=structure.relative_to(workspace_root).as_posix(),
        confidence_json=confidence.relative_to(workspace_root).as_posix(),
    )
    with db_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=viewer_db.DATABASE_COLUMNS)
        writer.writeheader()
        writer.writerow(row)

    monkeypatch.setattr(import_db, "WORKSPACE_ROOT", workspace_root)
    monkeypatch.setattr(import_db, "VIEWER_ROOT", viewer_root)
    count = import_db.import_rows(db_path, published_root, {"Flu HA trimer ligand screen"})

    assert count == 1
    imported_csv = published_root / "flu-ha-trimer-ligand-screen" / "flu-trimer-ligand-experiments" / "viewer.csv"
    rows = viewer_db.read_rows(imported_csv)
    assert rows[0]["folder"] == "published/flu-ha-trimer-ligand-screen/flu-trimer-ligand-experiments"
    assert rows[0]["structure_path"].startswith(
        "published/flu-ha-trimer-ligand-screen/flu-trimer-ligand-experiments/"
    )
