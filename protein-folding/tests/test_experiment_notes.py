from __future__ import annotations

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


notes = load_module("pf_experiment_notes", "scripts/utils/experiment_notes.py")


def test_export_experiment_notes_copies_curated_and_generates_fallback(tmp_path: Path) -> None:
    viewer_root = tmp_path / "viewer"
    workspace_root = tmp_path / "workspace"
    (workspace_root / "flu_trimer_ligand_experiments").mkdir(parents=True)
    (workspace_root / "flu_trimer_ligand_experiments" / "summary.md").write_text(
        "# Curated Flu Summary\n"
    )
    rows = [
        {
            "experiment": "Flu HA trimer ligand screen",
            "experiment_slug": "flu-ha-trimer-ligand-screen",
            "target_id": "a",
            "assembly": "trimer",
            "config": "cfg1",
            "sampling_steps": "40",
            "diffusion_samples": "1",
            "max_parallel_samples": "1",
            "recycling_steps": "3",
            "sample_rank": "0",
            "folder": "published/flu/a",
            "note": "",
        },
        {
            "experiment": "Custom Experiment",
            "experiment_slug": "custom-experiment",
            "target_id": "b",
            "assembly": "monomer",
            "config": "cfg2",
            "sampling_steps": "20",
            "diffusion_samples": "2",
            "max_parallel_samples": "1",
            "recycling_steps": "1",
            "sample_rank": "0",
            "folder": "published/custom/b",
            "note": "generated metadata",
        },
    ]

    notes.export_experiment_notes(rows, viewer_root, workspace_root)

    curated = viewer_root / "experiments" / "flu-ha-trimer-ligand-screen" / "summary.md"
    generated = viewer_root / "experiments" / "custom-experiment" / "summary.md"
    assert curated.read_text() == "# Curated Flu Summary\n"
    assert "This summary was generated automatically" in generated.read_text()
