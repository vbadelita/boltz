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


job_spec = load_module("pf_job_spec", "scripts/utils/job_spec.py")


def test_job_hash_is_deterministic_for_equivalent_json(tmp_path: Path) -> None:
    yaml_path = tmp_path / "target.yaml"
    yaml_path.write_text("version: 1\nsequences:\n  - protein:\n      sequence: MPEPTIDE\n")

    job_a = {
        "experiment_name": "Example experiment",
        "experiment_slug": "example-experiment",
        "target_id": "target_001",
        "target_label": "Target 001",
        "yaml_path": str(yaml_path),
        "program": "boltz",
        "args": {"sampling_steps": 40, "diffusion_samples": 1},
        "env": {"BOLTZ_SRC": "boltz"},
    }
    job_b = {
        "target_label": "Target 001",
        "target_id": "target_001",
        "experiment_slug": "example-experiment",
        "program": "boltz",
        "yaml_path": str(yaml_path),
        "env": {"BOLTZ_SRC": "boltz"},
        "experiment_name": "Example experiment",
        "args": {"diffusion_samples": 1, "sampling_steps": 40},
        "created_at": "ignored",
    }

    normalized_a = job_spec.ensure_job_hash(job_a)
    normalized_b = job_spec.ensure_job_hash(job_b)

    assert normalized_a["job_hash"] == normalized_b["job_hash"]


def test_job_hash_changes_when_yaml_changes(tmp_path: Path) -> None:
    yaml_a = tmp_path / "a.yaml"
    yaml_b = tmp_path / "b.yaml"
    yaml_a.write_text("sequence: AAAA\n")
    yaml_b.write_text("sequence: BBBB\n")

    job_a = job_spec.ensure_job_hash(
        {
            "experiment_name": "Example experiment",
            "target_id": "target_001",
            "target_label": "Target 001",
            "yaml_path": str(yaml_a),
            "program": "boltz",
            "args": {},
        }
    )
    job_b = job_spec.ensure_job_hash(
        {
            "experiment_name": "Example experiment",
            "target_id": "target_001",
            "target_label": "Target 001",
            "yaml_path": str(yaml_b),
            "program": "boltz",
            "args": {},
        }
    )

    assert job_a["job_hash"] != job_b["job_hash"]
