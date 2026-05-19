from __future__ import annotations

import importlib.util
import json
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


run_queue = load_module("pf_run_queue", "scripts/run_queue.py")
job_spec = load_module("pf_job_spec_run_queue", "scripts/utils/job_spec.py")


def test_process_pending_job_skips_when_hash_already_succeeded(tmp_path: Path) -> None:
    queue_root = tmp_path / "queue"
    for name in ("pending", "running", "succeeded", "failed", "by-hash"):
        (queue_root / name).mkdir(parents=True, exist_ok=True)

    yaml_path = tmp_path / "target.yaml"
    yaml_path.write_text("sequence: AAAA\n")
    job = job_spec.ensure_job_hash(
        {
            "experiment_name": "Example experiment",
            "target_id": "target_001",
            "target_label": "Target 001",
            "yaml_path": str(yaml_path),
            "program": "boltz",
            "args": {},
        }
    )
    pending_path = queue_root / "pending" / "job.json"
    pending_path.write_text(json.dumps(job))

    index_path = run_queue.hash_index_path(job["job_hash"], queue_root)
    run_queue.write_hash_status(index_path, {"job_hash": job["job_hash"], "status": "succeeded"})

    processed = run_queue.process_pending_job(
        pending_path,
        queue_root,
        hosts_config={},
        dry_run=True,
        allow_failed_rerun=False,
    )

    assert processed is True
    assert not pending_path.exists()
    assert (queue_root / "succeeded" / "job.json").exists()


def test_process_pending_job_claims_and_marks_success(tmp_path: Path, monkeypatch) -> None:
    queue_root = tmp_path / "queue"
    for name in ("pending", "running", "succeeded", "failed", "by-hash"):
        (queue_root / name).mkdir(parents=True, exist_ok=True)

    yaml_path = tmp_path / "target.yaml"
    yaml_path.write_text("sequence: AAAA\n")
    job = job_spec.ensure_job_hash(
        {
            "experiment_name": "Example experiment",
            "target_id": "target_001",
            "target_label": "Target 001",
            "yaml_path": str(yaml_path),
            "program": "boltz",
            "args": {},
        }
    )
    pending_path = queue_root / "pending" / "job.json"
    pending_path.write_text(json.dumps(job))

    monkeypatch.setattr(run_queue, "REPO_ROOT", tmp_path)
    fake_run_dir = tmp_path / "experiments" / "example-experiment" / "runs" / "run-001"
    fake_run_dir.mkdir(parents=True, exist_ok=True)

    def fake_run_job(job_path, queue_root, hosts_config, dry_run):
        return "ok", fake_run_dir

    monkeypatch.setattr(run_queue, "run_job", fake_run_job)

    processed = run_queue.process_pending_job(
        pending_path,
        queue_root,
        hosts_config={},
        dry_run=False,
        allow_failed_rerun=False,
    )

    assert processed is True
    assert (queue_root / "succeeded" / "job.json").exists()
    status = run_queue.read_hash_status(run_queue.hash_index_path(job["job_hash"], queue_root))
    assert status is not None
    assert status["status"] == "succeeded"
