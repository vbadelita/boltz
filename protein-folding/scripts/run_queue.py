#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.common import HOSTS_CONFIG_PATH, QUEUE_DIRS, REPO_ROOT, ensure_repo_layout
from utils.job_spec import dump_json, ensure_job_hash, load_json, validate_job_hash
from utils.viewer_db import build_viewer_rows, write_rows

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]


def iso_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one or more protein-folding queue jobs.")
    parser.add_argument("--queue-root", type=Path, default=QUEUE_DIRS["pending"].parent)
    parser.add_argument("--once", action="store_true", help="Claim and run at most one job.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the command that would run without executing it.",
    )
    parser.add_argument(
        "--allow-failed-rerun",
        action="store_true",
        help="Ignore existing failed markers for the same job hash.",
    )
    return parser.parse_args()


def load_hosts_config() -> dict[str, Any]:
    if not HOSTS_CONFIG_PATH.exists():
        return {}
    return tomllib.loads(HOSTS_CONFIG_PATH.read_text())


def hash_index_path(job_hash: str, queue_root: Path) -> Path:
    digest = job_hash.split(":", 1)[1]
    prefix = digest[:2]
    return queue_root / "by-hash" / prefix / f"{digest}.json"


def read_hash_status(index_path: Path) -> dict[str, Any] | None:
    if not index_path.exists():
        return None
    return load_json(index_path)


def write_hash_status(index_path: Path, payload: dict[str, Any]) -> None:
    dump_json(index_path, payload)


def claim_job(pending_path: Path, running_dir: Path) -> Path | None:
    destination = running_dir / pending_path.name
    try:
        pending_path.replace(destination)
    except FileNotFoundError:
        return None
    return destination


def discover_predictions(run_dir: Path, target_id: str) -> list[dict[str, Any]]:
    predictions_dir = run_dir / "predictions"
    rows = []
    if not predictions_dir.exists():
        return rows
    for structure_path in sorted(predictions_dir.glob(f"{target_id}_model_*.cif")):
        sample_suffix = structure_path.stem.rsplit("_model_", 1)[1]
        confidence_path = predictions_dir / f"confidence_{target_id}_model_{sample_suffix}.json"
        metrics: dict[str, Any] = {}
        if confidence_path.exists():
            metrics = json.loads(confidence_path.read_text())
        rows.append(
            {
                "sample_rank": sample_suffix,
                "status": "ok",
                "structure_path": structure_path,
                "confidence_json": confidence_path,
                "confidence_score": metrics.get("confidence_score", ""),
                "complex_plddt": metrics.get("complex_plddt", ""),
                "ptm": metrics.get("ptm", ""),
                "iptm": metrics.get("iptm", ""),
                "ligand_iptm": metrics.get("ligand_iptm", ""),
                "protein_iptm": metrics.get("protein_iptm", ""),
                "complex_iplddt": metrics.get("complex_iplddt", ""),
                "complex_pde": metrics.get("complex_pde", ""),
                "complex_ipde": metrics.get("complex_ipde", ""),
            }
        )
    return rows


def program_command(
    program: str,
    job: dict[str, Any],
    run_dir: Path,
    hosts_config: dict[str, Any],
) -> tuple[list[str], dict[str, str]]:
    args = {str(key): value for key, value in dict(job.get("args", {})).items()}
    env = os.environ.copy()
    env.update({str(key): str(value) for key, value in dict(job.get("env", {})).items()})
    yaml_path = Path(str(job["yaml_path"])).resolve()
    python_bin = env.get("PYTHON_BIN", sys.executable)
    if program == "boltz":
        command = [python_bin, "-m", "boltz.main", "predict", str(yaml_path)]
    elif program == "lmi4boltz":
        lmi_src = env.get("LMI4BOLTZ_SRC") or str(
            Path(hosts_config.get("compute", {}).get("lmi4boltz_src", "lmi4boltz")).resolve()
        )
        existing_pythonpath = env.get("PYTHONPATH", "")
        src_path = str(Path(lmi_src) / "src")
        env["PYTHONPATH"] = (
            f"{src_path}{os.pathsep}{existing_pythonpath}" if existing_pythonpath else src_path
        )
        command = [python_bin, "-m", "boltz.main", "predict", str(yaml_path)]
    else:
        msg = f"Unsupported program: {program}"
        raise ValueError(msg)

    flag_map = {
        "model": "--model",
        "accelerator": "--accelerator",
        "sampling_steps": "--sampling_steps",
        "diffusion_samples": "--diffusion_samples",
        "max_parallel_samples": "--max_parallel_samples",
        "recycling_steps": "--recycling_steps",
        "chunk_size_transition_z": "--chunk_size_transition_z",
        "chunk_size_transition_msa": "--chunk_size_transition_msa",
        "chunk_size_outer_product": "--chunk_size_outer_product",
        "chunk_size_tri_attn": "--chunk_size_tri_attn",
        "triangle_mult_gate_nchunks": "--triangle_mult_gate_nchunks",
        "chunk_size_threshold": "--chunk_size_threshold",
    }
    for key, value in args.items():
        if key == "use_msa_server":
            if value:
                command.append("--use_msa_server")
            continue
        if key == "override":
            if value:
                command.append("--override")
            continue
        if key == "out_dir":
            continue
        flag = flag_map.get(key)
        if flag is None:
            continue
        command.extend([flag, str(value)])
    command.extend(["--out_dir", str(run_dir)])
    return command, env


def write_run_json(run_dir: Path, job: dict[str, Any], command: list[str]) -> None:
    payload = {
        "job_hash": job["job_hash"],
        "job_id": job.get("job_id", ""),
        "program": job["program"],
        "command": command,
        "args": job.get("args", {}),
        "experiment_slug": job["experiment_slug"],
        "experiment_name": job["experiment_name"],
        "target_id": job["target_id"],
        "target_label": job["target_label"],
        "assembly": job.get("assembly", ""),
        "kind": job.get("kind", ""),
        "length": job.get("length", ""),
        "started_at": iso_now(),
    }
    dump_json(run_dir / "run.json", payload)


def update_run_json(run_dir: Path, **updates: Any) -> None:
    path = run_dir / "run.json"
    payload = load_json(path)
    payload.update(updates)
    dump_json(path, payload)


def run_job(
    job_path: Path,
    queue_root: Path,
    hosts_config: dict[str, Any],
    dry_run: bool,
) -> tuple[str, Path]:
    job = ensure_job_hash(load_json(job_path))
    validate_job_hash(job)
    experiment_slug = str(job["experiment_slug"])
    target_slug = str(job["target_slug"])
    run_id = f"{target_slug}-{job['job_hash'].split(':', 1)[1][:12]}"
    run_dir = REPO_ROOT / "experiments" / experiment_slug / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(job["yaml_path"], run_dir / Path(str(job["yaml_path"])).name)

    command, env = program_command(str(job["program"]), job, run_dir, hosts_config)
    write_run_json(run_dir, job, command)
    log_path = run_dir / "logs" / "runner.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    if dry_run:
        update_run_json(
            run_dir,
            finished_at=iso_now(),
            status="dry_run",
            note="Command not executed.",
        )
        return "dry_run", run_dir

    with log_path.open("w") as handle:
        result = subprocess.run(  # noqa: S603
            command,
            cwd=REPO_ROOT,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    if result.returncode != 0:
        update_run_json(
            run_dir,
            finished_at=iso_now(),
            status="failed",
            exit_code=result.returncode,
            log_path=log_path.relative_to(REPO_ROOT).as_posix(),
        )
        return "failed", run_dir

    predictions = discover_predictions(run_dir, str(job["target_id"]))
    if not predictions:
        update_run_json(
            run_dir,
            finished_at=iso_now(),
            status="failed",
            exit_code=0,
            note="No predictions discovered after command completed.",
            log_path=log_path.relative_to(REPO_ROOT).as_posix(),
        )
        return "failed", run_dir

    viewer_rows = build_viewer_rows(run_dir, job, predictions)
    write_rows(run_dir / "viewer.csv", viewer_rows)
    update_run_json(
        run_dir,
        finished_at=iso_now(),
        status="ok",
        exit_code=0,
        log_path=log_path.relative_to(REPO_ROOT).as_posix(),
    )
    return "ok", run_dir


def maybe_skip_job(
    job_path: Path,
    queue_root: Path,
    allow_failed_rerun: bool,
) -> tuple[bool, str]:
    job = ensure_job_hash(load_json(job_path))
    index_path = hash_index_path(str(job["job_hash"]), queue_root)
    status = read_hash_status(index_path)
    if not status:
        return False, ""
    if status.get("status") == "succeeded":
        return True, "already_succeeded"
    if status.get("status") == "running":
        return True, "already_running"
    if status.get("status") == "failed" and not allow_failed_rerun:
        return True, "already_failed"
    return False, ""


def write_skip_marker(queue_root: Path, job: dict[str, Any], status: str, run_dir: str = "") -> None:
    index_path = hash_index_path(str(job["job_hash"]), queue_root)
    write_hash_status(
        index_path,
        {
            "job_hash": job["job_hash"],
            "status": status,
            "experiment_slug": job["experiment_slug"],
            "target_id": job["target_id"],
            "run_dir": run_dir,
            "updated_at": iso_now(),
        },
    )


def process_pending_job(
    pending_path: Path,
    queue_root: Path,
    hosts_config: dict[str, Any],
    dry_run: bool,
    allow_failed_rerun: bool,
) -> bool:
    skip, reason = maybe_skip_job(pending_path, queue_root, allow_failed_rerun)
    job = ensure_job_hash(load_json(pending_path))
    if skip:
        destination = queue_root / "succeeded" / pending_path.name
        if reason == "already_failed":
            destination = queue_root / "failed" / pending_path.name
        pending_path.replace(destination)
        return True

    claimed = claim_job(pending_path, queue_root / "running")
    if claimed is None:
        return False

    write_skip_marker(queue_root, job, "running")
    status, run_dir = run_job(claimed, queue_root, hosts_config, dry_run)
    if status == "dry_run":
        destination = queue_root / "pending" / claimed.name
        claimed.replace(destination)
        write_skip_marker(queue_root, job, "pending", run_dir.relative_to(REPO_ROOT).as_posix())
        return True

    if status == "ok":
        destination = queue_root / "succeeded" / claimed.name
        claimed.replace(destination)
        write_skip_marker(queue_root, job, "succeeded", run_dir.relative_to(REPO_ROOT).as_posix())
        return True

    destination = queue_root / "failed" / claimed.name
    claimed.replace(destination)
    write_skip_marker(queue_root, job, "failed", run_dir.relative_to(REPO_ROOT).as_posix())
    return True


def main() -> int:
    args = parse_args()
    queue_root = args.queue_root.resolve()
    ensure_repo_layout()
    for name in ("pending", "running", "succeeded", "failed", "by-hash"):
        (queue_root / name).mkdir(parents=True, exist_ok=True)

    hosts_config = load_hosts_config()
    processed = 0
    for pending_path in sorted((queue_root / "pending").glob("*.json")):
        if process_pending_job(
            pending_path,
            queue_root,
            hosts_config,
            args.dry_run,
            args.allow_failed_rerun,
        ):
            processed += 1
            if args.once:
                break
    print(f"Processed {processed} job(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
