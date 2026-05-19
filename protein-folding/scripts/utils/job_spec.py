from __future__ import annotations

import hashlib
import json
import re
from copy import deepcopy
from pathlib import Path
from typing import Any

IGNORED_HASH_FIELDS = {
    "job_hash",
    "job_id",
    "created_at",
    "updated_at",
    "queue_path",
    "queue_status",
    "run_dir",
    "log_path",
    "attempt",
    "retry_count",
    "last_error",
}


def load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        msg = f"Expected JSON object in {path}"
        raise ValueError(msg)
    return data


def dump_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def slugify(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return normalized or "item"


def yaml_digest(yaml_path: Path) -> str:
    return hashlib.sha256(yaml_path.read_bytes()).hexdigest()


def _strip_ignored(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _strip_ignored(item)
            for key, item in sorted(value.items())
            if key not in IGNORED_HASH_FIELDS
        }
    if isinstance(value, list):
        return [_strip_ignored(item) for item in value]
    if isinstance(value, Path):
        return value.as_posix()
    return value


def canonical_job_payload(job: dict[str, Any]) -> dict[str, Any]:
    payload = deepcopy(job)
    yaml_path_value = payload.get("yaml_path")
    if yaml_path_value:
        payload["yaml_path"] = Path(yaml_path_value).as_posix()
    return _strip_ignored(payload)


def compute_job_hash(job: dict[str, Any]) -> str:
    canonical = canonical_job_payload(job)
    encoded = json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()
    digest = hashlib.sha256(encoded).hexdigest()
    return f"sha256:{digest}"


def ensure_job_hash(job: dict[str, Any]) -> dict[str, Any]:
    payload = deepcopy(job)
    yaml_path_value = payload.get("yaml_path")
    if not yaml_path_value:
        msg = "Job spec is missing yaml_path"
        raise ValueError(msg)
    yaml_path = Path(yaml_path_value)
    payload["yaml_digest"] = yaml_digest(yaml_path)
    payload["yaml_path"] = yaml_path.as_posix()
    payload["target_slug"] = payload.get("target_slug") or slugify(
        str(payload.get("target_id") or payload.get("target_label") or "target")
    )
    payload["experiment_slug"] = payload.get("experiment_slug") or slugify(
        str(payload.get("experiment_name") or "experiment")
    )
    payload["job_hash"] = compute_job_hash(payload)
    return payload


def validate_job_hash(job: dict[str, Any]) -> None:
    expected = ensure_job_hash(job)["job_hash"]
    actual = job.get("job_hash")
    if actual != expected:
        msg = f"Job hash mismatch: expected {expected}, found {actual}"
        raise ValueError(msg)
