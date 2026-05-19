from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable

try:
    from .common import REPO_ROOT, repo_relative
except ImportError:  # pragma: no cover
    from utils.common import REPO_ROOT, repo_relative  # type: ignore[no-redef]

DATABASE_COLUMNS = [
    "experiment",
    "experiment_slug",
    "folder",
    "target_id",
    "target_slug",
    "target_label",
    "assembly",
    "kind",
    "length",
    "msa_depth",
    "msa_raw_depth",
    "msa_csv",
    "config",
    "sampling_steps",
    "diffusion_samples",
    "max_parallel_samples",
    "recycling_steps",
    "sampling_steps_affinity",
    "diffusion_samples_affinity",
    "sample_rank",
    "status",
    "confidence_score",
    "complex_plddt",
    "ptm",
    "iptm",
    "ligand_iptm",
    "protein_iptm",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
    "affinity_pred_value",
    "affinity_probability_binary",
    "structure_path",
    "confidence_json",
    "affinity_json",
    "note",
]


def base_row(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in DATABASE_COLUMNS}
    row.update({key: str(value) for key, value in overrides.items() if key in row})
    return row


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DATABASE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def row_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("experiment_slug", ""),
        row.get("target_id", ""),
        row.get("sample_rank", ""),
    )


def normalize_existing_row(source: dict[str, str]) -> dict[str, str]:
    row = base_row(**source)
    if not row["experiment_slug"] and row["experiment"]:
        try:
            from .job_spec import slugify
        except ImportError:  # pragma: no cover
            from utils.job_spec import slugify  # type: ignore[no-redef]

        row["experiment_slug"] = slugify(row["experiment"])
    if not row["target_slug"] and row["target_id"]:
        try:
            from .job_spec import slugify
        except ImportError:  # pragma: no cover
            from utils.job_spec import slugify  # type: ignore[no-redef]

        row["target_slug"] = slugify(row["target_id"])
    return row


def build_viewer_rows(
    run_dir: Path,
    job: dict[str, object],
    predictions: list[dict[str, object]],
) -> list[dict[str, str]]:
    run_json_path = run_dir / "run.json"
    run_data = json.loads(run_json_path.read_text()) if run_json_path.exists() else {}
    rows = []
    for prediction in predictions:
        structure_path = Path(str(prediction["structure_path"]))
        confidence_path = Path(str(prediction["confidence_json"]))
        row = base_row(
            experiment=str(job["experiment_name"]),
            experiment_slug=str(job["experiment_slug"]),
            folder=repo_relative(run_dir),
            target_id=str(job["target_id"]),
            target_slug=str(job["target_slug"]),
            target_label=str(job["target_label"]),
            assembly=str(job.get("assembly", "")),
            kind=str(job.get("kind", "")),
            length=str(job.get("length", "")),
            config=str(run_data.get("program", job.get("program", ""))),
            sampling_steps=str(run_data.get("args", {}).get("sampling_steps", "")),
            diffusion_samples=str(run_data.get("args", {}).get("diffusion_samples", "")),
            max_parallel_samples=str(
                run_data.get("args", {}).get("max_parallel_samples", "")
            ),
            recycling_steps=str(run_data.get("args", {}).get("recycling_steps", "")),
            sample_rank=str(prediction.get("sample_rank", "")),
            status=str(prediction.get("status", "ok")),
            confidence_score=str(prediction.get("confidence_score", "")),
            complex_plddt=str(prediction.get("complex_plddt", "")),
            ptm=str(prediction.get("ptm", "")),
            iptm=str(prediction.get("iptm", "")),
            ligand_iptm=str(prediction.get("ligand_iptm", "")),
            protein_iptm=str(prediction.get("protein_iptm", "")),
            complex_iplddt=str(prediction.get("complex_iplddt", "")),
            complex_pde=str(prediction.get("complex_pde", "")),
            complex_ipde=str(prediction.get("complex_ipde", "")),
            structure_path=repo_relative(structure_path),
            confidence_json=repo_relative(confidence_path),
            note=str(prediction.get("note", "")),
        )
        rows.append(row)
    return rows


def discover_published_viewers(published_root: Path) -> list[Path]:
    return sorted(published_root.glob("*/*/viewer.csv"))
