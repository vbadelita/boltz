from __future__ import annotations

# ruff: noqa: D103,E501,INP001,I001

import csv
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results_viewer"
DATABASE = OUT_DIR / "database.csv"


FIELDS = [
    "experiment",
    "target_id",
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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def empty_row() -> dict[str, str]:
    return dict.fromkeys(FIELDS, "")


def eve_rows() -> list[dict[str, str]]:
    path = ROOT / "test_eve_experiments" / "combined_metadata.csv"
    if not path.exists():
        return []
    rows = []
    for source in read_csv(path):
        if source["status"] != "ok":
            continue
        row = empty_row()
        row.update(
            {
                "experiment": "EVE structures",
                "target_id": source["job_id"],
                "target_label": source["fasta_header"],
                "assembly": source["assembly"],
                "kind": source["kind"],
                "length": source["length"],
                "msa_depth": source["msa_depth"],
                "msa_raw_depth": source["msa_raw_depth"],
                "config": source["condition"],
                "sampling_steps": source["sampling_steps"],
                "diffusion_samples": source["diffusion_samples"],
                "max_parallel_samples": source["max_parallel_samples"],
                "recycling_steps": source["recycling_steps"],
                "sample_rank": source["model_rank"],
                "status": source["status"],
                "confidence_score": source["confidence_score"],
                "complex_plddt": source["complex_plddt"],
                "ptm": source["ptm"],
                "iptm": source["iptm"],
                "complex_iplddt": source["complex_iplddt"],
                "complex_pde": source["complex_pde"],
                "complex_ipde": source["complex_ipde"],
                "structure_path": source["structure_path"],
                "confidence_json": source["confidence_json"],
            }
        )
        rows.append(row)
    return rows


def influenza_experiment(run_name: str, target: str) -> str | None:
    if run_name.startswith("h5_"):
        return "Influenza H5 affinity matrix"
    if "flu_sia_affinity_msa" in run_name or "no_msa" in run_name:
        return "Influenza MSA ladder"
    if "a23" in run_name or "a26" in run_name:
        return "Influenza glycan linkage"
    if run_name.startswith(("flu_", "boltz_results_flu")):
        return "Influenza SIA/NGC ladder"
    if target.startswith("flu"):
        return "Influenza other"
    return None


def config_values(name: str) -> dict[str, str]:
    values = {
        "sampling_steps": "",
        "diffusion_samples": "",
        "max_parallel_samples": "",
        "recycling_steps": "",
    }
    match = re.search(r"_s(\d+)", name)
    if match:
        values["sampling_steps"] = match.group(1)
    match = re.search(r"_n(\d+)", name)
    if match:
        values["diffusion_samples"] = match.group(1)
    match = re.search(r"_p(\d+)", name)
    if match:
        values["max_parallel_samples"] = match.group(1)
    match = re.search(r"_r(\d+)", name)
    if match:
        values["recycling_steps"] = match.group(1)
    return values


def affinity_for(pred_dir: Path, target: str) -> tuple[str, dict[str, Any]]:
    affinity = pred_dir / f"affinity_{target}.json"
    if affinity.exists():
        return rel_path(affinity), load_json(affinity)
    matches = sorted(pred_dir.glob("affinity_*.json"))
    if matches:
        return rel_path(matches[0]), load_json(matches[0])
    return "", {}


def influenza_rows() -> list[dict[str, str]]:
    rows = []
    for conf_path in sorted((ROOT / "runs").glob("*/boltz_results_*/predictions/*/confidence_*_model_*.json")):
        run_dir = conf_path.relative_to(ROOT).parts[1]
        pred_dir = conf_path.parent
        target = pred_dir.name
        experiment = influenza_experiment(run_dir, target)
        if experiment is None:
            continue
        rank_match = re.search(r"_model_(\d+)\.json$", conf_path.name)
        rank = rank_match.group(1) if rank_match else "0"
        structure = pred_dir / f"{target}_model_{rank}.cif"
        if not structure.exists():
            continue
        confidence = load_json(conf_path)
        affinity_path, affinity = affinity_for(pred_dir, target)
        row = empty_row()
        row.update(
            {
                "experiment": experiment,
                "target_id": f"{run_dir}/{target}",
                "target_label": target,
                "assembly": "complex" if affinity else "structure",
                "kind": "influenza",
                "config": run_dir,
                "sample_rank": rank,
                "status": "ok",
                "confidence_score": str(confidence.get("confidence_score", "")),
                "complex_plddt": str(confidence.get("complex_plddt", "")),
                "ptm": str(confidence.get("ptm", "")),
                "iptm": str(confidence.get("iptm", "")),
                "ligand_iptm": str(confidence.get("ligand_iptm", "")),
                "protein_iptm": str(confidence.get("protein_iptm", "")),
                "complex_iplddt": str(confidence.get("complex_iplddt", "")),
                "complex_pde": str(confidence.get("complex_pde", "")),
                "complex_ipde": str(confidence.get("complex_ipde", "")),
                "affinity_pred_value": str(affinity.get("affinity_pred_value", "")),
                "affinity_probability_binary": str(
                    affinity.get("affinity_probability_binary", "")
                ),
                "structure_path": rel_path(structure),
                "confidence_json": rel_path(conf_path),
                "affinity_json": affinity_path,
            }
        )
        row.update(config_values(run_dir))
        rows.append(row)
    return rows


def flu_trimer_rows() -> list[dict[str, str]]:
    path = ROOT / "flu_trimer_ligand_experiments" / "metadata.csv"
    if not path.exists():
        return []

    rows = []
    for source in read_csv(path):
        if source["status"] != "ok":
            continue
        row = empty_row()
        row.update(
            {
                "experiment": "Flu HA trimer ligand screen",
                "target_id": source["job_id"],
                "target_label": f"{source['protein_label']} + {source['ligand_label']}",
                "assembly": source["assembly"],
                "kind": source["protein_group"],
                "length": source["total_polymer_residues"],
                "msa_depth": source["msa_depth"],
                "msa_raw_depth": source["msa_raw_depth"],
                "msa_csv": source["msa_csv"],
                "config": source["config"],
                "sampling_steps": source["sampling_steps"],
                "diffusion_samples": source["diffusion_samples"],
                "max_parallel_samples": source["max_parallel_samples"],
                "recycling_steps": source["recycling_steps"],
                "sample_rank": source["sample_rank"],
                "status": source["status"],
                "confidence_score": source["confidence_score"],
                "complex_plddt": source["complex_plddt"],
                "ptm": source["ptm"],
                "iptm": source["iptm"],
                "ligand_iptm": source["ligand_iptm"],
                "protein_iptm": source["protein_iptm"],
                "complex_iplddt": source["complex_iplddt"],
                "complex_pde": source["complex_pde"],
                "complex_ipde": source["complex_ipde"],
                "structure_path": source["structure_path"],
                "confidence_json": source["confidence_json"],
                "note": source.get("note", ""),
            }
        )
        rows.append(row)
    return rows


def flu_single_ligand_affinity_rows() -> list[dict[str, str]]:
    path = ROOT / "flu_trimer_single_ligand_affinity" / "metadata.csv"
    if not path.exists():
        return []

    rows = []
    for source in read_csv(path):
        if source["status"] != "ok":
            continue
        row = empty_row()
        row.update(
            {
                "experiment": "Flu HA trimer single-ligand affinity",
                "target_id": source["job_id"],
                "target_label": f"{source['protein_label']} + {source['ligand_label']}",
                "assembly": source["assembly"],
                "kind": source["protein_group"],
                "length": source["total_polymer_residues"],
                "msa_depth": source["msa_depth"],
                "msa_raw_depth": source["msa_raw_depth"],
                "msa_csv": source["msa_csv"],
                "config": source["config"],
                "sampling_steps": source["sampling_steps"],
                "diffusion_samples": source["diffusion_samples"],
                "max_parallel_samples": source["max_parallel_samples"],
                "recycling_steps": source["recycling_steps"],
                "sampling_steps_affinity": source["sampling_steps_affinity"],
                "diffusion_samples_affinity": source["diffusion_samples_affinity"],
                "sample_rank": source["sample_rank"],
                "status": source["status"],
                "confidence_score": source["confidence_score"],
                "complex_plddt": source["complex_plddt"],
                "ptm": source["ptm"],
                "iptm": source["iptm"],
                "ligand_iptm": source["ligand_iptm"],
                "protein_iptm": source["protein_iptm"],
                "complex_iplddt": source["complex_iplddt"],
                "complex_pde": source["complex_pde"],
                "complex_ipde": source["complex_ipde"],
                "affinity_pred_value": source["affinity_pred_value"],
                "affinity_probability_binary": source["affinity_probability_binary"],
                "structure_path": source["structure_path"],
                "confidence_json": source["confidence_json"],
                "affinity_json": source["affinity_json"],
                "note": source.get("note", ""),
            }
        )
        rows.append(row)
    return rows


def rel_path(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = (
        eve_rows()
        + influenza_rows()
        + flu_trimer_rows()
        + flu_single_ligand_affinity_rows()
    )
    rows.sort(
        key=lambda row: (
            row["experiment"],
            row["target_label"],
            row["assembly"],
            row["config"],
            int(row["sample_rank"] or 0),
        )
    )
    with DATABASE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
