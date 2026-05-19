from __future__ import annotations

# ruff: noqa: D103,INP001

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "test_eve_experiments"


FIELDS = [
    "record_index",
    "job_id",
    "fasta_header",
    "assembly",
    "kind",
    "length",
    "msa_depth",
    "msa_raw_depth",
    "condition",
    "sampling_steps",
    "diffusion_samples",
    "max_parallel_samples",
    "recycling_steps",
    "model_rank",
    "status",
    "confidence_score",
    "complex_plddt",
    "ptm",
    "iptm",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
    "structure_path",
    "confidence_json",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def monomer_rows() -> list[dict[str, str]]:
    rows = []
    for row in read_csv(OUT / "metadata.csv"):
        if row["status"] != "ok":
            continue
        normalized = {field: "" for field in FIELDS}
        normalized.update(
            {
                "record_index": row["record_index"],
                "job_id": row["job_id"],
                "fasta_header": row["fasta_header"],
                "assembly": "monomer",
                "kind": "monomer",
                "length": row["length"],
                "msa_depth": row["msa_depth"],
                "msa_raw_depth": row["msa_raw_depth"],
                "condition": row["condition"],
                "sampling_steps": row["sampling_steps"],
                "diffusion_samples": row["diffusion_samples"],
                "max_parallel_samples": row["max_parallel_samples"],
                "recycling_steps": row["recycling_steps"],
                "model_rank": row["model_rank"],
                "status": row["status"],
                "confidence_score": row["confidence_score"],
                "complex_plddt": row["complex_plddt"],
                "ptm": row["ptm"],
                "iptm": row["iptm"],
                "complex_iplddt": row["complex_iplddt"],
                "complex_pde": row["complex_pde"],
                "complex_ipde": row["complex_ipde"],
                "structure_path": row["structure_path"],
                "confidence_json": row["confidence_json"],
            }
        )
        rows.append(normalized)
    return rows


def multimer_rows() -> list[dict[str, str]]:
    rows = []
    metadata = OUT / "multimers" / "metadata.csv"
    if not metadata.exists():
        return rows
    for row in read_csv(metadata):
        if row["status"] != "ok":
            continue
        label = f"{row['fasta_header']} ({row['stoichiometry']}-mer)"
        normalized = {field: "" for field in FIELDS}
        normalized.update(
            {
                "record_index": row["record_index"],
                "job_id": row["job_id"],
                "fasta_header": label,
                "assembly": f"{row['stoichiometry']}-mer",
                "kind": row["kind"],
                "length": row["total_residues"],
                "msa_depth": "",
                "msa_raw_depth": row["msa_raw_depth"],
                "condition": row["config"],
                "sampling_steps": row["sampling_steps"],
                "diffusion_samples": row["diffusion_samples"],
                "max_parallel_samples": row["max_parallel_samples"],
                "recycling_steps": row["recycling_steps"],
                "model_rank": row["sample_rank"],
                "status": row["status"],
                "confidence_score": row["confidence_score"],
                "complex_plddt": row["complex_plddt"],
                "ptm": row["ptm"],
                "iptm": row["iptm"],
                "complex_iplddt": row["complex_iplddt"],
                "complex_pde": row["complex_pde"],
                "complex_ipde": row["complex_ipde"],
                "structure_path": row["structure_path"],
                "confidence_json": row["confidence_json"],
            }
        )
        rows.append(normalized)
    return rows


def main() -> None:
    rows = monomer_rows() + multimer_rows()
    rows.sort(
        key=lambda row: (
            int(row["record_index"]),
            0 if row["assembly"] == "monomer" else 1,
            row["assembly"],
            row["condition"],
            int(row["model_rank"] or 0),
        )
    )
    with (OUT / "combined_metadata.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
