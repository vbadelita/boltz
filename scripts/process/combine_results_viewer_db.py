#!/usr/bin/env python3
"""Combine multiple experiment metadata files into results_viewer/database.csv."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_PATH = ROOT / "results_viewer" / "database.csv"
DATABASE_COLUMNS = [
    "experiment",
    "folder",
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


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def top_level_folder(path: str) -> str:
    return Path(path).parts[0] if path else ""


def normalize_folder_label(folder: str) -> str:
    if not folder:
        return ""
    folder_path = Path(folder)
    if folder_path.is_absolute():
        try:
            return folder_path.relative_to(ROOT).as_posix()
        except ValueError:
            return folder_path.as_posix()
    return folder_path.as_posix()


def base_row(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in DATABASE_COLUMNS}
    row.update({key: value for key, value in overrides.items() if key in row})
    return row


def convert_existing_results(
    path: Path, sourced_experiments: set[str] | None = None
) -> list[dict[str, str]]:
    rows = []
    explicitly_sourced = [
        "flu_trimer_ligand_experiments",
        "flu_trimer_single_ligand_affinity",
        "test_eve_experiments",
    ]
    for source in read_rows(path):
        structure_root = top_level_folder(source.get("structure_path", ""))
        folder = normalize_folder_label(source.get("folder") or structure_root)
        if not folder:
            continue
        if structure_root in explicitly_sourced:
            continue
        if any(
            folder == explicit or folder.startswith(f"{explicit}/")
            for explicit in explicitly_sourced
        ):
            continue
        if sourced_experiments and source.get("experiment", "") in sourced_experiments:
            continue
        row = base_row(**source)
        row["folder"] = structure_root or folder
        rows.append(row)
    return rows


def collect_run_viewer_paths() -> list[Path]:
    """Return run-local viewer CSVs that should be folded into the shared DB."""

    return sorted((ROOT / "runs").glob("*/viewer.csv"))


def convert_run_viewer(path: Path) -> list[dict[str, str]]:
    """Normalize rows from a run-local viewer CSV."""

    rows = []
    for source in read_rows(path):
        row = base_row(**source)
        if not row["folder"]:
            structure_root = top_level_folder(row.get("structure_path", ""))
            row["folder"] = structure_root
        rows.append(row)
    return rows


def convert_flu_metadata(path: Path, experiment: str) -> list[dict[str, str]]:
    folder = path.parent.relative_to(ROOT).as_posix()
    rows = []
    for source in read_rows(path):
        row = base_row(**source)
        row["experiment"] = experiment
        row["folder"] = folder
        row["target_id"] = source["job_id"]
        row["target_label"] = source["job_id"]
        row["kind"] = source.get("protein_group") or source.get("assembly", "")
        row["length"] = source.get("total_polymer_residues") or source.get(
            "monomer_length", ""
        )
        rows.append(row)
    return rows


def convert_eve_combined(path: Path) -> list[dict[str, str]]:
    folder = path.parent.relative_to(ROOT).as_posix()
    rows = []
    for source in read_rows(path):
        row = base_row(
            experiment="EVE structures",
            folder=folder,
            target_id=source["job_id"],
            target_label=source["fasta_header"],
            assembly=source.get("assembly", "monomer"),
            kind=source.get("kind", "monomer"),
            length=source.get("length", ""),
            msa_depth=source.get("msa_depth", ""),
            msa_raw_depth=source.get("msa_raw_depth", ""),
            config=source.get("condition", ""),
            sampling_steps=source.get("sampling_steps", ""),
            diffusion_samples=source.get("diffusion_samples", ""),
            max_parallel_samples=source.get("max_parallel_samples", ""),
            recycling_steps=source.get("recycling_steps", ""),
            sample_rank=source.get("model_rank", ""),
            status=source.get("status", ""),
            confidence_score=source.get("confidence_score", ""),
            complex_plddt=source.get("complex_plddt", ""),
            ptm=source.get("ptm", ""),
            iptm=source.get("iptm", ""),
            complex_iplddt=source.get("complex_iplddt", ""),
            complex_pde=source.get("complex_pde", ""),
            complex_ipde=source.get("complex_ipde", ""),
            structure_path=source.get("structure_path", ""),
            confidence_json=source.get("confidence_json", ""),
        )
        rows.append(row)
    return rows


def convert_eve_multimers(path: Path) -> list[dict[str, str]]:
    folder = path.parent.relative_to(ROOT).as_posix()
    rows = []
    for source in read_rows(path):
        row = base_row(
            experiment="EVE multimers",
            folder=folder,
            target_id=source["job_id"],
            target_label=source["fasta_header"],
            assembly=f'{source.get("stoichiometry", "")}-mer',
            kind=source.get("kind", ""),
            length=source.get("total_residues", ""),
            msa_depth=source.get("msa_raw_depth", ""),
            msa_raw_depth=source.get("msa_raw_depth", ""),
            config=source.get("config", ""),
            sampling_steps=source.get("sampling_steps", ""),
            diffusion_samples=source.get("diffusion_samples", ""),
            max_parallel_samples=source.get("max_parallel_samples", ""),
            recycling_steps=source.get("recycling_steps", ""),
            sample_rank=source.get("sample_rank", ""),
            status=source.get("status", ""),
            confidence_score=source.get("confidence_score", ""),
            complex_plddt=source.get("complex_plddt", ""),
            ptm=source.get("ptm", ""),
            iptm=source.get("iptm", ""),
            complex_iplddt=source.get("complex_iplddt", ""),
            complex_pde=source.get("complex_pde", ""),
            complex_ipde=source.get("complex_ipde", ""),
            structure_path=source.get("structure_path", ""),
            confidence_json=source.get("confidence_json", ""),
            note=(
                f'source_job_id={source.get("source_job_id", "")}; '
                f'chains={source.get("chain_count", "")}; '
                f'monomer_length={source.get("monomer_length", "")}'
            ),
        )
        rows.append(row)
    return rows


def row_sort_key(row: dict[str, str]) -> tuple[str, str, str, int]:
    try:
        sample_rank = int(row.get("sample_rank", "0") or 0)
    except ValueError:
        sample_rank = 0
    return (
        row.get("folder", ""),
        row.get("experiment", ""),
        row.get("target_id", ""),
        sample_rank,
    )


def main() -> int:
    run_viewer_paths = collect_run_viewer_paths()
    sourced_experiments = {
        row.get("experiment", "")
        for path in run_viewer_paths
        for row in read_rows(path)
        if row.get("experiment", "")
    }
    rows: list[dict[str, str]] = []
    rows.extend(
        convert_existing_results(
            ROOT / "results_viewer" / "database.csv",
            sourced_experiments=sourced_experiments,
        )
    )
    rows.extend(
        convert_flu_metadata(
            ROOT / "flu_trimer_ligand_experiments" / "metadata.csv",
            "Flu HA trimer ligand screen",
        )
    )
    rows.extend(
        convert_flu_metadata(
            ROOT / "flu_trimer_single_ligand_affinity" / "metadata.csv",
            "Flu HA trimer single-ligand affinity",
        )
    )
    rows.extend(
        convert_eve_combined(ROOT / "test_eve_experiments" / "combined_metadata.csv")
    )
    rows.extend(
        convert_eve_multimers(
            ROOT / "test_eve_experiments" / "multimers" / "metadata.csv"
        )
    )
    for path in run_viewer_paths:
        rows.extend(convert_run_viewer(path))
    rows.sort(key=row_sort_key)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DATABASE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUTPUT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
