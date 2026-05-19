from __future__ import annotations

# ruff: noqa: ANN401,D101,D103,E501,INP001,I001,PERF401,S603

import argparse
import csv
import json
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "flu_trimer_single_ligand_affinity"

A23_SMILES = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H]2O)[C@@H](CO)O[C@H]1O"
A26_SMILES = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H](O)[C@H](O)[C@H]2O)[C@@H](CO)O[C@H]1O"


@dataclass(frozen=True)
class Protein:
    key: str
    label: str
    group: str
    sequence: str


@dataclass(frozen=True)
class Ligand:
    key: str
    label: str
    field: str
    value: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run HA trimer + one affinity ligand screens.")
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--boltz-cmd", default="uv run boltz")
    parser.add_argument("--accelerator", default="gpu")
    parser.add_argument("--model", default="boltz2")
    parser.add_argument("--sampling-steps", type=int, default=80)
    parser.add_argument("--diffusion-samples", type=int, default=3)
    parser.add_argument("--max-parallel-samples", type=int, default=1)
    parser.add_argument("--recycling-steps", type=int, default=5)
    parser.add_argument("--sampling-steps-affinity", type=int, default=200)
    parser.add_argument("--diffusion-samples-affinity", type=int, default=5)
    parser.add_argument("--panel", choices=["full", "restricted"], default="full")
    parser.add_argument("--use-msa-server", action="store_true", default=True)
    parser.add_argument("--only-summary", action="store_true")
    parser.add_argument("--override", action="store_true")
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq).upper()))
            name = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    if name is not None:
        records.append((name, "".join(seq).upper()))
    return records


def selected_proteins() -> list[Protein]:
    rows: list[Protein] = []
    for header, sequence in read_fasta(ROOT / "examples" / "proteins_h1.fasta"):
        lower = header.lower()
        if "duck" in lower:
            rows.append(Protein("duck_bavaria_h1", header, "H1", sequence))
        elif "memphis" in lower:
            rows.append(Protein("memphis_h1", header, "H1", sequence))
    for header, sequence in read_fasta(ROOT / "examples" / "proteins_h5.fasta"):
        lower = header.lower()
        if "texas" in lower:
            rows.append(Protein("texas_h5", header, "H5", sequence))
        else:
            rows.append(Protein("california_h5", header, "H5", sequence))
    return rows


def h1_ligands() -> list[Ligand]:
    return [
        Ligand("a23", "alpha2-3 SIA glycan", "smiles", A23_SMILES),
        Ligand("a26", "alpha2-6 SIA glycan", "smiles", A26_SMILES),
    ]


def h5_ligands() -> list[Ligand]:
    return [
        Ligand("sia", "SIA", "ccd", "SIA"),
        Ligand("ngc", "NGC", "ccd", "NGC"),
    ]


def all_ligands() -> list[Ligand]:
    return [*h1_ligands(), *h5_ligands()]


def pairs(panel: str) -> list[tuple[Protein, Ligand]]:
    out = []
    for protein in selected_proteins():
        if panel == "restricted":
            ligands = h1_ligands() if protein.group == "H1" else h5_ligands()
        else:
            ligands = all_ligands()
        out.extend((protein, ligand) for ligand in ligands)
    return out


def job_id(protein: Protein, ligand: Ligand) -> str:
    return f"{protein.key}_{ligand.key}_ha_trimer_affinity"


def write_yaml(path: Path, protein: Protein, ligand: Ligand) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ligand_line = (
        f'      smiles: "{ligand.value}"'
        if ligand.field == "smiles"
        else f"      ccd: {ligand.value}"
    )
    path.write_text(
        "\n".join(
            [
                "version: 1",
                "sequences:",
                "  - protein:",
                "      id: [A, B, C]",
                f"      sequence: {protein.sequence}",
                "  - ligand:",
                "      id: D",
                ligand_line,
                "properties:",
                "  - affinity:",
                "      binder: D",
                "",
            ]
        )
    )


def config_name(args: argparse.Namespace) -> str:
    return (
        f"ha_trimer_1lig_aff_s{args.sampling_steps}_n{args.diffusion_samples}"
        f"_p{args.max_parallel_samples}_r{args.recycling_steps}"
        f"_as{args.sampling_steps_affinity}_an{args.diffusion_samples_affinity}"
    )


def command_for(args: argparse.Namespace, input_yaml: Path, run_out: Path) -> list[str]:
    cmd = [
        *shlex.split(args.boltz_cmd),
        "predict",
        str(input_yaml),
        "--model",
        args.model,
        "--accelerator",
        args.accelerator,
        "--diffusion_samples",
        str(args.diffusion_samples),
        "--max_parallel_samples",
        str(args.max_parallel_samples),
        "--sampling_steps",
        str(args.sampling_steps),
        "--recycling_steps",
        str(args.recycling_steps),
        "--sampling_steps_affinity",
        str(args.sampling_steps_affinity),
        "--diffusion_samples_affinity",
        str(args.diffusion_samples_affinity),
        "--out_dir",
        str(run_out),
    ]
    if args.use_msa_server:
        cmd.append("--use_msa_server")
    if args.override:
        cmd.append("--override")
    return cmd


def run_jobs(args: argparse.Namespace, work: list[tuple[Protein, Ligand]]) -> None:
    config = config_name(args)
    for protein, ligand in work:
        jid = job_id(protein, ligand)
        input_yaml = args.out_dir / "inputs" / f"{jid}.yaml"
        run_out = args.out_dir / "runs" / config / jid
        log_path = args.out_dir / "logs" / config / f"{jid}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        existing_affinity = first_existing_affinity(args.out_dir, jid)
        if existing_affinity is not None and not args.override:
            log_path.write_text(f"Skipping existing affinity prediction: {existing_affinity}\n")
            continue
        cmd = command_for(args, input_yaml, run_out)
        with log_path.open("w") as log_file:
            log_file.write("$ " + shlex.join(cmd) + "\n\n")
            result = subprocess.run(
                cmd,
                cwd=ROOT,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=False,
            )
            log_file.write(f"\nexit_code={result.returncode}\n")


def first_existing_affinity(out_dir: Path, jid: str) -> Path | None:
    matches = sorted(
        out_dir.glob(
            f"runs/ha_trimer_1lig_aff_*/{jid}/boltz_results_{jid}/predictions/{jid}/affinity_{jid}.json"
        )
    )
    return matches[0] if matches else None


def collect_rows(args: argparse.Namespace, work: list[tuple[Protein, Ligand]]) -> list[dict[str, Any]]:
    rows = []
    for config in configs_to_collect(args):
        values = config_values(config)
        for protein, ligand in work:
            jid = job_id(protein, ligand)
            pred_dir = (
                args.out_dir
                / "runs"
                / config
                / jid
                / f"boltz_results_{jid}"
                / "predictions"
                / jid
            )
            depth = msa_depth(args, jid, config)
            raw_depth = msa_raw_depth(args, jid, config)
            msa_csv = msa_path(args, jid, config)
            affinity_path = pred_dir / f"affinity_{jid}.json"
            affinity = json.loads(affinity_path.read_text()) if affinity_path.exists() else {}
            conf_paths = sorted(pred_dir.glob(f"confidence_{jid}_model_*.json"))
            if not conf_paths:
                status = log_status(args, jid, config)
                if status == "missing":
                    continue
                rows.append(
                    base_row(
                        protein,
                        ligand,
                        config,
                        values,
                        "",
                        status,
                        {},
                        affinity,
                        "",
                        "",
                        "",
                        depth,
                        raw_depth,
                        msa_csv,
                    )
                )
                continue
            for conf_path in conf_paths:
                rank = re.search(r"_model_(\d+)\.json$", conf_path.name).group(1)
                structure = pred_dir / f"{jid}_model_{rank}.cif"
                confidence = json.loads(conf_path.read_text())
                rows.append(
                    base_row(
                        protein,
                        ligand,
                        config,
                        values,
                        rank,
                        "ok" if structure.exists() else "no_structure",
                        confidence,
                        affinity,
                        rel_path(structure) if structure.exists() else "",
                        rel_path(conf_path),
                        rel_path(affinity_path) if affinity_path.exists() else "",
                        depth,
                        raw_depth,
                        msa_csv,
                    )
                )
    return rows


def configs_to_collect(args: argparse.Namespace) -> list[str]:
    runs_dir = args.out_dir / "runs"
    if not runs_dir.exists():
        return [config_name(args)]
    configs = sorted(
        path.name
        for path in runs_dir.iterdir()
        if path.is_dir() and path.name.startswith("ha_trimer_1lig_aff_")
    )
    return configs or [config_name(args)]


def config_values(config: str) -> dict[str, str]:
    patterns = {
        "sampling_steps": r"_s(\d+)",
        "diffusion_samples": r"_n(\d+)",
        "max_parallel_samples": r"_p(\d+)",
        "recycling_steps": r"_r(\d+)",
        "sampling_steps_affinity": r"_as(\d+)",
        "diffusion_samples_affinity": r"_an(\d+)",
    }
    return {
        key: (match.group(1) if (match := re.search(pattern, config)) else "")
        for key, pattern in patterns.items()
    }


def msa_path(args: argparse.Namespace, jid: str, config: str) -> Path:
    return args.out_dir / "runs" / config / jid / f"boltz_results_{jid}" / "msa" / f"{jid}_0.csv"


def processed_msa_path(args: argparse.Namespace, jid: str, config: str) -> Path:
    return args.out_dir / "runs" / config / jid / f"boltz_results_{jid}" / "processed" / "msa" / f"{jid}_0.npz"


def msa_depth(args: argparse.Namespace, jid: str, config: str) -> int | str:
    processed = processed_msa_path(args, jid, config)
    if processed.exists():
        with np.load(processed) as data:
            return int(data["sequences"].shape[0])
    raw_depth = msa_raw_depth(args, jid, config)
    if raw_depth == "":
        return ""
    return min(int(raw_depth), 8192)


def msa_raw_depth(args: argparse.Namespace, jid: str, config: str) -> int | str:
    path = msa_path(args, jid, config)
    if not path.exists():
        return ""
    with path.open() as handle:
        return max(0, sum(1 for _ in handle) - 1)


def base_row(
    protein: Protein,
    ligand: Ligand,
    config: str,
    values: dict[str, str],
    rank: str,
    status: str,
    confidence: dict[str, Any],
    affinity: dict[str, Any],
    structure_path: str,
    confidence_json: str,
    affinity_json: str,
    msa_depth_value: int | str,
    msa_raw_depth_value: int | str,
    msa_csv: Path,
) -> dict[str, Any]:
    return {
        "job_id": job_id(protein, ligand),
        "protein": protein.key,
        "protein_label": protein.label,
        "protein_group": protein.group,
        "ligand": ligand.key,
        "ligand_label": ligand.label,
        "ligand_field": ligand.field,
        "assembly": "HA trimer + 1 affinity ligand",
        "ha_chains": 3,
        "ligand_copies": 1,
        "monomer_length": len(protein.sequence),
        "total_polymer_residues": len(protein.sequence) * 3,
        "config": config,
        "sampling_steps": values["sampling_steps"],
        "diffusion_samples": values["diffusion_samples"],
        "max_parallel_samples": values["max_parallel_samples"],
        "recycling_steps": values["recycling_steps"],
        "sampling_steps_affinity": values["sampling_steps_affinity"],
        "diffusion_samples_affinity": values["diffusion_samples_affinity"],
        "sample_rank": rank,
        "status": status,
        "msa_depth": msa_depth_value,
        "msa_raw_depth": msa_raw_depth_value,
        "msa_csv": rel_path(msa_csv) if msa_csv.exists() else "",
        "confidence_score": confidence.get("confidence_score", ""),
        "complex_plddt": confidence.get("complex_plddt", ""),
        "ptm": confidence.get("ptm", ""),
        "iptm": confidence.get("iptm", ""),
        "ligand_iptm": confidence.get("ligand_iptm", ""),
        "protein_iptm": confidence.get("protein_iptm", ""),
        "complex_iplddt": confidence.get("complex_iplddt", ""),
        "complex_pde": confidence.get("complex_pde", ""),
        "complex_ipde": confidence.get("complex_ipde", ""),
        "affinity_pred_value": affinity.get("affinity_pred_value", ""),
        "affinity_probability_binary": affinity.get("affinity_probability_binary", ""),
        "affinity_pred_value1": affinity.get("affinity_pred_value1", ""),
        "affinity_probability_binary1": affinity.get("affinity_probability_binary1", ""),
        "affinity_pred_value2": affinity.get("affinity_pred_value2", ""),
        "affinity_probability_binary2": affinity.get("affinity_probability_binary2", ""),
        "structure_path": structure_path,
        "confidence_json": confidence_json,
        "affinity_json": affinity_json,
        "note": "",
    }


def log_status(args: argparse.Namespace, jid: str, config: str) -> str:
    log_path = args.out_dir / "logs" / config / f"{jid}.log"
    if not log_path.exists():
        return "missing"
    text = log_path.read_text(errors="replace").lower()
    if "ran out of memory" in text:
        return "oom"
    if "number of failed examples: 1" in text:
        return "failed"
    if "exit_code=0" in text:
        return "no_output"
    return "failed"


FIELDS = [
    "job_id",
    "protein",
    "protein_label",
    "protein_group",
    "ligand",
    "ligand_label",
    "ligand_field",
    "assembly",
    "ha_chains",
    "ligand_copies",
    "monomer_length",
    "total_polymer_residues",
    "config",
    "sampling_steps",
    "diffusion_samples",
    "max_parallel_samples",
    "recycling_steps",
    "sampling_steps_affinity",
    "diffusion_samples_affinity",
    "sample_rank",
    "status",
    "msa_depth",
    "msa_raw_depth",
    "msa_csv",
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
    "affinity_pred_value1",
    "affinity_probability_binary1",
    "affinity_pred_value2",
    "affinity_probability_binary2",
    "structure_path",
    "confidence_json",
    "affinity_json",
    "note",
]


def write_metadata(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    with (out_dir / "metadata.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def copy_structures(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    structures_dir = out_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        if row["status"] != "ok" or not row["structure_path"]:
            continue
        src = ROOT / row["structure_path"]
        dst_dir = structures_dir / row["config"]
        dst_dir.mkdir(parents=True, exist_ok=True)
        if src.exists():
            shutil.copy2(src, dst_dir / src.name)


def write_summary(args: argparse.Namespace, rows: list[dict[str, Any]]) -> None:
    best_rows = best_by_job(rows)
    lines = [
        "# Flu HA Trimer Single-Ligand Affinity Screen",
        "",
        "## Overview",
        "",
        f"- Conditions: {len({row['job_id'] for row in rows})}",
        f"- Completed conditions: {sum(1 for row in best_rows.values() if row['status'] == 'ok')}",
        f"- Completed affinity predictions: {sum(1 for row in best_rows.values() if row['affinity_json'])}",
        f"- Successful configs: {', '.join(f'`{config}`' for config in sorted({row['config'] for row in best_rows.values()}))}",
        "- HA modeled as a homotrimer: chains `A`, `B`, `C`.",
        "- Ligand modeled as a single affinity binder: chain `D`.",
        "",
        "## Summary",
        "",
        "| Protein | Ligand | Config | Status | Affinity | Binder prob | MSA depth | pLDDT | ipTM | ligand ipTM |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(best_rows.values(), key=lambda r: (r["protein_group"], r["protein"], r["ligand"])):
        lines.append(
            "| "
            f"{row['protein_label']} | {row['ligand_label']} | `{row['config']}` | {row['status']} | "
            f"{fmt(row['affinity_pred_value'])} | {fmt(row['affinity_probability_binary'])} | "
            f"{row['msa_depth'] or '-'} | {fmt(row['complex_plddt'])} | "
            f"{fmt(row['iptm'])} | {fmt(row['ligand_iptm'])} |"
        )
    lines += [
        "",
        "## Files",
        "",
        f"- Metadata CSV: `{rel_path(args.out_dir / 'metadata.csv')}`",
        f"- Structures: `{rel_path(args.out_dir / 'structures')}/`",
        f"- Raw outputs: `{rel_path(args.out_dir / 'runs')}/`",
        f"- Logs: `{rel_path(args.out_dir / 'logs')}/`",
    ]
    (args.out_dir / "summary.md").write_text("\n".join(lines) + "\n")


def best_by_job(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    best = {}
    for row in rows:
        if row["sample_rank"] not in {"", "0"}:
            continue
        current = best.get(row["job_id"])
        if current is None or row_score(row) > row_score(current):
            best[row["job_id"]] = row
    return best


def row_score(row: dict[str, Any]) -> tuple[int, int, float]:
    has_affinity = int(bool(row["affinity_json"]))
    is_ok = int(row["status"] == "ok")
    plddt = float(row["complex_plddt"] or 0)
    return has_affinity, is_ok, plddt


def fmt(value: Any) -> str:
    if value == "":
        return "-"
    return f"{float(value):.3f}"


def rel_path(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def main() -> None:
    args = parse_args()
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    work = pairs(args.panel)
    for protein, ligand in work:
        write_yaml(args.out_dir / "inputs" / f"{job_id(protein, ligand)}.yaml", protein, ligand)

    if not args.only_summary:
        run_jobs(args, work)

    rows = collect_rows(args, work)
    write_metadata(args.out_dir, rows)
    copy_structures(args.out_dir, rows)
    write_summary(args, rows)


if __name__ == "__main__":
    main()
