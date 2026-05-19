from __future__ import annotations

# ruff: noqa: D101,D102,D103,E501,INP001,I001,S603

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


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = ROOT / "test_EVE_structures.fasta"
DEFAULT_OUT = ROOT / "test_eve_experiments" / "multimers"
MONOMER_OUT = ROOT / "test_eve_experiments"


@dataclass(frozen=True)
class FastaRecord:
    index: int
    header: str
    sequence: str
    job_id: str


@dataclass(frozen=True)
class MultimerTarget:
    record: FastaRecord
    kind: str
    stoichiometry: int
    chain_ids: list[str]

    @property
    def job_id(self) -> str:
        return f"{self.record.job_id}_{self.kind}_{self.stoichiometry}mer"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run EVE homo-multimer Boltz trials.")
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--boltz-cmd", default="uv run boltz")
    parser.add_argument("--accelerator", default="gpu")
    parser.add_argument("--model", default="boltz2")
    parser.add_argument("--sampling-steps", type=int, default=40)
    parser.add_argument("--diffusion-samples", type=int, default=1)
    parser.add_argument("--max-parallel-samples", type=int, default=1)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--only-summary", action="store_true")
    parser.add_argument("--override", action="store_true")
    return parser.parse_args()


def parse_fasta(path: Path) -> list[FastaRecord]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq).upper()))
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq).upper()))

    return [
        FastaRecord(index, name, sequence, make_job_id(index, name))
        for index, (name, sequence) in enumerate(records, 1)
    ]


def make_job_id(index: int, header: str) -> str:
    text = header.split("[", 1)[0].strip()
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    text = text[:60].strip("_") or f"record_{index:02d}"
    return f"eve_{index:02d}_{text}"


def select_targets(records: list[FastaRecord]) -> list[MultimerTarget]:
    targets: list[MultimerTarget] = []
    for record in records:
        header = record.header.lower()
        if "gag" in header:
            targets.append(
                MultimerTarget(
                    record=record,
                    kind="gag",
                    stoichiometry=6,
                    chain_ids=list("ABCDEF"),
                )
            )
        if (
            " env " in f" {header} "
            or " env[" in header
            or "gpr80" in header
            or "gp60" in header
        ):
            targets.append(
                MultimerTarget(
                    record=record,
                    kind="env",
                    stoichiometry=3,
                    chain_ids=list("ABC"),
                )
            )
        if "sheplvp1" in header:
            targets.append(
                MultimerTarget(
                    record=record,
                    kind="sheplvp1",
                    stoichiometry=3,
                    chain_ids=list("ABC"),
                )
            )
    return targets


def baseline_msa_path(record: FastaRecord) -> Path:
    return (
        MONOMER_OUT
        / "runs"
        / "h5_start_s40_n6_p2_r3"
        / record.job_id
        / f"boltz_results_{record.job_id}"
        / "msa"
        / f"{record.job_id}_0.csv"
    )


def write_yaml(target: MultimerTarget, path: Path) -> None:
    msa_path = baseline_msa_path(target.record)
    ids = ", ".join(target.chain_ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "version: 1",
                "sequences:",
                "  - protein:",
                f"      id: [{ids}]",
                f"      sequence: {target.record.sequence}",
                f"      msa: {msa_path.resolve()}",
                "",
            ]
        )
    )


def command_for(args: argparse.Namespace, input_yaml: Path, out_dir: Path) -> list[str]:
    return [
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
        "--out_dir",
        str(out_dir),
        *(["--override"] if args.override else []),
    ]


def run_targets(args: argparse.Namespace, targets: list[MultimerTarget]) -> None:
    config = config_name(args)
    for target in targets:
        input_yaml = args.out_dir / "inputs" / f"{target.job_id}.yaml"
        run_out = args.out_dir / "runs" / config / target.job_id
        pred_dir = run_out / f"boltz_results_{target.job_id}" / "predictions" / target.job_id
        log_path = args.out_dir / "logs" / config / f"{target.job_id}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        if pred_dir.exists() and not args.override:
            log_path.write_text(f"Skipping existing prediction directory: {pred_dir}\n")
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


def collect_rows(args: argparse.Namespace, targets: list[MultimerTarget]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    config = config_name(args)
    for target in targets:
        pred_dir = (
            args.out_dir
            / "runs"
            / config
            / target.job_id
            / f"boltz_results_{target.job_id}"
            / "predictions"
            / target.job_id
        )
        conf_paths = sorted(pred_dir.glob(f"confidence_{target.job_id}_model_*.json"))
        if not conf_paths:
            rows.append(base_row(args, target, "", log_status(args, target), {}, "", ""))
            continue
        for conf_path in conf_paths:
            match = re.search(r"_model_(\d+)\.json$", conf_path.name)
            rank = match.group(1) if match else ""
            structure = pred_dir / f"{target.job_id}_model_{rank}.cif"
            confidence = json.loads(conf_path.read_text())
            rows.append(
                base_row(
                    args,
                    target,
                    rank,
                    "ok" if structure.exists() else "no_structure",
                    confidence,
                    rel_path(structure) if structure.exists() else "",
                    rel_path(conf_path),
                )
            )
    return rows


def base_row(
    args: argparse.Namespace,
    target: MultimerTarget,
    rank: str,
    status: str,
    confidence: dict[str, Any],
    structure_path: str,
    confidence_path: str,
) -> dict[str, Any]:
    return {
        "record_index": target.record.index,
        "job_id": target.job_id,
        "source_job_id": target.record.job_id,
        "fasta_header": target.record.header,
        "kind": target.kind,
        "stoichiometry": target.stoichiometry,
        "chain_count": len(target.chain_ids),
        "monomer_length": len(target.record.sequence),
        "total_residues": len(target.record.sequence) * len(target.chain_ids),
        "msa_raw_depth": msa_raw_depth(target.record),
        "config": config_name(args),
        "sampling_steps": args.sampling_steps,
        "diffusion_samples": args.diffusion_samples,
        "max_parallel_samples": args.max_parallel_samples,
        "recycling_steps": args.recycling_steps,
        "sample_rank": rank,
        "status": status,
        "confidence_score": confidence.get("confidence_score", ""),
        "complex_plddt": confidence.get("complex_plddt", ""),
        "ptm": confidence.get("ptm", ""),
        "iptm": confidence.get("iptm", ""),
        "complex_iplddt": confidence.get("complex_iplddt", ""),
        "complex_pde": confidence.get("complex_pde", ""),
        "complex_ipde": confidence.get("complex_ipde", ""),
        "structure_path": structure_path,
        "confidence_json": confidence_path,
    }


def log_status(args: argparse.Namespace, target: MultimerTarget) -> str:
    log_path = args.out_dir / "logs" / config_name(args) / f"{target.job_id}.log"
    if not log_path.exists():
        return "missing"
    text = log_path.read_text(errors="replace")
    if "ran out of memory" in text.lower():
        return "oom"
    if "Number of failed examples: 1" in text:
        return "failed"
    if "exit_code=0" in text:
        return "no_output"
    return "failed"


def msa_raw_depth(record: FastaRecord) -> int | str:
    path = baseline_msa_path(record)
    if not path.exists():
        return ""
    with path.open() as handle:
        return max(0, sum(1 for _ in handle) - 1)


def write_metadata(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "record_index",
        "job_id",
        "source_job_id",
        "fasta_header",
        "kind",
        "stoichiometry",
        "chain_count",
        "monomer_length",
        "total_residues",
        "msa_raw_depth",
        "config",
        "sampling_steps",
        "diffusion_samples",
        "max_parallel_samples",
        "recycling_steps",
        "sample_rank",
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
    with (out_dir / "metadata.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
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


def best_rows(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    best: dict[str, dict[str, Any]] = {}
    for row in rows:
        if row["status"] != "ok" or row["complex_plddt"] == "":
            continue
        current = best.get(row["job_id"])
        if current is None or float(row["complex_plddt"]) > float(current["complex_plddt"]):
            best[row["job_id"]] = row
    return best


def write_summary(args: argparse.Namespace, targets: list[MultimerTarget], rows: list[dict[str, Any]]) -> None:
    best = best_rows(rows)
    ok_targets = sum(1 for target in targets if target.job_id in best)
    lines = [
        "# Test EVE Multimer Experiment",
        "",
        "## Overview",
        "",
        f"- Targets: {len(targets)}",
        f"- Completed targets: {ok_targets}",
        f"- Config: `{config_name(args)}`",
        f"- Sampling steps: {args.sampling_steps}",
        f"- Diffusion samples: {args.diffusion_samples}",
        f"- Max parallel samples: {args.max_parallel_samples}",
        f"- Recycling steps: {args.recycling_steps}",
        "",
        "## Best Result Per Multimer",
        "",
        "| Target | Kind | Stoichiometry | Total residues | MSA raw | Status | Best rank | complex_plddt | ptm | confidence_score |",
        "|---|---|---:|---:|---:|---|---:|---:|---:|---:|",
    ]
    for target in targets:
        row = best.get(target.job_id)
        if row is None:
            status_rows = [row for row in rows if row["job_id"] == target.job_id]
            status = status_rows[0]["status"] if status_rows else "missing"
            lines.append(
                f"| {md_escape(target.record.header)} | {target.kind} | {target.stoichiometry} | {len(target.record.sequence) * target.stoichiometry} | {msa_raw_depth(target.record) or '-'} | {status} | - | - | - | - |"
            )
            continue
        lines.append(
            "| "
            f"{md_escape(target.record.header)} | {target.kind} | {target.stoichiometry} | "
            f"{row['total_residues']} | {row['msa_raw_depth'] or '-'} | {row['status']} | "
            f"{row['sample_rank']} | {float(row['complex_plddt']):.3f} | "
            f"{float(row['ptm']):.3f} | {float(row['confidence_score']):.3f} |"
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


def config_name(args: argparse.Namespace) -> str:
    return (
        f"multimer_s{args.sampling_steps}_n{args.diffusion_samples}"
        f"_p{args.max_parallel_samples}_r{args.recycling_steps}"
    )


def rel_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def md_escape(value: str) -> str:
    return value.replace("|", "\\|")


def write_repro(args: argparse.Namespace) -> None:
    script = ROOT / "experiments" / "run_test_eve_multimers.py"
    cmd = [
        "uv",
        "run",
        "python",
        str(script.relative_to(ROOT)),
        "--fasta",
        str(args.fasta),
        "--out-dir",
        str(args.out_dir),
        "--boltz-cmd",
        args.boltz_cmd,
        "--accelerator",
        args.accelerator,
        "--model",
        args.model,
        "--sampling-steps",
        str(args.sampling_steps),
        "--diffusion-samples",
        str(args.diffusion_samples),
        "--max-parallel-samples",
        str(args.max_parallel_samples),
        "--recycling-steps",
        str(args.recycling_steps),
    ]
    if args.override:
        cmd.append("--override")
    path = args.out_dir / "repro.sh"
    path.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n\n"
        f"cd {shlex.quote(str(ROOT))}\n"
        f"{shlex.join(cmd)}\n"
    )
    path.chmod(0o755)


def main() -> None:
    args = parse_args()
    args.fasta = args.fasta.resolve()
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    records = parse_fasta(args.fasta)
    targets = select_targets(records)
    for target in targets:
        write_yaml(target, args.out_dir / "inputs" / f"{target.job_id}.yaml")

    if not args.only_summary:
        run_targets(args, targets)

    rows = collect_rows(args, targets)
    write_metadata(args.out_dir, rows)
    copy_structures(args.out_dir, rows)
    write_summary(args, targets, rows)
    write_repro(args)


if __name__ == "__main__":
    main()
