from __future__ import annotations

# ruff: noqa: C901,D101,D103,E501,INP001,I001,PLR2004,S603

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
DEFAULT_FASTA = ROOT / "test_EVE_structures.fasta"
DEFAULT_OUT = ROOT / "test_eve_experiments"


@dataclass(frozen=True)
class FastaRecord:
    index: int
    header: str
    sequence: str
    job_id: str


@dataclass(frozen=True)
class Condition:
    name: str
    sampling_steps: int
    diffusion_samples: int
    max_parallel_samples: int
    recycling_steps: int
    use_msa_server: bool
    include_missing: bool = True


BASELINE = Condition(
    name="h5_start_s40_n6_p2_r3",
    sampling_steps=40,
    diffusion_samples=6,
    max_parallel_samples=2,
    recycling_steps=3,
    use_msa_server=True,
)

GRID_CONDITIONS = [
    Condition(
        name="more_steps_s80_n6_p2_r3",
        sampling_steps=80,
        diffusion_samples=6,
        max_parallel_samples=2,
        recycling_steps=3,
        use_msa_server=False,
    ),
    Condition(
        name="more_samples_s40_n8_p2_r3",
        sampling_steps=40,
        diffusion_samples=8,
        max_parallel_samples=2,
        recycling_steps=3,
        use_msa_server=False,
    ),
    Condition(
        name="more_recycling_s40_n6_p2_r5",
        sampling_steps=40,
        diffusion_samples=6,
        max_parallel_samples=2,
        recycling_steps=5,
        use_msa_server=False,
    ),
]

FAILED_RETRY_CONDITIONS = [
    Condition(
        name="retry_s200_n6_p2_r5",
        sampling_steps=200,
        diffusion_samples=6,
        max_parallel_samples=2,
        recycling_steps=5,
        use_msa_server=False,
        include_missing=False,
    )
]

ALL_200_CONDITIONS = [
    Condition(
        name="all_s200_n6_p2_r5",
        sampling_steps=200,
        diffusion_samples=6,
        max_parallel_samples=2,
        recycling_steps=5,
        use_msa_server=False,
        include_missing=False,
    )
]

HIGH_STEP_CONDITIONS = [*FAILED_RETRY_CONDITIONS, *ALL_200_CONDITIONS]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Boltz2 on test_EVE_structures.fasta and summarize results."
    )
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--boltz-cmd",
        default="uv run boltz",
        help="Command prefix used to invoke Boltz.",
    )
    parser.add_argument("--accelerator", default="gpu")
    parser.add_argument("--model", default="boltz2")
    parser.add_argument(
        "--low-plddt",
        type=float,
        default=0.70,
        help="Baseline complex_plddt below this value is targeted by the grid.",
    )
    parser.add_argument(
        "--grid-pilot-size",
        type=int,
        default=3,
        help="If fewer targets are below --low-plddt, grid-search this many lowest baseline targets.",
    )
    parser.add_argument(
        "--only-summary",
        action="store_true",
        help="Do not run Boltz; only rebuild metadata.csv, structures/, and summary.md.",
    )
    parser.add_argument(
        "--override",
        action="store_true",
        help="Pass --override to Boltz and recompute existing predictions.",
    )
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

    parsed = []
    for index, (name, sequence) in enumerate(records, 1):
        parsed.append(FastaRecord(index, name, sequence, make_job_id(index, name)))
    return parsed


def make_job_id(index: int, header: str) -> str:
    text = header.split("[", 1)[0].strip()
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    text = text[:60].strip("_") or f"record_{index:02d}"
    return f"eve_{index:02d}_{text}"


def write_yaml(record: FastaRecord, path: Path, msa_path: Path | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {record.sequence}",
    ]
    if msa_path is not None:
        lines.append(f"      msa: {msa_path}")
    path.write_text("\n".join(lines) + "\n")


def command_for(
    args: argparse.Namespace,
    condition: Condition,
    input_yaml: Path,
    out_dir: Path,
) -> list[str]:
    cmd = shlex.split(args.boltz_cmd)
    cmd += [
        "predict",
        str(input_yaml),
        "--model",
        args.model,
        "--accelerator",
        args.accelerator,
        "--diffusion_samples",
        str(condition.diffusion_samples),
        "--max_parallel_samples",
        str(condition.max_parallel_samples),
        "--sampling_steps",
        str(condition.sampling_steps),
        "--recycling_steps",
        str(condition.recycling_steps),
        "--out_dir",
        str(out_dir),
    ]
    if condition.use_msa_server:
        cmd.append("--use_msa_server")
    if args.override:
        cmd.append("--override")
    return cmd


def run_condition(
    args: argparse.Namespace,
    records: list[FastaRecord],
    condition: Condition,
    input_dir: Path,
    out_dir: Path,
    logs_dir: Path,
) -> None:
    for record in records:
        input_yaml = input_dir / f"{record.job_id}.yaml"
        run_out = out_dir / condition.name / record.job_id
        log_path = logs_dir / condition.name / f"{record.job_id}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        prediction_dir = (
            run_out
            / f"boltz_results_{record.job_id}"
            / "predictions"
            / record.job_id
        )
        if prediction_dir.exists() and not args.override:
            log_path.write_text(
                f"Skipping existing prediction directory: {prediction_dir}\n"
            )
            continue

        cmd = command_for(args, condition, input_yaml, run_out)
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


def baseline_msa_path(out_dir: Path, record: FastaRecord) -> Path:
    return (
        out_dir
        / "runs"
        / BASELINE.name
        / record.job_id
        / f"boltz_results_{record.job_id}"
        / "msa"
        / f"{record.job_id}_0.csv"
    )


def prediction_dir(out_dir: Path, condition: str, record: FastaRecord) -> Path:
    return (
        out_dir
        / "runs"
        / condition
        / record.job_id
        / f"boltz_results_{record.job_id}"
        / "predictions"
        / record.job_id
    )


def msa_depth(out_dir: Path, record: FastaRecord) -> int | str:
    processed = (
        out_dir
        / "runs"
        / BASELINE.name
        / record.job_id
        / f"boltz_results_{record.job_id}"
        / "processed"
        / "msa"
        / f"{record.job_id}_0.npz"
    )
    if processed.exists():
        with np.load(processed) as data:
            return int(data["sequences"].shape[0])
    raw_depth = msa_raw_depth(out_dir, record)
    if raw_depth == "":
        return ""
    return min(int(raw_depth), 8192)


def msa_raw_depth(out_dir: Path, record: FastaRecord) -> int | str:
    path = baseline_msa_path(out_dir, record)
    if not path.exists():
        return ""
    with path.open() as handle:
        return max(0, sum(1 for _ in handle) - 1)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def collect_rows(out_dir: Path, records: list[FastaRecord]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    conditions = [BASELINE, *GRID_CONDITIONS, *HIGH_STEP_CONDITIONS]
    for record in records:
        depth = msa_depth(out_dir, record)
        raw_depth = msa_raw_depth(out_dir, record)
        for condition in conditions:
            pred_dir = prediction_dir(out_dir, condition.name, record)
            conf_paths = sorted(pred_dir.glob(f"confidence_{record.job_id}_model_*.json"))
            if not conf_paths:
                if not condition.include_missing:
                    continue
                rows.append(
                    {
                        "record_index": record.index,
                        "job_id": record.job_id,
                        "fasta_header": record.header,
                        "length": len(record.sequence),
                        "msa_depth": depth,
                        "msa_raw_depth": raw_depth,
                        "condition": condition.name,
                        "sampling_steps": condition.sampling_steps,
                        "diffusion_samples": condition.diffusion_samples,
                        "max_parallel_samples": condition.max_parallel_samples,
                        "recycling_steps": condition.recycling_steps,
                        "model_rank": "",
                        "status": "missing",
                        "confidence_score": "",
                        "complex_plddt": "",
                        "ptm": "",
                        "iptm": "",
                        "complex_iplddt": "",
                        "complex_pde": "",
                        "complex_ipde": "",
                        "structure_path": "",
                        "confidence_json": "",
                    }
                )
                continue

            for conf_path in conf_paths:
                match = re.search(r"_model_(\d+)\.json$", conf_path.name)
                rank = int(match.group(1)) if match else -1
                confidence = load_json(conf_path)
                structure = pred_dir / f"{record.job_id}_model_{rank}.cif"
                rows.append(
                    {
                        "record_index": record.index,
                        "job_id": record.job_id,
                        "fasta_header": record.header,
                        "length": len(record.sequence),
                        "msa_depth": depth,
                        "msa_raw_depth": raw_depth,
                        "condition": condition.name,
                        "sampling_steps": condition.sampling_steps,
                        "diffusion_samples": condition.diffusion_samples,
                        "max_parallel_samples": condition.max_parallel_samples,
                        "recycling_steps": condition.recycling_steps,
                        "model_rank": rank,
                        "status": "ok" if structure.exists() else "no_structure",
                        "confidence_score": confidence.get("confidence_score", ""),
                        "complex_plddt": confidence.get("complex_plddt", ""),
                        "ptm": confidence.get("ptm", ""),
                        "iptm": confidence.get("iptm", ""),
                        "complex_iplddt": confidence.get("complex_iplddt", ""),
                        "complex_pde": confidence.get("complex_pde", ""),
                        "complex_ipde": confidence.get("complex_ipde", ""),
                        "structure_path": rel_path(structure) if structure.exists() else "",
                        "confidence_json": rel_path(conf_path),
                    }
                )
    return rows


def write_metadata(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    path = out_dir / "metadata.csv"
    fields = [
        "record_index",
        "job_id",
        "fasta_header",
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
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def copy_structures(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    structures_dir = out_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        if row["status"] != "ok" or row["structure_path"] == "":
            continue
        src = ROOT / row["structure_path"]
        dst_dir = structures_dir / row["condition"]
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = dst_dir / src.name
        if src.exists():
            shutil.copy2(src, dst)


def best_rows(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    best: dict[str, dict[str, Any]] = {}
    for row in rows:
        if row["status"] != "ok" or row["complex_plddt"] == "":
            continue
        current = best.get(row["job_id"])
        if current is None or float(row["complex_plddt"]) > float(current["complex_plddt"]):
            best[row["job_id"]] = row
    return best


def select_grid_records(
    rows: list[dict[str, Any]],
    records: list[FastaRecord],
    low_plddt: float,
    grid_pilot_size: int,
) -> list[FastaRecord]:
    baseline_rows = [
        row
        for row in rows
        if row["condition"] == BASELINE.name
        and row["status"] == "ok"
        and row["complex_plddt"] != ""
    ]
    best_baseline: dict[str, dict[str, Any]] = {}
    for row in baseline_rows:
        current = best_baseline.get(row["job_id"])
        if current is None or float(row["complex_plddt"]) > float(current["complex_plddt"]):
            best_baseline[row["job_id"]] = row

    selected_ids = {
        record.job_id
        for record in records
        if record.job_id not in best_baseline
        or float(best_baseline[record.job_id]["complex_plddt"]) < low_plddt
    }
    ranked = sorted(
        best_baseline.values(),
        key=lambda row: (float(row["complex_plddt"]), -int(row["length"])),
    )
    for row in ranked[:grid_pilot_size]:
        selected_ids.add(row["job_id"])
    return [record for record in records if record.job_id in selected_ids]


def quality_label(row: dict[str, Any] | None) -> str:
    if row is None:
        return "failed"
    plddt = float(row["complex_plddt"])
    if plddt >= 0.70:
        return "folded"
    if plddt >= 0.50:
        return "low confidence"
    return "failed"


def write_summary(
    out_dir: Path,
    rows: list[dict[str, Any]],
    records: list[FastaRecord],
    grid_records: list[FastaRecord],
    fasta: Path,
) -> None:
    best = best_rows(rows)
    folded = sum(1 for record in records if quality_label(best.get(record.job_id)) == "folded")
    low = sum(
        1 for record in records if quality_label(best.get(record.job_id)) == "low confidence"
    )
    failed = len(records) - folded - low
    grid_ids = {record.job_id for record in grid_records}
    retry_conditions = {condition.name for condition in FAILED_RETRY_CONDITIONS}
    high_step_conditions = {condition.name for condition in HIGH_STEP_CONDITIONS}
    pre_high_step_conditions = {
        BASELINE.name,
        *{condition.name for condition in GRID_CONDITIONS},
    }

    lines = [
        "# Test EVE Boltz Experiment",
        "",
        "## Overview",
        "",
        f"- Input FASTA: `{rel_path(fasta)}`",
        f"- Targets: {len(records)}",
        f"- Folded at complex_plddt >= 0.70: {folded}",
        f"- Low confidence, 0.50 <= complex_plddt < 0.70: {low}",
        f"- Failed or very low confidence: {failed}",
        "",
        "The starting point was the local H5 setting: "
        "`diffusion_samples=6`, `max_parallel_samples=2`, "
        "`sampling_steps=40`, `recycling_steps=3`. The grid then increased "
        "sampling steps, diffusion samples, or recycling for the lowest-confidence "
        "baseline targets.",
        "",
        "## Best Result Per Target",
        "",
        "| # | Target | Length | MSA used | MSA raw | Best status | Best config | Best sample rank | complex_plddt | ptm | confidence_score |",
        "|---:|---|---:|---:|---:|---|---|---:|---:|---:|---:|",
    ]
    for record in records:
        row = best.get(record.job_id)
        depth = msa_depth(out_dir, record)
        if row is None:
            lines.append(
                f"| {record.index} | {md_escape(record.header)} | {len(record.sequence)} | {depth or '-'} | {msa_raw_depth(out_dir, record) or '-'} | failed | - | - | - | - | - |"
            )
            continue
        lines.append(
            "| "
            f"{record.index} | {md_escape(record.header)} | {len(record.sequence)} | "
            f"{row['msa_depth'] or '-'} | {row['msa_raw_depth'] or '-'} | "
            f"{quality_label(row)} | {row['condition']} | {row['model_rank']} | "
            f"{float(row['complex_plddt']):.3f} | {float(row['ptm']):.3f} | "
            f"{float(row['confidence_score']):.3f} |"
        )

    lines += [
        "",
        "## Grid Search Targets",
        "",
        "| Target | Baseline pLDDT | Best grid pLDDT | Best grid config | Delta |",
        "|---|---:|---:|---|---:|",
    ]
    for record in records:
        if record.job_id not in grid_ids:
            continue
        baseline = best_for(rows, record.job_id, BASELINE.name)
        grid_best = best_for_any(rows, record.job_id, {c.name for c in GRID_CONDITIONS})
        baseline_plddt = float(baseline["complex_plddt"]) if baseline else None
        grid_plddt = float(grid_best["complex_plddt"]) if grid_best else None
        delta = (
            f"{grid_plddt - baseline_plddt:.3f}"
            if baseline_plddt is not None and grid_plddt is not None
            else "-"
        )
        lines.append(
            "| "
            f"{md_escape(record.header)} | {fmt_optional(baseline_plddt)} | "
            f"{fmt_optional(grid_plddt)} | "
            f"{grid_best['condition'] if grid_best else '-'} | {delta} |"
        )

    lines += [
        "",
        "## Failed-Target Retry",
        "",
        "The very low-confidence targets were retried with "
        "`sampling_steps=200`, `diffusion_samples=6`, "
        "`max_parallel_samples=2`, `recycling_steps=5`.",
        "",
        "| Target | Previous best pLDDT | Retry best pLDDT | Delta |",
        "|---|---:|---:|---:|",
    ]
    for record in records:
        retry_best = best_for_any(rows, record.job_id, retry_conditions)
        if retry_best is None:
            continue
        previous_best = best_for_any(
            rows,
            record.job_id,
            {BASELINE.name, *{condition.name for condition in GRID_CONDITIONS}},
        )
        previous_plddt = float(previous_best["complex_plddt"]) if previous_best else None
        retry_plddt = float(retry_best["complex_plddt"])
        delta = (
            f"{retry_plddt - previous_plddt:.3f}"
            if previous_plddt is not None
            else "-"
        )
        lines.append(
            "| "
            f"{md_escape(record.header)} | {fmt_optional(previous_plddt)} | "
            f"{retry_plddt:.3f} | {delta} |"
        )

    lines += [
        "",
        "## 200-Step Check",
        "",
        "All targets were checked with `sampling_steps=200`, "
        "`diffusion_samples=6`, `max_parallel_samples=2`, "
        "`recycling_steps=5`. The two very low-confidence targets use the "
        "`retry_s200_n6_p2_r5` outputs above; the rest use `all_s200_n6_p2_r5`.",
        "",
        "| Target | Previous best pLDDT | Best 200-step pLDDT | 200-step config | Delta |",
        "|---|---:|---:|---|---:|",
    ]
    for record in records:
        high_step_best = best_for_any(rows, record.job_id, high_step_conditions)
        if high_step_best is None:
            continue
        previous_best = best_for_any(rows, record.job_id, pre_high_step_conditions)
        previous_plddt = float(previous_best["complex_plddt"]) if previous_best else None
        high_step_plddt = float(high_step_best["complex_plddt"])
        delta = (
            f"{high_step_plddt - previous_plddt:.3f}"
            if previous_plddt is not None
            else "-"
        )
        lines.append(
            "| "
            f"{md_escape(record.header)} | {fmt_optional(previous_plddt)} | "
            f"{high_step_plddt:.3f} | {high_step_best['condition']} | {delta} |"
        )

    lines += [
        "",
        "## Files",
        "",
        f"- Metadata CSV: `{rel_path(out_dir / 'metadata.csv')}`",
        f"- Collected structures: `{rel_path(out_dir / 'structures')}/`",
        f"- Raw Boltz outputs: `{rel_path(out_dir / 'runs')}/`",
        f"- Logs: `{rel_path(out_dir / 'logs')}/`",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")


def best_for(rows: list[dict[str, Any]], job_id: str, condition: str) -> dict[str, Any] | None:
    return best_for_any(rows, job_id, {condition})


def best_for_any(
    rows: list[dict[str, Any]], job_id: str, conditions: set[str]
) -> dict[str, Any] | None:
    best: dict[str, Any] | None = None
    for row in rows:
        if (
            row["job_id"] != job_id
            or row["condition"] not in conditions
            or row["status"] != "ok"
            or row["complex_plddt"] == ""
        ):
            continue
        if best is None or float(row["complex_plddt"]) > float(best["complex_plddt"]):
            best = row
    return best


def fmt_optional(value: float | None) -> str:
    return "-" if value is None else f"{value:.3f}"


def md_escape(value: str) -> str:
    return value.replace("|", "\\|")


def rel_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def write_repro(out_dir: Path, args: argparse.Namespace) -> None:
    script = ROOT / "experiments" / "run_test_eve_experiments.py"
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
        "--low-plddt",
        str(args.low_plddt),
        "--grid-pilot-size",
        str(args.grid_pilot_size),
    ]
    if args.override:
        cmd.append("--override")
    text = "#!/usr/bin/env bash\nset -euo pipefail\n\n"
    text += "cd " + shlex.quote(str(ROOT)) + "\n"
    text += shlex.join(cmd) + "\n"
    path = out_dir / "repro.sh"
    path.write_text(text)
    path.chmod(0o755)


def main() -> None:
    args = parse_args()
    args.fasta = args.fasta.resolve()
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    records = parse_fasta(args.fasta)
    inputs_dir = args.out_dir / "inputs"
    baseline_inputs = inputs_dir / "baseline"
    for record in records:
        write_yaml(record, baseline_inputs / f"{record.job_id}.yaml")

    if not args.only_summary:
        run_condition(
            args=args,
            records=records,
            condition=BASELINE,
            input_dir=baseline_inputs,
            out_dir=args.out_dir / "runs",
            logs_dir=args.out_dir / "logs",
        )

    rows = collect_rows(args.out_dir, records)
    grid_records = select_grid_records(rows, records, args.low_plddt, args.grid_pilot_size)

    grid_inputs = inputs_dir / "grid"
    for record in grid_records:
        msa_path = baseline_msa_path(args.out_dir, record)
        if msa_path.exists():
            write_yaml(record, grid_inputs / f"{record.job_id}.yaml", msa_path.resolve())

    all_200_inputs = inputs_dir / "all_s200"
    for record in records:
        msa_path = baseline_msa_path(args.out_dir, record)
        if msa_path.exists():
            write_yaml(record, all_200_inputs / f"{record.job_id}.yaml", msa_path.resolve())

    if not args.only_summary:
        runnable_grid_records = [
            record for record in grid_records if (grid_inputs / f"{record.job_id}.yaml").exists()
        ]
        for condition in GRID_CONDITIONS:
            run_condition(
                args=args,
                records=runnable_grid_records,
                condition=condition,
                input_dir=grid_inputs,
                out_dir=args.out_dir / "runs",
                logs_dir=args.out_dir / "logs",
            )

        rows = collect_rows(args.out_dir, records)
        current_best = best_rows(rows)
        retry_records = [
            record
            for record in records
            if quality_label(current_best.get(record.job_id)) == "failed"
            and (grid_inputs / f"{record.job_id}.yaml").exists()
        ]
        for condition in FAILED_RETRY_CONDITIONS:
            run_condition(
                args=args,
                records=retry_records,
                condition=condition,
                input_dir=grid_inputs,
                out_dir=args.out_dir / "runs",
                logs_dir=args.out_dir / "logs",
            )

        retry_condition_names = {condition.name for condition in FAILED_RETRY_CONDITIONS}
        all_200_records = [
            record
            for record in records
            if (all_200_inputs / f"{record.job_id}.yaml").exists()
            and not any(
                prediction_dir(args.out_dir, condition_name, record).exists()
                for condition_name in retry_condition_names
            )
        ]
        for condition in ALL_200_CONDITIONS:
            run_condition(
                args=args,
                records=all_200_records,
                condition=condition,
                input_dir=all_200_inputs,
                out_dir=args.out_dir / "runs",
                logs_dir=args.out_dir / "logs",
            )

    rows = collect_rows(args.out_dir, records)
    write_metadata(args.out_dir, rows)
    copy_structures(args.out_dir, rows)
    write_summary(args.out_dir, rows, records, grid_records, args.fasta)
    write_repro(args.out_dir, args)


if __name__ == "__main__":
    main()
