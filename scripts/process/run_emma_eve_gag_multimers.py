#!/usr/bin/env python3
# ruff: noqa: D202, PLR0915, S603
"""Run reconstructed Emma EVE gag multimers as sequential Boltz jobs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np

DEFAULT_ALIGNDATA = Path("emma_eve_gag/aligndata.json")
DEFAULT_INPUT_FASTA = Path("emma_eve_gag/gag_sample_101.fasta")
DEFAULT_EXPERIMENT = "Emma EVE gag reconstructed multimers"
DEFAULT_CONFIG = "emma_eve_gag_multimers_s40_n6_p1_r5"
DEFAULT_LINKED_EXPERIMENT = "Emma EVE gag reconstructed multimers linked capsid"
DEFAULT_LINKED_CONFIG = "emma_eve_gag_multimers_linked_capsid_s40_n6_p1_r5"
DEFAULT_VIEWER_DB = None
COMPLEX_OOM_THRESHOLD = 1800
DEFAULT_MPLCONFIGDIR = Path("/tmp/matplotlib-boltz")
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
CHAIN_IDS = tuple("ABCDEFGHIJKLMNOPQRSTUVWXYZ")


@dataclass(frozen=True)
class DomainHit:
    """A contiguous query-space domain hit extracted from CD-search output."""

    domain_name: str
    start: int
    end: int

    @property
    def length(self) -> int:
        """Return the hit length in residues."""

        return self.end - self.start + 1


@dataclass(frozen=True)
class MultimerTarget:
    """One reconstructed multimer prepared as an individual Boltz input."""

    index: int
    source_header: str
    target_id: str
    input_yaml: Path
    class_name: str
    assembly_label: str
    total_length: int
    note: str


class _GracefulStopError(Exception):
    """Raised when the user interrupts the batch."""


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the multimer batch runner."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligndata", type=Path, default=DEFAULT_ALIGNDATA)
    parser.add_argument("--input-fasta", type=Path, default=DEFAULT_INPUT_FASTA)
    parser.add_argument("--experiment", default=DEFAULT_EXPERIMENT)
    parser.add_argument("--config-name", default=DEFAULT_CONFIG)
    parser.add_argument("--run-root", type=Path, default=None)
    parser.add_argument("--viewer-db", type=Path, default=DEFAULT_VIEWER_DB)
    parser.add_argument("--cache", default="~/.boltz")
    parser.add_argument("--model", default="boltz2", choices=["boltz1", "boltz2"])
    parser.add_argument("--devices", type=int, default=1)
    parser.add_argument("--accelerator", default="gpu")
    parser.add_argument("--sampling-steps", type=int, default=40)
    parser.add_argument("--diffusion-samples", type=int, default=6)
    parser.add_argument("--max-parallel-samples", type=int, default=1)
    parser.add_argument("--recycling-steps", type=int, default=5)
    parser.add_argument("--num-workers", type=int, default=2)
    parser.add_argument(
        "--preprocessing-threads",
        type=int,
        default=os.cpu_count() or 1,
    )
    parser.add_argument("--max-msa-seqs", type=int, default=8192)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--use-msa-server", action="store_true", default=True)
    parser.add_argument(
        "--no-use-msa-server",
        dest="use_msa_server",
        action="store_false",
    )
    parser.add_argument("--msa-server-url", default="https://api.colabfold.com")
    parser.add_argument("--msa-pairing-strategy", default="greedy")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--only-target", action="append", default=[])
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--refresh-viewer-only", action="store_true")
    parser.add_argument("--combine-class2-capsid", action="store_true")
    parser.add_argument("--python-executable", type=Path, default=None)
    parser.add_argument("--use-lmi4boltz", action="store_true")
    parser.add_argument("--lmi4boltz-src", type=Path, default=Path("lmi4boltz/src"))
    parser.add_argument("--mplconfigdir", type=Path, default=DEFAULT_MPLCONFIGDIR)
    parser.add_argument("--chunk-size-transition-z", type=int, default=None)
    parser.add_argument("--chunk-size-transition-msa", type=int, default=None)
    parser.add_argument("--chunk-size-outer-product", type=int, default=None)
    parser.add_argument("--chunk-size-tri-attn", type=int, default=None)
    parser.add_argument("--triangle-mult-gate-nchunks", type=int, default=None)
    parser.add_argument("--chunk-size-threshold", type=int, default=None)
    parser.add_argument("--use-bfloat16", action="store_true")
    return parser.parse_args()


def parse_simple_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse a minimal FASTA file into `(header, sequence)` pairs."""

    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence_parts)))
            header = line[1:].strip()
            sequence_parts = []
            continue
        sequence_parts.append(line)
    if header is not None:
        records.append((header, "".join(sequence_parts)))
    return records


def load_cdsearch_objects(path: Path) -> list[dict]:
    """Load concatenated JSON objects emitted by the CD-search batch tool."""

    decoder = json.JSONDecoder()
    text = path.read_text()
    objects: list[dict] = []
    index = 0
    while True:
        while index < len(text) and text[index].isspace():
            index += 1
        if index >= len(text):
            return objects
        obj, index = decoder.raw_decode(text, index)
        objects.append(obj)


def extract_title(seq_annot: dict) -> str | None:
    """Return the record title if present."""

    for desc in seq_annot.get("desc", []):
        title = desc.get("title")
        if title:
            return title
    return None


def title_to_fasta_header(title: str) -> str:
    """Strip the `Q#n - >` prefix from a CD-search title."""

    _, _, suffix = title.partition(" - >")
    return suffix if suffix else title.removeprefix(">")


def denseg_query_span(denseg: dict) -> tuple[int, int]:
    """Collapse a dense-seg alignment into a contiguous 0-based query span."""

    dim = denseg["dim"]
    segments: list[tuple[int, int]] = []
    for seg_index, seg_len in enumerate(denseg["lens"]):
        query_start = denseg["starts"][seg_index * dim]
        if query_start >= 0:
            segments.append((query_start, query_start + seg_len))
    if not segments:
        msg = "denseg alignment does not contain query residues"
        raise ValueError(msg)
    return segments[0][0], segments[-1][1] - 1


def collect_raw_domain_hits(path: Path) -> dict[str, dict[str, DomainHit]]:
    """Collect raw Gag domain hits keyed by FASTA header and domain name."""

    wanted = {"Gag_MA", "Gag_p30", "Gag_p10", "Gag_p24", "Gag_p24_C"}
    hits_by_header: dict[str, dict[str, DomainHit]] = {}
    for obj in load_cdsearch_objects(path):
        seq_annot = obj.get("Seq_annot", {})
        title = extract_title(seq_annot)
        if title is None:
            continue
        fasta_header = title_to_fasta_header(title)
        domain_hits = hits_by_header.setdefault(fasta_header, {})
        for align in seq_annot.get("data", {}).get("align", []):
            labels = [
                item["str"]
                for item in align.get("id", [])
                if isinstance(item, dict) and "str" in item
            ]
            matching_domain = next((label for label in labels if label in wanted), None)
            if matching_domain is None or matching_domain in domain_hits:
                continue
            start, end = denseg_query_span(align["segs"]["denseg"])
            domain_hits[matching_domain] = DomainHit(
                domain_name=matching_domain,
                start=start,
                end=end,
            )
    return {header: hits for header, hits in hits_by_header.items() if hits}


def slugify(text: str, limit: int = 48) -> str:
    """Convert a free-text header into a stable filesystem slug."""

    slug = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
    slug = re.sub(r"_+", "_", slug)
    if not slug:
        return "unnamed"
    return slug[:limit].strip("_") or "unnamed"


def sequence_lookup(fasta_path: Path) -> dict[str, str]:
    """Return a mapping from FASTA header to sequence."""

    return {header: sequence for header, sequence in parse_simple_fasta(fasta_path)}


def slice_sequence(sequence: str, hit: DomainHit) -> str:
    """Slice the annotated domain span from the full sequence."""

    return sequence[hit.start : hit.end + 1]


def alloc_chain_ids(count: int, offset: int) -> tuple[list[str], int]:
    """Allocate deterministic chain ids from the alphabet."""

    next_offset = offset + count
    return list(CHAIN_IDS[offset:next_offset]), next_offset


def yaml_entity_block(ids: list[str], sequence: str) -> str:
    """Render one protein entity block for a Boltz YAML input."""

    if len(ids) == 1:
        id_repr = ids[0]
    else:
        id_repr = f"[{', '.join(ids)}]"
    return (
        "  - protein:\n"
        f"      id: {id_repr}\n"
        f"      sequence: {sequence}\n"
    )


def multimer_components_for_header(
    header: str,
    sequence: str,
    domain_hits: dict[str, DomainHit],
    combine_class2_capsid: bool,
) -> tuple[str, list[tuple[str, str, int]], str]:
    """Return class, components, and note for one eligible source sequence."""

    if {"Gag_MA", "Gag_p30"} <= set(domain_hits):
        return (
            "class_1",
            [
                ("matrix", slice_sequence(sequence, domain_hits["Gag_MA"]), 3),
                ("capsid", slice_sequence(sequence, domain_hits["Gag_p30"]), 6),
            ],
            "class_1: 3x Gag_MA + 6x Gag_p30",
        )
    if {"Gag_p10", "Gag_p24", "Gag_p24_C"} <= set(domain_hits):
        if combine_class2_capsid:
            capsid_start = domain_hits["Gag_p24"].start
            capsid_end = domain_hits["Gag_p24_C"].end
            linked_capsid = sequence[capsid_start : capsid_end + 1]
            return (
                "class_2",
                [
                    ("matrix", slice_sequence(sequence, domain_hits["Gag_p10"]), 3),
                    ("capsid", linked_capsid, 6),
                ],
                "class_2: linked 6x capsid span from Gag_p24 start through Gag_p24_C end, including linker; 3x Gag_p10 + 6x linked_capsid",
            )
        return (
            "class_2",
            [
                ("matrix", slice_sequence(sequence, domain_hits["Gag_p10"]), 3),
                ("capsid_n_terminal", slice_sequence(sequence, domain_hits["Gag_p24"]), 6),
                (
                    "capsid_c_terminal",
                    slice_sequence(sequence, domain_hits["Gag_p24_C"]),
                    6,
                ),
            ],
            "class_2: assumed 3x Gag_p10 + 6x Gag_p24 + 6x Gag_p24_C",
        )
    msg = f"Header is not a complete class-1 or class-2 target: {header}"
    raise ValueError(msg)


def maybe_reduce_stoichiometry(
    components: list[tuple[str, str, int]],
) -> tuple[list[tuple[str, str, int]], str]:
    """Reduce to one copy of each component if the full complex is too large."""

    total_length = sum(len(sequence) * copies for _, sequence, copies in components)
    if total_length <= COMPLEX_OOM_THRESHOLD:
        return components, "full_stoichiometry"
    reduced = [(label, sequence, 1) for label, sequence, _ in components]
    return reduced, "reduced_to_1x_each_due_to_total_length"


def build_yaml_text(components: list[tuple[str, str, int]]) -> tuple[str, str, int]:
    """Render the Boltz YAML for one multimer input."""

    lines = ["version: 1", "sequences:"]
    assembly_parts: list[str] = []
    chain_offset = 0
    total_length = 0
    for _, sequence, copies in components:
        ids, chain_offset = alloc_chain_ids(copies, chain_offset)
        lines.append(yaml_entity_block(ids, sequence).rstrip())
        assembly_parts.append(f"{copies}x{len(sequence)}")
        total_length += len(sequence) * copies
    return "\n".join(lines) + "\n", " + ".join(assembly_parts), total_length


def build_targets(
    aligndata_path: Path,
    fasta_path: Path,
    inputs_dir: Path,
    combine_class2_capsid: bool = False,
) -> list[MultimerTarget]:
    """Create multimer targets for every sequence that contains all required parts."""

    sequences = sequence_lookup(fasta_path)
    raw_hits = collect_raw_domain_hits(aligndata_path)
    seen_ids: set[str] = set()
    targets: list[MultimerTarget] = []
    index = 1
    for header, sequence in sequences.items():
        if header not in raw_hits:
            continue
        domain_hits = raw_hits[header]
        if not (
            {"Gag_MA", "Gag_p30"} <= set(domain_hits)
            or {"Gag_p10", "Gag_p24", "Gag_p24_C"} <= set(domain_hits)
        ):
            continue
        class_name, components, note = multimer_components_for_header(
            header, sequence, domain_hits, combine_class2_capsid
        )
        adjusted_components, size_mode = maybe_reduce_stoichiometry(components)
        yaml_text, assembly_label, total_length = build_yaml_text(adjusted_components)
        target_id = f"{class_name}_{index:03d}_{slugify(header)}"
        suffix = 2
        base_target_id = target_id
        while target_id in seen_ids:
            target_id = f"{base_target_id}_{suffix}"
            suffix += 1
        seen_ids.add(target_id)
        input_yaml = inputs_dir / f"{target_id}.yaml"
        input_yaml.parent.mkdir(parents=True, exist_ok=True)
        input_yaml.write_text(yaml_text)
        targets.append(
            MultimerTarget(
                index=index,
                source_header=header,
                target_id=target_id,
                input_yaml=input_yaml,
                class_name=class_name,
                assembly_label=assembly_label,
                total_length=total_length,
                note=f"{note}; size_mode={size_mode}",
            )
        )
        index += 1
    return targets


def write_targets_manifest(targets: list[MultimerTarget], manifest_path: Path) -> None:
    """Persist the prepared target manifest for later inspection."""

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "index",
                "target_id",
                "source_header",
                "class_name",
                "assembly_label",
                "total_length",
                "input_yaml",
                "note",
            ],
        )
        writer.writeheader()
        for target in targets:
            writer.writerow(
                {
                    "index": target.index,
                    "target_id": target.target_id,
                    "source_header": target.source_header,
                    "class_name": target.class_name,
                    "assembly_label": target.assembly_label,
                    "total_length": target.total_length,
                    "input_yaml": target.input_yaml.as_posix(),
                    "note": target.note,
                }
            )


def prediction_root(run_root: Path, target_id: str) -> Path:
    """Return the root Boltz output directory for one target."""

    return run_root / target_id / f"boltz_results_{target_id}"


def predictions_dir(run_root: Path, target_id: str) -> Path:
    """Return the structure prediction directory for one target."""

    return prediction_root(run_root, target_id) / "predictions" / target_id


def processed_record_path(run_root: Path, target_id: str) -> Path:
    """Return the processed-record JSON path for one target."""

    return (
        prediction_root(run_root, target_id)
        / "processed"
        / "records"
        / f"{target_id}.json"
    )


def processed_msa_path(run_root: Path, target_id: str) -> Path:
    """Return the processed MSA `.npz` path for the first chain."""

    return (
        prediction_root(run_root, target_id)
        / "processed"
        / "msa"
        / f"{target_id}_0.npz"
    )


def count_completed_samples(predictions_path: Path, target_id: str) -> int:
    """Count sample ranks that have both confidence JSON and structure output."""

    if not predictions_path.exists():
        return 0
    complete = 0
    for confidence_path in sorted(
        predictions_path.glob(f"confidence_{target_id}_model_*.json")
    ):
        match = re.search(r"_model_(\d+)\.json$", confidence_path.name)
        if match is None:
            continue
        model_index = int(match.group(1))
        structure_path = predictions_path / f"{target_id}_model_{model_index}.cif"
        if structure_path.exists():
            complete += 1
    return complete


def should_skip_target(run_root: Path, target: MultimerTarget, diffusion_samples: int) -> bool:
    """Return whether all requested samples already exist for the target."""

    return (
        count_completed_samples(
            predictions_dir(run_root, target.target_id), target.target_id
        )
        >= diffusion_samples
    )


def atomic_write_csv(
    path: Path, rows: list[dict[str, str]], fieldnames: list[str]
) -> None:
    """Write a CSV atomically to avoid half-written viewer files."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        temp_name = handle.name
    Path(temp_name).replace(path)


def relative_to_repo(path: Path) -> str:
    """Return a repository-relative POSIX path for viewer links."""

    resolved = path.resolve()
    repo_root = Path.cwd().resolve()
    try:
        return resolved.relative_to(repo_root).as_posix()
    except ValueError:
        return resolved.as_posix()


def top_level_folder(path: str) -> str:
    """Return the first path segment for a repository-relative path."""

    return Path(path).parts[0] if path else ""


def load_existing_database(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Load the current viewer database if present."""

    if not path.exists():
        return DATABASE_COLUMNS, []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or DATABASE_COLUMNS
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def load_record_metadata(run_root: Path, target: MultimerTarget) -> tuple[int | None, int | None]:
    """Load total length and first-chain MSA depth from Boltz processed outputs."""

    record_path = processed_record_path(run_root, target.target_id)
    msa_path = processed_msa_path(run_root, target.target_id)
    length = None
    msa_depth = None
    if record_path.exists():
        record = json.loads(record_path.read_text())
        chains = record.get("chains") or []
        if chains:
            length = sum(int(chain["num_residues"]) for chain in chains)
    if msa_path.exists():
        with np.load(msa_path) as msa_data:
            msa_depth = int(msa_data["sequences"].shape[0])
    return length, msa_depth


def build_viewer_rows(
    run_root: Path,
    targets: list[MultimerTarget],
    experiment: str,
    config_name: str,
    sampling_steps: int,
    diffusion_samples: int,
    max_parallel_samples: int,
    recycling_steps: int,
) -> list[dict[str, str]]:
    """Build viewer rows from every completed sample currently on disk."""

    rows: list[dict[str, str]] = []
    for target in targets:
        pred_dir = predictions_dir(run_root, target.target_id)
        if not pred_dir.exists():
            continue
        inferred_length, msa_depth = load_record_metadata(run_root, target)
        length = inferred_length or target.total_length
        for confidence_path in sorted(
            pred_dir.glob(f"confidence_{target.target_id}_model_*.json")
        ):
            match = re.search(r"_model_(\d+)\.json$", confidence_path.name)
            if match is None:
                continue
            sample_rank = int(match.group(1))
            structure_path = pred_dir / f"{target.target_id}_model_{sample_rank}.cif"
            if not structure_path.exists():
                continue
            confidence = json.loads(confidence_path.read_text())
            rows.append(
                {
                    "experiment": experiment,
                    "folder": top_level_folder(relative_to_repo(structure_path)),
                    "target_id": target.target_id,
                    "target_label": target.source_header,
                    "assembly": target.assembly_label,
                    "kind": "multimer",
                    "length": str(length),
                    "msa_depth": "" if msa_depth is None else str(msa_depth),
                    "msa_raw_depth": "" if msa_depth is None else str(msa_depth),
                    "msa_csv": "",
                    "config": config_name,
                    "sampling_steps": str(sampling_steps),
                    "diffusion_samples": str(diffusion_samples),
                    "max_parallel_samples": str(max_parallel_samples),
                    "recycling_steps": str(recycling_steps),
                    "sampling_steps_affinity": "",
                    "diffusion_samples_affinity": "",
                    "sample_rank": str(sample_rank),
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
                    "affinity_pred_value": "",
                    "affinity_probability_binary": "",
                    "structure_path": relative_to_repo(structure_path),
                    "confidence_json": relative_to_repo(confidence_path),
                    "affinity_json": "",
                    "note": f"{relative_to_repo(target.input_yaml)}; {target.note}",
                }
            )
    rows.sort(key=lambda row: (row["target_id"], int(row["sample_rank"])))
    return rows


def refresh_viewer_database(
    viewer_db: Path,
    run_root: Path,
    targets: list[MultimerTarget],
    experiment: str,
    config_name: str,
    sampling_steps: int,
    diffusion_samples: int,
    max_parallel_samples: int,
    recycling_steps: int,
) -> int:
    """Replace this experiment's rows inside the target viewer database."""

    fieldnames, existing_rows = load_existing_database(viewer_db)
    kept_rows = [row for row in existing_rows if row.get("experiment") != experiment]
    new_rows = build_viewer_rows(
        run_root=run_root,
        targets=targets,
        experiment=experiment,
        config_name=config_name,
        sampling_steps=sampling_steps,
        diffusion_samples=diffusion_samples,
        max_parallel_samples=max_parallel_samples,
        recycling_steps=recycling_steps,
    )
    atomic_write_csv(viewer_db, kept_rows + new_rows, list(fieldnames))
    return len(new_rows)


def shlex_quote(text: str) -> str:
    """Quote one shell argument for the repro script."""

    return shlex.quote(text)


def write_repro_script(run_root: Path, argv: list[str]) -> None:
    """Write a small repro script next to the run outputs."""

    repro_path = run_root / "repro.sh"
    command = " ".join(shlex_quote(arg) for arg in argv)
    repro_path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "# Reproduce the Emma EVE gag multimer batch from the repository root.\n\n"
        f"{command}\n"
    )
    repro_path.chmod(0o755)


def build_boltz_command(
    args: argparse.Namespace,
    target: MultimerTarget,
    override: bool,
) -> list[str]:
    """Build the Boltz CLI command for one target."""

    if args.use_lmi4boltz:
        python_executable = args.python_executable or (
            Path(".venv/bin/python") if Path(".venv/bin/python").exists() else Path(sys.executable)
        )
        command = [
            str(python_executable),
            "-m",
            "boltz.main",
            "predict",
            str(target.input_yaml),
        ]
    else:
        command = [
            "uv",
            "run",
            "boltz",
            "predict",
            str(target.input_yaml),
        ]
    command.extend(
        [
            "--out_dir",
            str((args.run_root or Path("runs") / args.config_name) / target.target_id),
            "--cache",
            args.cache,
            "--devices",
            str(args.devices),
            "--accelerator",
            args.accelerator,
            "--recycling_steps",
            str(args.recycling_steps),
            "--sampling_steps",
            str(args.sampling_steps),
            "--diffusion_samples",
            str(args.diffusion_samples),
            "--max_parallel_samples",
            str(args.max_parallel_samples),
            "--num_workers",
            str(args.num_workers),
            "--preprocessing-threads",
            str(args.preprocessing_threads),
            "--max_msa_seqs",
            str(args.max_msa_seqs),
            "--model",
            args.model,
            "--write_full_pae",
            "--write_full_pde",
            "--output_format",
            "mmcif",
        ]
    )
    if args.use_msa_server:
        command.extend(
            [
                "--use_msa_server",
                "--msa_server_url",
                args.msa_server_url,
                "--msa_pairing_strategy",
                args.msa_pairing_strategy,
            ]
        )
    if args.seed is not None:
        command.extend(["--seed", str(args.seed + target.index - 1)])
    if args.chunk_size_transition_z is not None:
        command.extend(["--chunk_size_transition_z", str(args.chunk_size_transition_z)])
    if args.chunk_size_transition_msa is not None:
        command.extend(["--chunk_size_transition_msa", str(args.chunk_size_transition_msa)])
    if args.chunk_size_outer_product is not None:
        command.extend(["--chunk_size_outer_product", str(args.chunk_size_outer_product)])
    if args.chunk_size_tri_attn is not None:
        command.extend(["--chunk_size_tri_attn", str(args.chunk_size_tri_attn)])
    if args.triangle_mult_gate_nchunks is not None:
        command.extend(
            [
                "--triangle_mult_gate_nchunks",
                str(args.triangle_mult_gate_nchunks),
            ]
        )
    if args.chunk_size_threshold is not None:
        command.extend(["--chunk_size_threshold", str(args.chunk_size_threshold)])
    if args.use_bfloat16:
        command.append("--use_bfloat16")
    if override:
        command.append("--override")
    return command


def build_boltz_env(args: argparse.Namespace) -> dict[str, str]:
    """Build the environment for the Boltz subprocess."""

    env = os.environ.copy()
    if args.use_lmi4boltz:
        lmi4boltz_src = args.lmi4boltz_src.resolve()
        existing_pythonpath = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = (
            f"{lmi4boltz_src}{os.pathsep}{existing_pythonpath}"
            if existing_pythonpath
            else str(lmi4boltz_src)
        )
        env["MPLCONFIGDIR"] = str(args.mplconfigdir)
    return env


def install_signal_handlers(stop_flag: dict[str, bool]) -> None:
    """Install signal handlers that request a clean stop after the target ends."""

    def _handler(signum: int, _frame: object) -> None:
        stop_flag["stop"] = True
        _warn(f"\nReceived signal {signum}; stopping after the current target.")

    signal.signal(signal.SIGINT, _handler)
    signal.signal(signal.SIGTERM, _handler)


def prepare_run(args: argparse.Namespace) -> tuple[Path, Path, list[MultimerTarget]]:
    """Prepare multimer inputs and bookkeeping files for the batch."""

    run_root = args.run_root or Path("runs") / args.config_name
    viewer_db = args.viewer_db or run_root / "viewer.csv"
    inputs_dir = run_root / "inputs"
    targets = build_targets(
        args.aligndata,
        args.input_fasta,
        inputs_dir,
        combine_class2_capsid=args.combine_class2_capsid,
    )
    if args.only_target:
        wanted = set(args.only_target)
        targets = [target for target in targets if target.target_id in wanted]
    if args.limit is not None:
        targets = targets[: args.limit]
    write_targets_manifest(targets, run_root / "targets.csv")
    write_repro_script(run_root, sys.argv)
    return run_root, viewer_db, targets


def _info(message: str) -> None:
    """Write a normal progress message."""

    sys.stdout.write(f"{message}\n")
    sys.stdout.flush()


def _warn(message: str) -> None:
    """Write a warning or status message to stderr."""

    sys.stderr.write(f"{message}\n")
    sys.stderr.flush()


def _format_duration(seconds: float) -> str:
    """Format a duration in a compact human-readable form."""

    total_seconds = max(0, int(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{seconds:02d}s"
    if minutes:
        return f"{minutes}m{seconds:02d}s"
    return f"{seconds}s"


def _count_finished_targets(
    run_root: Path,
    targets: list[MultimerTarget],
    diffusion_samples: int,
) -> int:
    """Count targets that already have the full requested sample set."""

    return sum(
        should_skip_target(run_root, target, diffusion_samples) for target in targets
    )


def _count_finished_samples(run_root: Path, targets: list[MultimerTarget]) -> int:
    """Count completed structure samples across the whole batch."""

    return sum(
        count_completed_samples(
            predictions_dir(run_root, target.target_id), target.target_id
        )
        for target in targets
    )


def run_batch(
    args: argparse.Namespace,
    run_root: Path,
    viewer_db: Path,
    targets: list[MultimerTarget],
) -> int:
    """Run the prepared targets sequentially with viewer refreshes between them."""

    stop_flag = {"stop": False}
    install_signal_handlers(stop_flag)
    completed_targets = 0
    batch_started_at = time.monotonic()
    total_targets = len(targets)
    total_requested_samples = total_targets * args.diffusion_samples
    for position, target in enumerate(targets, start=1):
        pred_dir = predictions_dir(run_root, target.target_id)
        partial_count = count_completed_samples(pred_dir, target.target_id)
        finished_targets_before = _count_finished_targets(
            run_root, targets, args.diffusion_samples
        )
        finished_samples_before = _count_finished_samples(run_root, targets)
        batch_elapsed = _format_duration(time.monotonic() - batch_started_at)
        if should_skip_target(run_root, target, args.diffusion_samples):
            _info(
                f"[{position}/{total_targets}] skip {target.target_id}: "
                f"target {finished_targets_before}/{total_targets}, "
                f"samples {finished_samples_before}/{total_requested_samples}, "
                f"elapsed {batch_elapsed}, "
                f"already complete ({partial_count}/{args.diffusion_samples})"
            )
            completed_targets += 1
            continue

        override = pred_dir.exists()
        mode = "resume-from-partial" if override else "new"
        _info(
            f"[{position}/{total_targets}] run {target.target_id}: {mode}, "
            f"target {finished_targets_before}/{total_targets}, "
            f"samples {finished_samples_before}/{total_requested_samples}, "
            f"elapsed {batch_elapsed}, "
            f"models present {partial_count}/{args.diffusion_samples}"
        )

        log_path = run_root / "logs" / f"{target.target_id}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        command = build_boltz_command(args, target, override=override)
        target_started_at = time.monotonic()
        next_heartbeat_at = target_started_at + 60.0

        with log_path.open("a") as log_handle:
            log_handle.write(f"$ {' '.join(command)}\n")
            log_handle.flush()
            process = subprocess.Popen(
                command,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                env=build_boltz_env(args),
            )
            while True:
                try:
                    return_code = process.wait(timeout=1.0)
                    break
                except subprocess.TimeoutExpired:
                    now = time.monotonic()
                    if now >= next_heartbeat_at:
                        current_samples = count_completed_samples(
                            pred_dir, target.target_id
                        )
                        _info(
                            f"  still running {target.target_id}: "
                            "target elapsed "
                            f"{_format_duration(now - target_started_at)}, "
                            "batch elapsed "
                            f"{_format_duration(now - batch_started_at)}, "
                            "models finished "
                            f"{current_samples}/{args.diffusion_samples}, "
                            f"log {relative_to_repo(log_path)}"
                        )
                        next_heartbeat_at = now + 60.0
                    if stop_flag["stop"]:
                        process.send_signal(signal.SIGINT)
                        return_code = process.wait()
                        break

        rows_written = refresh_viewer_database(
            viewer_db=viewer_db,
            run_root=run_root,
            targets=targets,
            experiment=args.experiment,
            config_name=args.config_name,
            sampling_steps=args.sampling_steps,
            diffusion_samples=args.diffusion_samples,
            max_parallel_samples=args.max_parallel_samples,
            recycling_steps=args.recycling_steps,
        )
        finished_targets_after = _count_finished_targets(
            run_root, targets, args.diffusion_samples
        )
        finished_samples_after = _count_finished_samples(run_root, targets)
        _info(
            f"  viewer refreshed: {rows_written} rows for {args.experiment}; "
            f"progress target {finished_targets_after}/{total_targets}, "
            f"samples {finished_samples_after}/{total_requested_samples}, "
            f"target time {_format_duration(time.monotonic() - target_started_at)}"
        )

        if return_code != 0:
            _warn(f"  target failed: see {log_path}")
            if stop_flag["stop"]:
                raise _GracefulStopError
            return return_code

        completed_targets += 1
        if stop_flag["stop"]:
            raise _GracefulStopError

    _info(
        f"Finished {completed_targets}/{total_targets} targets in "
        f"{_format_duration(time.monotonic() - batch_started_at)}."
    )
    return 0


def main() -> int:
    """Prepare or execute the Emma EVE gag reconstructed multimer batch."""

    args = parse_args()
    if not args.aligndata.exists():
        msg = f"aligndata not found: {args.aligndata}"
        raise FileNotFoundError(msg)
    if not args.input_fasta.exists():
        msg = f"Input FASTA not found: {args.input_fasta}"
        raise FileNotFoundError(msg)

    run_root, viewer_db, targets = prepare_run(args)
    if not targets:
        _warn("No targets selected.")
        return 1

    if args.refresh_viewer_only:
        rows_written = refresh_viewer_database(
            viewer_db=viewer_db,
            run_root=run_root,
            targets=targets,
            experiment=args.experiment,
            config_name=args.config_name,
            sampling_steps=args.sampling_steps,
            diffusion_samples=args.diffusion_samples,
            max_parallel_samples=args.max_parallel_samples,
            recycling_steps=args.recycling_steps,
        )
        _info(f"Viewer database refreshed with {rows_written} rows.")
        return 0

    rows_written = refresh_viewer_database(
        viewer_db=viewer_db,
        run_root=run_root,
        targets=targets,
        experiment=args.experiment,
        config_name=args.config_name,
        sampling_steps=args.sampling_steps,
        diffusion_samples=args.diffusion_samples,
        max_parallel_samples=args.max_parallel_samples,
        recycling_steps=args.recycling_steps,
    )
    _info(
        f"Prepared {len(targets)} multimer targets from {args.input_fasta} "
        f"and found {rows_written} existing viewer rows."
    )

    if args.prepare_only:
        return 0

    try:
        return run_batch(args, run_root, viewer_db, targets)
    except _GracefulStopError:
        refresh_viewer_database(
            viewer_db=viewer_db,
            run_root=run_root,
            targets=targets,
            experiment=args.experiment,
            config_name=args.config_name,
            sampling_steps=args.sampling_steps,
            diffusion_samples=args.diffusion_samples,
            max_parallel_samples=args.max_parallel_samples,
            recycling_steps=args.recycling_steps,
        )
        _info("Stopped cleanly after refreshing the viewer database.")
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
