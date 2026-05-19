from __future__ import annotations

from collections import defaultdict
from pathlib import Path

try:
    from .common import VIEWER_ROOT
except ImportError:  # pragma: no cover
    from utils.common import VIEWER_ROOT  # type: ignore[no-redef]

SOURCE_SUMMARIES = {
    "Flu HA trimer ligand screen": "flu_trimer_ligand_experiments/summary.md",
    "Flu HA trimer single-ligand affinity": "flu_trimer_single_ligand_affinity/summary.md",
    "EVE structures": "test_eve_experiments/summary.md",
    "EVE multimers": "test_eve_experiments/multimers/summary.md",
}


def _sorted_unique(rows: list[dict[str, str]], field: str) -> list[str]:
    return sorted({row.get(field, "").strip() for row in rows if row.get(field, "").strip()})


def generated_summary(experiment: str, experiment_slug: str, rows: list[dict[str, str]]) -> str:
    targets = _sorted_unique(rows, "target_id")
    assemblies = _sorted_unique(rows, "assembly")
    configs = _sorted_unique(rows, "config")
    sampling = _sorted_unique(rows, "sampling_steps")
    diffusion = _sorted_unique(rows, "diffusion_samples")
    parallel = _sorted_unique(rows, "max_parallel_samples")
    recycling = _sorted_unique(rows, "recycling_steps")
    folders = _sorted_unique(rows, "folder")
    sample_ranks = _sorted_unique(rows, "sample_rank")
    notes = [row.get("note", "").strip() for row in rows if row.get("note", "").strip()]

    lines = [
        f"# {experiment}",
        "",
        "## Overview",
        "",
        f"- Experiment slug: `{experiment_slug}`",
        f"- Samples in viewer: {len(rows)}",
        f"- Unique targets: {len(targets)}",
        f"- Assemblies: {', '.join(assemblies) if assemblies else 'unknown'}",
        f"- Configs: {', '.join(f'`{item}`' for item in configs) if configs else 'unknown'}",
        "",
        "## Parameters",
        "",
        f"- Sampling steps: {', '.join(sampling) if sampling else 'unknown'}",
        f"- Diffusion samples: {', '.join(diffusion) if diffusion else 'unknown'}",
        f"- Max parallel samples: {', '.join(parallel) if parallel else 'unknown'}",
        f"- Recycling steps: {', '.join(recycling) if recycling else 'unknown'}",
        f"- Sample ranks present: {', '.join(sample_ranks) if sample_ranks else 'unknown'}",
        "",
        "## Published Bundles",
        "",
    ]
    for folder in folders:
        lines.append(f"- `{folder}`")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "This summary was generated automatically from the standalone viewer metadata because no curated experiment markdown was found for this experiment.",
        ]
    )
    if notes:
        lines.extend(["", "Representative note fields from the run metadata:", ""])
        for note in notes[:5]:
            lines.append(f"- {note}")
    return "\n".join(lines) + "\n"


def export_experiment_notes(
    rows: list[dict[str, str]],
    viewer_root: Path,
    workspace_root: Path,
) -> None:
    notes_root = viewer_root / "experiments"
    if notes_root.exists():
        for path in sorted(notes_root.glob("*")):
            if path.is_dir():
                for nested in sorted(path.glob("**/*"), reverse=True):
                    if nested.is_file():
                        nested.unlink()
                for nested in sorted(path.glob("**/*"), reverse=True):
                    if nested.is_dir():
                        nested.rmdir()
                path.rmdir()
            elif path.is_file():
                path.unlink()
    notes_root.mkdir(parents=True, exist_ok=True)

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row.get("experiment_slug", "")].append(row)

    for experiment_slug, experiment_rows in sorted(grouped.items()):
        if not experiment_slug:
            continue
        experiment = experiment_rows[0].get("experiment", experiment_slug)
        summary_path = notes_root / experiment_slug / "summary.md"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        source_rel = SOURCE_SUMMARIES.get(experiment)
        if source_rel:
            source_path = workspace_root / source_rel
            if source_path.exists():
                summary_path.write_text(source_path.read_text())
                continue
        summary_path.write_text(generated_summary(experiment, experiment_slug, experiment_rows))
