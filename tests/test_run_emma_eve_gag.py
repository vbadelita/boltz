from __future__ import annotations

# ruff: noqa: ANN001, D202, PT018, SLF001
import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

_MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "process" / "run_emma_eve_gag.py"
)
_SPEC = importlib.util.spec_from_file_location("run_emma_eve_gag", _MODULE_PATH)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)

DATABASE_COLUMNS = _MODULE.DATABASE_COLUMNS
_Target = _MODULE._Target
_build_targets = _MODULE.build_targets
_count_completed_samples = _MODULE.count_completed_samples
_refresh_viewer_database = _MODULE.refresh_viewer_database
_write_boltz_fastas = _MODULE.write_boltz_fastas


def test_build_targets_creates_stable_prefixed_ids(tmp_path) -> None:
    """Prepared inputs should keep stable numeric prefixes and unique slugs."""

    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">Alpha/Beta\nMPEPTIDE\n>Alpha Beta\nMAGAIN\n")

    targets = _build_targets(input_fasta, tmp_path / "inputs")

    assert [target.target_id for target in targets] == [
        "gag_001_alpha_beta",
        "gag_002_alpha_beta",
    ]


def test_count_completed_samples_ignores_missing_structure(tmp_path) -> None:
    """A confidence JSON alone should not count as a completed sample."""

    pred_dir = tmp_path / "predictions" / "gag_001_test"
    pred_dir.mkdir(parents=True)
    (pred_dir / "confidence_gag_001_test_model_0.json").write_text("{}")
    (pred_dir / "confidence_gag_001_test_model_1.json").write_text("{}")
    (pred_dir / "gag_001_test_model_1.cif").write_text("data")

    assert _count_completed_samples(pred_dir, "gag_001_test") == 1


def test_refresh_viewer_database_replaces_experiment_rows(tmp_path) -> None:
    """Refreshing the database should replace only this experiment's rows."""

    run_root = tmp_path / "runs" / "emma_eve_gag_s40_n6_p6_r5"
    target = _Target(
        index=1,
        source_header="Example gag",
        sequence="MPEPTIDE",
        target_id="gag_001_example_gag",
        input_fasta=run_root / "inputs" / "gag_001_example_gag.fasta",
    )
    _write_boltz_fastas([target])

    target_root = run_root / target.target_id / f"boltz_results_{target.target_id}"
    pred_dir = target_root / "predictions" / target.target_id
    pred_dir.mkdir(parents=True)
    (pred_dir / f"{target.target_id}_model_0.cif").write_text("data_example\n")
    (pred_dir / f"confidence_{target.target_id}_model_0.json").write_text(
        json.dumps(
            {
                "confidence_score": 0.75,
                "complex_plddt": 0.8,
                "ptm": 0.7,
                "iptm": 0.0,
                "ligand_iptm": 0.0,
                "protein_iptm": 0.0,
                "complex_iplddt": 0.0,
                "complex_pde": 1.2,
                "complex_ipde": 0.0,
            }
        )
    )

    record_dir = target_root / "processed"
    (record_dir / "records").mkdir(parents=True)
    (record_dir / "msa").mkdir(parents=True)
    (record_dir / "records" / f"{target.target_id}.json").write_text(
        json.dumps({"chains": [{"num_residues": 8}]})
    )
    np.savez(
        record_dir / "msa" / f"{target.target_id}_0.npz",
        sequences=np.zeros((42,), dtype=np.int8),
    )

    viewer_db = tmp_path / "results_viewer" / "database.csv"
    viewer_db.parent.mkdir(parents=True)
    with viewer_db.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DATABASE_COLUMNS)
        writer.writeheader()
        writer.writerow(
            {
                "experiment": "Old experiment",
                "target_id": "old",
                "target_label": "old",
                "assembly": "monomer",
                "kind": "monomer",
                "length": "1",
                "msa_depth": "",
                "msa_raw_depth": "",
                "msa_csv": "",
                "config": "old",
                "sampling_steps": "1",
                "diffusion_samples": "1",
                "max_parallel_samples": "1",
                "recycling_steps": "1",
                "sampling_steps_affinity": "",
                "diffusion_samples_affinity": "",
                "sample_rank": "0",
                "status": "ok",
                "confidence_score": "1",
                "complex_plddt": "1",
                "ptm": "1",
                "iptm": "0",
                "ligand_iptm": "0",
                "protein_iptm": "0",
                "complex_iplddt": "0",
                "complex_pde": "0",
                "complex_ipde": "0",
                "affinity_pred_value": "",
                "affinity_probability_binary": "",
                "structure_path": "old.cif",
                "confidence_json": "old.json",
                "affinity_json": "",
                "note": "",
            }
        )

    _refresh_viewer_database(
        viewer_db=viewer_db,
        run_root=run_root,
        targets=[target],
        experiment="Emma EVE gag",
        config_name="emma_eve_gag_s40_n6_p6_r5",
        sampling_steps=40,
        diffusion_samples=6,
        max_parallel_samples=6,
        recycling_steps=5,
    )

    with viewer_db.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == len({"Old experiment", "Emma EVE gag"})
    assert {row["experiment"] for row in rows} == {"Old experiment", "Emma EVE gag"}
    new_row = next(row for row in rows if row["experiment"] == "Emma EVE gag")
    assert new_row["target_id"] == target.target_id
    assert new_row["sample_rank"] == "0"
    assert new_row["msa_depth"] == "42"
