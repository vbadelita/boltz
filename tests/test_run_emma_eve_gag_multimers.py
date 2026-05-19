from __future__ import annotations

import importlib.util
import json
import sys
from argparse import Namespace
from pathlib import Path

_MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "process"
    / "run_emma_eve_gag_multimers.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "run_emma_eve_gag_multimers",
    _MODULE_PATH,
)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)

_DomainHit = _MODULE.DomainHit
_build_boltz_command = _MODULE.build_boltz_command
_build_boltz_env = _MODULE.build_boltz_env
_build_targets = _MODULE.build_targets
_build_yaml_text = _MODULE.build_yaml_text
_maybe_reduce_stoichiometry = _MODULE.maybe_reduce_stoichiometry
_multimer_components_for_header = _MODULE.multimer_components_for_header


def _align_record(title: str, domain_name: str, start: int, end: int) -> dict:
    return {
        "Seq_annot": {
            "desc": [{"title": title}],
            "data": {
                "align": [
                    {
                        "segs": {
                            "denseg": {
                                "dim": 2,
                                "numseg": 1,
                                "starts": [start, 0],
                                "lens": [end - start + 1],
                            }
                        },
                        "id": [{"str": domain_name}, {"str": f"{domain_name} superfamily"}],
                    }
                ]
            },
        }
    }


def test_maybe_reduce_stoichiometry_falls_back_for_large_complex() -> None:
    """Large complexes should reduce to one copy of each component."""

    adjusted, size_mode = _maybe_reduce_stoichiometry(
        [
            ("matrix", "A" * 100, 3),
            ("capsid", "B" * 260, 6),
        ]
    )

    assert size_mode == "reduced_to_1x_each_due_to_total_length"
    assert adjusted == [("matrix", "A" * 100, 1), ("capsid", "B" * 260, 1)]


def test_build_yaml_text_groups_identical_copies() -> None:
    """Grouped identical copies should render as Boltz YAML id lists."""

    yaml_text, assembly_label, total_length = _build_yaml_text(
        [
            ("matrix", "AAAA", 3),
            ("capsid", "BBBBB", 2),
        ]
    )

    assert "id: [A, B, C]" in yaml_text
    assert "id: [D, E]" in yaml_text
    assert assembly_label == "3x4 + 2x5"
    assert total_length == 22


def test_build_targets_creates_class_specific_yaml_inputs(tmp_path: Path) -> None:
    """Only complete class-1/class-2 inputs should be prepared."""

    fasta = tmp_path / "input.fasta"
    fasta.write_text(">class1\nABCDEFGHIJKL\n>class2\nMNOPQRSTUVWXYZ\n>partial\nABCDEFG\n")
    aligndata = tmp_path / "aligndata.json"
    records = [
        _align_record("Q#1 - >class1", "Gag_MA", 0, 2),
        _align_record("Q#1 - >class1", "Gag_p30", 3, 8),
        _align_record("Q#2 - >class2", "Gag_p10", 0, 1),
        _align_record("Q#2 - >class2", "Gag_p24", 2, 4),
        _align_record("Q#2 - >class2", "Gag_p24_C", 5, 7),
        _align_record("Q#3 - >partial", "Gag_p10", 0, 1),
    ]
    aligndata.write_text("".join(json.dumps(record) for record in records))

    targets = _build_targets(aligndata, fasta, tmp_path / "inputs")

    assert [target.class_name for target in targets] == ["class_1", "class_2"]
    assert "GHI" in targets[0].input_yaml.read_text()
    assert "id: [A, B, C]" in targets[1].input_yaml.read_text()


def test_linked_class2_capsid_includes_intervening_linker() -> None:
    """Linked class-2 mode should keep the full span from p24 start to p24_C end."""

    class_name, components, note = _multimer_components_for_header(
        header="class2",
        sequence="ABCDEFGHIJKLMNOP",
        domain_hits={
            "Gag_p10": _DomainHit("Gag_p10", 0, 1),
            "Gag_p24": _DomainHit("Gag_p24", 4, 6),
            "Gag_p24_C": _DomainHit("Gag_p24_C", 10, 12),
        },
        combine_class2_capsid=True,
    )

    assert class_name == "class_2"
    assert components == [("matrix", "AB", 3), ("capsid", "EFGHIJKLM", 6)]
    assert "including linker" in note


def test_build_boltz_command_supports_lmi4boltz(tmp_path: Path) -> None:
    """The multimer runner should be able to invoke the LMI fork in-place."""

    args = Namespace(
        run_root=tmp_path / "runs",
        config_name="cfg",
        cache="~/.boltz",
        devices=1,
        accelerator="gpu",
        recycling_steps=5,
        sampling_steps=40,
        diffusion_samples=6,
        max_parallel_samples=1,
        num_workers=2,
        preprocessing_threads=8,
        max_msa_seqs=8192,
        model="boltz2",
        use_msa_server=True,
        msa_server_url="https://api.colabfold.com",
        msa_pairing_strategy="greedy",
        seed=None,
        chunk_size_transition_z=32,
        chunk_size_transition_msa=16,
        chunk_size_outer_product=1,
        chunk_size_tri_attn=64,
        triangle_mult_gate_nchunks=4,
        chunk_size_threshold=1,
        use_bfloat16=False,
        use_lmi4boltz=True,
        python_executable=tmp_path / "venv" / "bin" / "python",
        lmi4boltz_src=tmp_path / "lmi4boltz" / "src",
        mplconfigdir=tmp_path / "mpl",
    )
    target = _MODULE.MultimerTarget(
        index=1,
        source_header="header",
        target_id="target",
        input_yaml=tmp_path / "target.yaml",
        class_name="class_2",
        assembly_label="3x10",
        total_length=30,
        note="note",
    )

    command = _build_boltz_command(args, target, override=True)
    env = _build_boltz_env(args)

    assert command[:4] == [str(args.python_executable), "-m", "boltz.main", "predict"]
    assert "--chunk_size_transition_z" in command
    assert "--chunk_size_tri_attn" in command
    assert command[-1] == "--override"
    assert env["PYTHONPATH"].startswith(str(args.lmi4boltz_src.resolve()))
    assert env["MPLCONFIGDIR"] == str(args.mplconfigdir)
