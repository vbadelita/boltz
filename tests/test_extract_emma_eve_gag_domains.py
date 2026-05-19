from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

_MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "process"
    / "extract_emma_eve_gag_domains.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "extract_emma_eve_gag_domains",
    _MODULE_PATH,
)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)

_build_output_records = _MODULE.build_output_records
_collect_domain_hits = _MODULE.collect_domain_hits
_write_output_fastas = _MODULE.write_output_fastas


def _make_align(start: int, end: int, domain_name: str) -> dict:
    length = end - start + 1
    return {
        "segs": {
            "denseg": {
                "dim": 2,
                "numseg": 1,
                "starts": [start, 0],
                "lens": [length],
            }
        },
        "id": [
            {"str": domain_name},
            {"str": f"{domain_name} superfamily"},
        ],
    }


def test_collect_domain_hits_maps_class_specific_domains(tmp_path: Path) -> None:
    """Matrix and capsid outputs should map from the expected CD-search labels."""

    aligndata = tmp_path / "aligndata.json"
    records = [
        {
            "Seq_annot": {
                "desc": [{"title": "Q#1 - >class2"}],
                "data": {
                    "align": [
                        _make_align(0, 2, "Gag_p10"),
                        _make_align(3, 5, "Gag_p24"),
                        _make_align(6, 8, "Gag_p24_C"),
                    ]
                },
            }
        },
        {
            "Seq_annot": {
                "desc": [{"title": "Q#2 - >class1"}],
                "data": {"align": [_make_align(1, 4, "Gag_MA")]},
            }
        },
    ]
    aligndata.write_text("".join(json.dumps(record) for record in records))

    hits_by_header = _collect_domain_hits(aligndata)

    assert hits_by_header["class2"]["matrix_protein"].domain_name == "Gag_p10"
    assert hits_by_header["class2"]["capsid_n_terminal"].start == 3
    assert hits_by_header["class2"]["capsid_c_terminal"].end == 8
    assert hits_by_header["class1"]["matrix_protein"].domain_name == "Gag_MA"


def test_build_output_records_slices_expected_subsequences() -> None:
    """Extracted records should use the annotated residue spans."""

    outputs = _build_output_records(
        fasta_records=[("seq1", "ABCDEFGHI"), ("seq2", "MNOPQRST")],
        hits_by_header={
            "seq1": {
                "matrix_protein": _MODULE.DomainHit("Gag_p10", 0, 2),
                "capsid_n_terminal": _MODULE.DomainHit("Gag_p24", 3, 5),
                "capsid_c_terminal": _MODULE.DomainHit("Gag_p24_C", 6, 8),
            },
            "seq2": {
                "matrix_protein": _MODULE.DomainHit("Gag_MA", 1, 4),
            },
        },
    )

    assert [sequence for _, sequence in outputs["matrix_protein"]] == ["ABC", "NOPQ"]
    assert [sequence for _, sequence in outputs["capsid_n_terminal"]] == ["DEF"]
    assert [sequence for _, sequence in outputs["capsid_c_terminal"]] == ["GHI"]


def test_write_output_fastas_creates_expected_files(tmp_path: Path) -> None:
    """Each output category should be written as its own FASTA file."""

    written_paths = _write_output_fastas(
        {
            "matrix_protein": [("seq1 | matrix", "ABC")],
            "capsid_n_terminal": [("seq1 | capsid_n", "DEF")],
            "capsid_c_terminal": [],
        },
        tmp_path,
    )

    assert written_paths["matrix_protein"].read_text() == ">seq1 | matrix\nABC\n"
    assert written_paths["capsid_n_terminal"].read_text() == ">seq1 | capsid_n\nDEF\n"
    assert written_paths["capsid_c_terminal"].read_text() == ""
