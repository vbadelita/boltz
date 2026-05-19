from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


FLU_LADDER = [
    (
        "20/r0",
        "runs/flu_sia_affinity_s20_n6_p2_r0/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_ngc_affinity_s20_n6_p2_r0/boltz_results_flu_ngc_affinity/predictions/flu_ngc_affinity/affinity_flu_ngc_affinity.json",
    ),
    (
        "40/r0",
        "runs/flu_sia_affinity_s40_n6_p2_r0/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_ngc_affinity_s40_n6_p2_r0/boltz_results_flu_ngc_affinity/predictions/flu_ngc_affinity/affinity_flu_ngc_affinity.json",
    ),
    (
        "40/r3",
        "runs/flu_sia_affinity_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_ngc_affinity_s40_n6_p2_r3/boltz_results_flu_ngc_affinity/predictions/flu_ngc_affinity/affinity_flu_ngc_affinity.json",
    ),
    (
        "200/r3",
        "runs/flu_sia_affinity_s200_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_ngc_affinity_s200_n6_p2_r3/boltz_results_flu_ngc_affinity/predictions/flu_ngc_affinity/affinity_flu_ngc_affinity.json",
    ),
]


GLYCAN_COMPARISON = [
    (
        "a2-3",
        "runs/flu_sia_a23_smiles_affinity_s40_n6_p2_r3/boltz_results_flu_sia_a23_smiles_affinity/predictions/flu_sia_a23_smiles_affinity/affinity_flu_sia_a23_smiles_affinity.json",
        "runs/flu_sia_a23_smiles_affinity_s40_n6_p2_r3/boltz_results_flu_sia_a23_smiles_affinity/predictions/flu_sia_a23_smiles_affinity/confidence_flu_sia_a23_smiles_affinity_model_0.json",
    ),
    (
        "a2-6",
        "runs/flu_sia_a26_smiles_affinity_s40_n6_p2_r3/boltz_results_flu_sia_a26_smiles_affinity/predictions/flu_sia_a26_smiles_affinity/affinity_flu_sia_a26_smiles_affinity.json",
        "runs/flu_sia_a26_smiles_affinity_s40_n6_p2_r3/boltz_results_flu_sia_a26_smiles_affinity/predictions/flu_sia_a26_smiles_affinity/confidence_flu_sia_a26_smiles_affinity_model_0.json",
    ),
]


def load_json(path: str) -> dict:
    return json.loads((ROOT / path).read_text())


def print_ladder() -> None:
    print("SIA vs NGC ladder")
    print("setting,sia_value,sia_prob,ngc_value,ngc_prob,preferred")
    for label, sia_path, ngc_path in FLU_LADDER:
        sia = load_json(sia_path)
        ngc = load_json(ngc_path)
        preferred = "SIA" if sia["affinity_pred_value"] < ngc["affinity_pred_value"] else "NGC"
        print(
            f"{label},"
            f"{sia['affinity_pred_value']:.6f},{sia['affinity_probability_binary']:.6f},"
            f"{ngc['affinity_pred_value']:.6f},{ngc['affinity_probability_binary']:.6f},"
            f"{preferred}"
        )


def print_glycan_comparison() -> None:
    print("\nGlycan comparison")
    print("ligand,affinity_pred_value,affinity_probability_binary,complex_plddt,ligand_iptm")
    for label, aff_path, conf_path in GLYCAN_COMPARISON:
        affinity = load_json(aff_path)
        confidence = load_json(conf_path)
        print(
            f"{label},"
            f"{affinity['affinity_pred_value']:.6f},"
            f"{affinity['affinity_probability_binary']:.6f},"
            f"{confidence['complex_plddt']:.6f},"
            f"{confidence['ligand_iptm']:.6f}"
        )


def main() -> None:
    print_ladder()
    print_glycan_comparison()


if __name__ == "__main__":
    main()
