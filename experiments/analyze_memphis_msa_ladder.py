from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


RUNS = [
    (
        "no_msa",
        "runs/flu_sia_affinity_no_msa_s40_n6_p2_r3/boltz_results_flu_sia_affinity_no_msa/predictions/flu_sia_affinity_no_msa/affinity_flu_sia_affinity_no_msa.json",
        "runs/flu_sia_affinity_no_msa_s40_n6_p2_r3/boltz_results_flu_sia_affinity_no_msa/predictions/flu_sia_affinity_no_msa/confidence_flu_sia_affinity_no_msa_model_0.json",
    ),
    (
        "msa_16",
        "runs/flu_sia_affinity_msa16_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_sia_affinity_msa16_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/confidence_flu_sia_affinity_model_0.json",
    ),
    (
        "msa_32",
        "runs/flu_sia_affinity_msa32_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_sia_affinity_msa32_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/confidence_flu_sia_affinity_model_0.json",
    ),
    (
        "msa_64",
        "runs/flu_sia_affinity_msa64_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_sia_affinity_msa64_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/confidence_flu_sia_affinity_model_0.json",
    ),
    (
        "msa_128",
        "runs/flu_sia_affinity_msa128_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_sia_affinity_msa128_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/confidence_flu_sia_affinity_model_0.json",
    ),
    (
        "msa_256",
        "runs/flu_sia_affinity_msa256_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_sia_affinity_msa256_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/confidence_flu_sia_affinity_model_0.json",
    ),
    (
        "msa_512",
        "runs/flu_sia_affinity_msa512_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_sia_affinity_msa512_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/confidence_flu_sia_affinity_model_0.json",
    ),
    (
        "msa_1024",
        "runs/flu_sia_affinity_msa1024_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/affinity_flu_sia_affinity.json",
        "runs/flu_sia_affinity_msa1024_s40_n6_p2_r3/boltz_results_flu_sia_affinity/predictions/flu_sia_affinity/confidence_flu_sia_affinity_model_0.json",
    ),
]


def load_json(path: str) -> dict:
    return json.loads((ROOT / path).read_text())


def main() -> None:
    print("condition,affinity_pred_value,affinity_probability_binary,complex_plddt,ligand_iptm")
    for label, aff_path, conf_path in RUNS:
        aff = load_json(aff_path)
        conf = load_json(conf_path)
        print(
            f"{label},"
            f"{aff['affinity_pred_value']:.6f},"
            f"{aff['affinity_probability_binary']:.6f},"
            f"{conf['complex_plddt']:.6f},"
            f"{conf['ligand_iptm']:.6f}"
        )


if __name__ == "__main__":
    main()
