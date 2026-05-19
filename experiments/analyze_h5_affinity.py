from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


RUNS = {
    "Texas + SIA": (
        "runs/h5_texas_sia_affinity_s40_n6_p2_r3/boltz_results_h5_texas_sia_affinity/predictions/h5_texas_sia_affinity/affinity_h5_texas_sia_affinity.json",
        "runs/h5_texas_sia_affinity_s40_n6_p2_r3/boltz_results_h5_texas_sia_affinity/predictions/h5_texas_sia_affinity/confidence_h5_texas_sia_affinity_model_0.json",
    ),
    "Texas + NGC": (
        "runs/h5_texas_ngc_affinity_s40_n6_p2_r3/boltz_results_h5_texas_ngc_affinity/predictions/h5_texas_ngc_affinity/affinity_h5_texas_ngc_affinity.json",
        "runs/h5_texas_ngc_affinity_s40_n6_p2_r3/boltz_results_h5_texas_ngc_affinity/predictions/h5_texas_ngc_affinity/confidence_h5_texas_ngc_affinity_model_0.json",
    ),
    "California + SIA": (
        "runs/h5_california_sia_affinity_s40_n6_p2_r3/boltz_results_h5_california_sia_affinity/predictions/h5_california_sia_affinity/affinity_h5_california_sia_affinity.json",
        "runs/h5_california_sia_affinity_s40_n6_p2_r3/boltz_results_h5_california_sia_affinity/predictions/h5_california_sia_affinity/confidence_h5_california_sia_affinity_model_0.json",
    ),
    "California + NGC": (
        "runs/h5_california_ngc_affinity_s40_n6_p2_r3/boltz_results_h5_california_ngc_affinity/predictions/h5_california_ngc_affinity/affinity_h5_california_ngc_affinity.json",
        "runs/h5_california_ngc_affinity_s40_n6_p2_r3/boltz_results_h5_california_ngc_affinity/predictions/h5_california_ngc_affinity/confidence_h5_california_ngc_affinity_model_0.json",
    ),
}


def load_json(path: str) -> dict:
    return json.loads((ROOT / path).read_text())


def main() -> None:
    rows: list[tuple[str, str, float, float, float, float]] = []
    for label, (aff_path, conf_path) in RUNS.items():
        affinity = load_json(aff_path)
        confidence = load_json(conf_path)
        protein, ligand = label.split(" + ")
        rows.append(
            (
                protein,
                ligand,
                affinity["affinity_pred_value"],
                affinity["affinity_probability_binary"],
                confidence["complex_plddt"],
                confidence["ligand_iptm"],
            )
        )

    print(
        "protein,ligand,affinity_pred_value,affinity_probability_binary,complex_plddt,ligand_iptm"
    )
    for row in rows:
        protein, ligand, value, prob, plddt, iptm = row
        print(f"{protein},{ligand},{value:.6f},{prob:.6f},{plddt:.6f},{iptm:.6f}")

    by_protein: dict[str, dict[str, tuple[float, float]]] = {}
    for protein, ligand, value, prob, _, _ in rows:
        by_protein.setdefault(protein, {})[ligand] = (value, prob)

    print("\nPreference summary")
    for protein, values in by_protein.items():
        sia_value, sia_prob = values["SIA"]
        ngc_value, ngc_prob = values["NGC"]
        preferred = "SIA" if sia_value < ngc_value else "NGC"
        print(
            f"{protein}: preferred={preferred}, "
            f"delta_value={sia_value - ngc_value:.6f}, "
            f"delta_prob={sia_prob - ngc_prob:.6f}"
        )


if __name__ == "__main__":
    main()
