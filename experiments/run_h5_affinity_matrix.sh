#!/usr/bin/env bash
set -euo pipefail

# Reproduce the 2x2 H5 affinity experiment matrix from the repository root.

COMMON_ARGS=(
  --model boltz2
  --accelerator gpu
  --use_msa_server
  --diffusion_samples 6
  --max_parallel_samples 2
  --sampling_steps 40
  --recycling_steps 3
  --override
)

uv run boltz predict examples/h5_texas_sia_affinity.yaml \
  "${COMMON_ARGS[@]}" \
  --out_dir ./runs/h5_texas_sia_affinity_s40_n6_p2_r3

uv run boltz predict examples/h5_texas_ngc_affinity.yaml \
  "${COMMON_ARGS[@]}" \
  --out_dir ./runs/h5_texas_ngc_affinity_s40_n6_p2_r3

uv run boltz predict examples/h5_california_sia_affinity.yaml \
  "${COMMON_ARGS[@]}" \
  --out_dir ./runs/h5_california_sia_affinity_s40_n6_p2_r3

uv run boltz predict examples/h5_california_ngc_affinity.yaml \
  "${COMMON_ARGS[@]}" \
  --out_dir ./runs/h5_california_ngc_affinity_s40_n6_p2_r3
