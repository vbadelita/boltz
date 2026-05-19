#!/usr/bin/env bash
set -euo pipefail

# Reproduce the Memphis SIA MSA-ablation experiment from the repository root.

COMMON_ARGS=(
  --model boltz2
  --accelerator gpu
  --diffusion_samples 6
  --max_parallel_samples 2
  --sampling_steps 40
  --recycling_steps 3
  --override
)

uv run boltz predict examples/flu_sia_affinity_no_msa.yaml \
  "${COMMON_ARGS[@]}" \
  --out_dir ./runs/flu_sia_affinity_no_msa_s40_n6_p2_r3

for n in 16 32 64 128 256 512 1024; do
  uv run boltz predict examples/flu_sia_affinity.yaml \
    "${COMMON_ARGS[@]}" \
    --use_msa_server \
    --subsample_msa \
    --num_subsampled_msa "${n}" \
    --max_msa_seqs "${n}" \
    --out_dir "./runs/flu_sia_affinity_msa${n}_s40_n6_p2_r3"
done
