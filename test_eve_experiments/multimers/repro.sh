#!/usr/bin/env bash
set -euo pipefail

cd /home/vlad/code/boltz
uv run python experiments/run_test_eve_multimers.py --fasta /home/vlad/code/boltz/test_EVE_structures.fasta --out-dir /home/vlad/code/boltz/test_eve_experiments/multimers --boltz-cmd 'uv run boltz' --accelerator gpu --model boltz2 --sampling-steps 40 --diffusion-samples 1 --max-parallel-samples 1 --recycling-steps 3
