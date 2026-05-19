#!/usr/bin/env bash
set -euo pipefail

cd /home/vlad/code/boltz
uv run python experiments/run_test_eve_experiments.py --fasta /home/vlad/code/boltz/test_EVE_structures.fasta --out-dir /home/vlad/code/boltz/test_eve_experiments --boltz-cmd 'uv run boltz' --accelerator gpu --model boltz2 --low-plddt 0.7 --grid-pilot-size 3
