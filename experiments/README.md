# Experiments

This folder is the entrypoint for reproducing the local Boltz influenza binding
analysis after a long gap.

## Start Here

Read these in order:

1. [Scientific report](../docs/flu_binding_report.md)
2. This file
3. The relevant `runs/<name>/repro.sh` for any run you want to inspect or rerun

Core locations:

- Inputs: `examples/`
- Saved outputs: `runs/`
- Reproducible run scripts: `experiments/*.sh`
- Analysis summaries: `experiments/*.py`

Current stable inference setting used for most comparisons:

```text
--diffusion_samples 6
--max_parallel_samples 2
--sampling_steps 40
--recycling_steps 3
```

## Main Findings

- Boltz2 runs locally on the 32 GB GPU for these HA-ligand experiments.
- For the Memphis HA case, the model consistently prefers `SIA` over `NGC`.
- For the H5 two-sequence comparison, Boltz still prefers `SIA` for both
  sequences and does not recover the expected ligand preference switch.
- For the `alpha2-3` versus `alpha2-6` glycan test, the current SMILES-based
  workflow gives the biologically wrong direction and should not be trusted.
- MSA depth matters a lot. No-MSA mode is substantially worse, while most gains
  are recovered by roughly `128-256` sequences.

## Known Limitations

- The local ligand YAML path supports `ccd:` and `smiles:`, but not native
  `sdf:` input.
- SMILES input does not preserve the original 3D ligand geometry.
- The affinity score is more believable for simple ligands than for large,
  flexible glycans.
- `run_memphis_msa_ladder.sh` reruns the full ladder if invoked.

## H5 Affinity Matrix

Run the four H5 protein/ligand combinations:

```bash
./experiments/run_h5_affinity_matrix.sh
```

This reproduces:

- `Texas x SIA`
- `Texas x NGC`
- `California x SIA`
- `California x NGC`

using the shared setting:

```text
diffusion_samples=6
max_parallel_samples=2
sampling_steps=40
recycling_steps=3
```

Summarize those outputs:

```bash
uv run python experiments/analyze_h5_affinity.py
```

Expected interpretation:

- Both H5 proteins come out `SIA`-favored.
- The model does not capture the expected switch between the two sequences.

## Flu Affinity Summaries

Summarize the earlier flu runs, including the hyperparameter ladder and the
SMILES-based `alpha2-3` versus `alpha2-6` comparison:

```bash
uv run python experiments/analyze_flu_affinity.py
```

Expected interpretation:

- `SIA > NGC` is stable across the tested inference settings.
- The `alpha2-3 > alpha2-6` outcome is not biologically trustworthy in the
  current workflow.

## Memphis MSA Ablation

Run the Memphis `SIA` no-MSA control plus the shallow-MSA ladder:

```bash
./experiments/run_memphis_msa_ladder.sh
```

This includes:

- `no_msa`
- `msa_16`
- `msa_32`
- `msa_64`
- `msa_128`
- `msa_256`
- `msa_512`
- `msa_1024`

Summarize those outputs:

```bash
uv run python experiments/analyze_memphis_msa_ladder.py
```

Expected interpretation:

- `no_msa` is clearly worse than MSA-supported runs.
- Performance improves strongly with depth up to roughly `128-256`.
- `256`, `512`, and `1024` are in the plateau regime.

## Notes

- All Boltz outputs live under `runs/`.
- Each run directory should contain a `repro.sh` file with the exact command.
- The analysis scripts read the saved JSON outputs directly and do not rerun
  inference.
- If you want to rerun only one condition, prefer the individual
  `runs/<name>/repro.sh` script over the batch runners in this folder.
