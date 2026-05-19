# Scientific Report: Boltz2 Experiments on Influenza HA Ligand Binding

## Objective

We evaluated whether local `boltz2` inference on a 32 GB GPU can model influenza hemagglutinin (HA) binding to glycans and recover expected relative binding preferences. The main biological question was whether the model distinguishes `alpha2-3` versus `alpha2-6` sialylated glycans in a way consistent with known host-adaptation biology.

## Methods

All experiments were run locally with `boltz2` under BF16 mixed precision. We first established a conservative inference configuration and then performed a small hyperparameter search over:

- `diffusion_samples`
- `max_parallel_samples`
- `sampling_steps`
- `recycling_steps`

We compared:

1. Simple CCD ligands: `SIA` versus `NGC`
2. Full glycan ligands encoded as single SMILES molecules, derived from experimentally resolved `alpha2-3` and `alpha2-6` sialylated glycans

For affinity-enabled runs, Boltz reports:

- `affinity_pred_value` (interpreted as `log10(IC50 in uM)`, lower is stronger)
- `affinity_probability_binary` (binder probability)
- structural confidence terms such as `complex_plddt` and `ligand_iptm`

## Results

### 1. Feasibility

Local inference was successful on 32 GB VRAM. A minimal flu HA plus ligand case ran comfortably, and even `diffusion_samples=6` with moderate refinement was practical.

### 2. Hyperparameter sensitivity

Prediction quality was highly sensitive to very cheap settings. A single-sample, low-step run produced poor confidence. Increasing to multiple diffusion samples and at least moderate sampling steps dramatically improved structural confidence.

However, beyond approximately:

- `diffusion_samples=6`
- `sampling_steps=40`
- `recycling_steps=3`

returns largely plateaued. Increasing to `sampling_steps=200` did not materially change the conclusions. A practical recommended setting is therefore:

`--diffusion_samples 6 --max_parallel_samples 2 --sampling_steps 40 --recycling_steps 3`

### 3. SIA versus NGC

Across all tested settings, Boltz consistently preferred `SIA` over `NGC`.

At the highest tested setting:

- `SIA`: `affinity_pred_value ~= 0.538`, binder probability `~= 0.606`
- `NGC`: `affinity_pred_value ~= 0.960`, binder probability `~= 0.505`

This ranking was stable across the inference ladder. Protein confidence remained high and nearly unchanged across runs, with protein/complex pLDDT generally around `0.89-0.90`.

### 4. alpha2-3 versus alpha2-6 glycans

We extracted isomeric SMILES from two SDF structures representing `alpha2-3` and `alpha2-6` sialylated glycans and ran affinity predictions using those full glycans as single ligands.

The model preferred `alpha2-3` over `alpha2-6`:

- `alpha2-3`: `affinity_pred_value ~= 0.912`, binder probability `~= 0.234`
- `alpha2-6`: `affinity_pred_value ~= 1.641`, binder probability `~= 0.224`

This is opposite to the expected biological preference for the system under study.

### 5. Effect of MSA depth

We also tested the Memphis `HA + SIA` example in true single-sequence mode and
across a shallow-MSA ladder using:

- `msa: empty` for the no-MSA control
- `--subsample_msa --num_subsampled_msa N --max_msa_seqs N`

while keeping inference fixed at:

`--diffusion_samples 6 --max_parallel_samples 2 --sampling_steps 40 --recycling_steps 3`

Results:

| Condition | affinity_pred_value | binder probability | complex_plddt | ligand_iptm |
| --- | ---: | ---: | ---: | ---: |
| `no_msa` | `1.294` | `0.456` | `0.856` | `0.974` |
| `msa_16` | `1.237` | `0.387` | `0.840` | `0.967` |
| `msa_32` | `1.061` | `0.401` | `0.832` | `0.941` |
| `msa_64` | `0.873` | `0.500` | `0.848` | `0.981` |
| `msa_128` | `0.695` | `0.559` | `0.880` | `0.981` |
| `msa_256` | `0.597` | `0.588` | `0.885` | `0.985` |
| `msa_512` | `0.602` | `0.603` | `0.883` | `0.984` |
| `msa_1024` | `0.557` | `0.606` | `0.889` | `0.984` |

The main pattern is that no-MSA mode substantially degrades the prediction,
while increasing MSA depth steadily improves the result up to roughly
`128-256` sequences. Beyond that, gains largely plateau.

## Interpretation

The `SIA > NGC` result appears internally robust within this pipeline. In contrast, the `alpha2-3 > alpha2-6` result should not be treated as biologically trustworthy.

The likely reasons are methodological:

1. The local YAML/SMILES path does not preserve the original experimental 3D glycan geometry. Ligands are rebuilt from SMILES.
2. The affinity head appears more naturally suited to small-molecule-like ligands than large flexible glycans.
3. Both glycan binder probabilities were low in absolute terms, indicating a weak-signal regime.
4. MSA depth materially affects quality; single-sequence mode is clearly weaker than even moderate-depth MSA inputs.

## Conclusions

1. Boltz2 runs reliably on local 32 GB hardware for these HA-ligand experiments.
2. Moderate inference settings are sufficient; more expensive settings do not materially improve conclusions.
3. The local pipeline gives a stable `SIA > NGC` ranking.
4. The Memphis `SIA` example is highly MSA-dependent; most of the recoverable signal is already present by roughly `128-256` MSA sequences, with only modest gains beyond that.
5. The current SMILES-based glycan workflow does not provide a trustworthy test of the biologically important `alpha2-3` versus `alpha2-6` question.
6. A more faithful next step would be to preserve the glycan's original 3D geometry, likely by adding native `sdf:` ligand support or by evaluating structure/interface geometry rather than relying on the current affinity score alone.
