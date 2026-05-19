# Test EVE Boltz Experiment

## Overview

- Input FASTA: `test_EVE_structures.fasta`
- Targets: 12
- Folded at complex_plddt >= 0.70: 4
- Low confidence, 0.50 <= complex_plddt < 0.70: 6
- Failed or very low confidence: 2

The starting point was the local H5 setting: `diffusion_samples=6`, `max_parallel_samples=2`, `sampling_steps=40`, `recycling_steps=3`. The grid then increased sampling steps, diffusion samples, or recycling for the lowest-confidence baseline targets.

## Best Result Per Target

| # | Target | Length | MSA used | MSA raw | Best status | Best config | Best sample rank | complex_plddt | ptm | confidence_score |
|---:|---|---:|---:|---:|---|---|---:|---:|---:|---:|
| 1 | AAD14589.1 polymerase [Ebola virus - Mayinga, Zaire, 1976] | 1157 | 268 | 271 | folded | h5_start_s40_n6_p2_r3 | 0 | 0.876 | 0.926 | 0.886 |
| 2 | shEPLVP1_SHeplvp1 | 159 | 355 | 408 | folded | all_s200_n6_p2_r5 | 0 | 0.784 | 0.659 | 0.759 |
| 3 | AAA67459.1 polyprotein [Canine parvovirus] | 668 | 542 | 561 | folded | h5_start_s40_n6_p2_r3 | 0 | 0.754 | 0.505 | 0.705 |
| 4 | DBA62700.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 3] | 591 | 3328 | 3355 | low confidence | all_s200_n6_p2_r5 | 2 | 0.543 | 0.328 | 0.500 |
| 5 | DBA62703.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 5] | 559 | 5037 | 5108 | low confidence | more_samples_s40_n8_p2_r3 | 1 | 0.593 | 0.335 | 0.542 |
| 6 | DBA62706.1 TPA_asm: Env [Marsupial endogenous betaretrovirus 5] | 489 | 227 | 230 | failed | more_samples_s40_n8_p2_r3 | 0 | 0.486 | 0.511 | 0.491 |
| 7 | sp\|Q9YNA8.3\|GAK19_HUMAN RecName: Full=Endogenous retrovirus group K member 19 Gag polyprotein; AltName: Full=HERV-K(C19) Gag protein; AltName: Full=HERV-K_19q11 provirus ancestral Gag polyprotein; Short=Gag polyprotein | 666 | 5062 | 5214 | low confidence | more_samples_s40_n8_p2_r3 | 0 | 0.512 | 0.327 | 0.475 |
| 8 | AAC82593.1 Gag [Human immunodeficiency virus 1] | 500 | 3791 | 3834 | low confidence | more_recycling_s40_n6_p2_r5 | 0 | 0.654 | 0.482 | 0.620 |
| 9 | AAC82596.1 Env [Human immunodeficiency virus 1] | 856 | 8192 | 10996 | low confidence | more_samples_s40_n8_p2_r3 | 0 | 0.599 | 0.593 | 0.598 |
| 10 | AAC82567.1 gPr80 [Moloney murine leukemia virus] | 665 | 2967 | 2990 | low confidence | all_s200_n6_p2_r5 | 0 | 0.694 | 0.487 | 0.652 |
| 11 | AAC82588.1 gp60 SU [Bovine leukemia virus] | 515 | 476 | 499 | failed | more_recycling_s40_n6_p2_r5 | 0 | 0.366 | 0.264 | 0.346 |
| 12 | NP_955618.1 p24 CA [Human T-cell leukemia virus type I] | 215 | 536 | 587 | folded | h5_start_s40_n6_p2_r3 | 0 | 0.863 | 0.582 | 0.807 |

## Grid Search Targets

| Target | Baseline pLDDT | Best grid pLDDT | Best grid config | Delta |
|---|---:|---:|---|---:|
| DBA62700.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 3] | 0.542 | 0.539 | more_samples_s40_n8_p2_r3 | -0.002 |
| DBA62703.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 5] | 0.590 | 0.593 | more_samples_s40_n8_p2_r3 | 0.004 |
| DBA62706.1 TPA_asm: Env [Marsupial endogenous betaretrovirus 5] | 0.486 | 0.486 | more_samples_s40_n8_p2_r3 | 0.000 |
| sp\|Q9YNA8.3\|GAK19_HUMAN RecName: Full=Endogenous retrovirus group K member 19 Gag polyprotein; AltName: Full=HERV-K(C19) Gag protein; AltName: Full=HERV-K_19q11 provirus ancestral Gag polyprotein; Short=Gag polyprotein | 0.509 | 0.512 | more_samples_s40_n8_p2_r3 | 0.003 |
| AAC82593.1 Gag [Human immunodeficiency virus 1] | 0.652 | 0.654 | more_recycling_s40_n6_p2_r5 | 0.002 |
| AAC82596.1 Env [Human immunodeficiency virus 1] | 0.594 | 0.599 | more_samples_s40_n8_p2_r3 | 0.005 |
| AAC82567.1 gPr80 [Moloney murine leukemia virus] | 0.689 | 0.693 | more_recycling_s40_n6_p2_r5 | 0.005 |
| AAC82588.1 gp60 SU [Bovine leukemia virus] | 0.357 | 0.366 | more_recycling_s40_n6_p2_r5 | 0.009 |

## Failed-Target Retry

The very low-confidence targets were retried with `sampling_steps=200`, `diffusion_samples=6`, `max_parallel_samples=2`, `recycling_steps=5`.

| Target | Previous best pLDDT | Retry best pLDDT | Delta |
|---|---:|---:|---:|
| DBA62706.1 TPA_asm: Env [Marsupial endogenous betaretrovirus 5] | 0.486 | 0.452 | -0.033 |
| AAC82588.1 gp60 SU [Bovine leukemia virus] | 0.366 | 0.358 | -0.008 |

## 200-Step Check

All targets were checked with `sampling_steps=200`, `diffusion_samples=6`, `max_parallel_samples=2`, `recycling_steps=5`. The two very low-confidence targets use the `retry_s200_n6_p2_r5` outputs above; the rest use `all_s200_n6_p2_r5`.

| Target | Previous best pLDDT | Best 200-step pLDDT | 200-step config | Delta |
|---|---:|---:|---|---:|
| AAD14589.1 polymerase [Ebola virus - Mayinga, Zaire, 1976] | 0.876 | 0.874 | all_s200_n6_p2_r5 | -0.003 |
| shEPLVP1_SHeplvp1 | 0.784 | 0.784 | all_s200_n6_p2_r5 | 0.000 |
| AAA67459.1 polyprotein [Canine parvovirus] | 0.754 | 0.752 | all_s200_n6_p2_r5 | -0.003 |
| DBA62700.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 3] | 0.542 | 0.543 | all_s200_n6_p2_r5 | 0.001 |
| DBA62703.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 5] | 0.593 | 0.588 | all_s200_n6_p2_r5 | -0.006 |
| DBA62706.1 TPA_asm: Env [Marsupial endogenous betaretrovirus 5] | 0.486 | 0.452 | retry_s200_n6_p2_r5 | -0.033 |
| sp\|Q9YNA8.3\|GAK19_HUMAN RecName: Full=Endogenous retrovirus group K member 19 Gag polyprotein; AltName: Full=HERV-K(C19) Gag protein; AltName: Full=HERV-K_19q11 provirus ancestral Gag polyprotein; Short=Gag polyprotein | 0.512 | 0.502 | all_s200_n6_p2_r5 | -0.010 |
| AAC82593.1 Gag [Human immunodeficiency virus 1] | 0.654 | 0.639 | all_s200_n6_p2_r5 | -0.016 |
| AAC82596.1 Env [Human immunodeficiency virus 1] | 0.599 | 0.591 | all_s200_n6_p2_r5 | -0.008 |
| AAC82567.1 gPr80 [Moloney murine leukemia virus] | 0.693 | 0.694 | all_s200_n6_p2_r5 | 0.000 |
| AAC82588.1 gp60 SU [Bovine leukemia virus] | 0.366 | 0.358 | retry_s200_n6_p2_r5 | -0.008 |
| NP_955618.1 p24 CA [Human T-cell leukemia virus type I] | 0.863 | 0.854 | all_s200_n6_p2_r5 | -0.009 |

## Files

- Metadata CSV: `test_eve_experiments/metadata.csv`
- Collected structures: `test_eve_experiments/structures/`
- Raw Boltz outputs: `test_eve_experiments/runs/`
- Logs: `test_eve_experiments/logs/`
