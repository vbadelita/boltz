# Test EVE Multimer Experiment

## Overview

- Targets: 9
- Completed targets: 3
- Config: `multimer_s40_n1_p1_r3`
- Sampling steps: 40
- Diffusion samples: 1
- Max parallel samples: 1
- Recycling steps: 3

## Best Result Per Multimer

| Target | Kind | Stoichiometry | Total residues | MSA raw | Status | Best rank | complex_plddt | ptm | confidence_score |
|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| shEPLVP1_SHeplvp1 | sheplvp1 | 3 | 477 | 408 | ok | 0 | 0.673 | 0.687 | 0.659 |
| DBA62700.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 3] | gag | 6 | 3546 | 3355 | oom | - | - | - | - |
| DBA62703.1 TPA_asm: Gag [Marsupial endogenous betaretrovirus 5] | gag | 6 | 3354 | 5108 | oom | - | - | - | - |
| DBA62706.1 TPA_asm: Env [Marsupial endogenous betaretrovirus 5] | env | 3 | 1467 | 230 | ok | 0 | 0.539 | 0.586 | 0.546 |
| sp\|Q9YNA8.3\|GAK19_HUMAN RecName: Full=Endogenous retrovirus group K member 19 Gag polyprotein; AltName: Full=HERV-K(C19) Gag protein; AltName: Full=HERV-K_19q11 provirus ancestral Gag polyprotein; Short=Gag polyprotein | gag | 6 | 3996 | 5214 | oom | - | - | - | - |
| AAC82593.1 Gag [Human immunodeficiency virus 1] | gag | 6 | 3000 | 3834 | oom | - | - | - | - |
| AAC82596.1 Env [Human immunodeficiency virus 1] | env | 3 | 2568 | 10996 | oom | - | - | - | - |
| AAC82567.1 gPr80 [Moloney murine leukemia virus] | env | 3 | 1995 | 2990 | oom | - | - | - | - |
| AAC82588.1 gp60 SU [Bovine leukemia virus] | env | 3 | 1545 | 499 | ok | 0 | 0.362 | 0.310 | 0.351 |

## Files

- Metadata CSV: `test_eve_experiments/multimers/metadata.csv`
- Structures: `test_eve_experiments/multimers/structures/`
- Raw outputs: `test_eve_experiments/multimers/runs/`
- Logs: `test_eve_experiments/multimers/logs/`
