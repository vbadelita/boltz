# Boltz-1 Training Resources Summary

This note summarizes what the Boltz-1 paper and the public repo explicitly disclose about the resources used to train Boltz-1, and what they do not disclose.

## What is explicitly disclosed

Source: `boltz1-paper.pdf`

- **Structure training length**
  - `68,000` training steps total
- **Global batch size**
  - `128`
- **Training stages**
  - First `53,000` steps:
    - crop size `384` tokens
    - `3456` atoms
    - sampled equally from:
      - PDB training data
      - OpenFold distillation data
  - Last `15,000` steps:
    - crop size `512` tokens
    - `4608` atoms
    - sampled from:
      - PDB only
- **OpenFold distillation dataset size**
  - approximately `270K` structures

## Simple derived quantity

Assuming the reported batch size is the effective global batch size:

- total cropped training samples seen:
  - `68,000 * 128 = 8,704,000`

This is a count of training samples/steps, not unique structures.

## What the paper says about compute

The paper does not provide an exact GPU count or wall-clock duration for Boltz-1 training.

It does state:

- AlphaFold3 trained a similar architecture for nearly:
  - `150,000` steps
  - batch size `256`
- and that this required approximately:
  - `4x` the computing time

So the paper's claim is that Boltz-1 required substantially less compute than the AlphaFold3 setup they compare against, but it does **not** quantify Boltz-1 training time directly.

## Hardware/resources explicitly acknowledged

The acknowledgments say that large portions of the GPU resources came from:

- **Genesis Therapeutics**
- **US Department of Energy / NERSC**

Specifically, the paper mentions use of:

- **NERSC** via the `GenAI@NERSC` award

## What is not disclosed

The paper does **not** clearly state:

- number of GPUs used for Boltz-1 training
- GPU model used for Boltz-1 training
- total GPU-hours
- total wall-clock duration in days or weeks
- exact distributed training setup

## Important caveat

The paper mentions **A100 80GB GPUs** in the context of evaluation/inference filtering for out-of-memory cases. That should not be read as an explicit statement of the Boltz-1 training hardware.

## Bottom line

What can be stated confidently:

- Boltz-1 structure training ran for `68k` steps at batch size `128`
- it used crops of `384/3456` first, then `512/4608`
- it trained on a mixture of PDB and OpenFold distillation, then PDB only
- it used substantial GPU resources from Genesis Therapeutics and NERSC
- the paper does **not** disclose exact GPU count, GPU type, GPU-hours, or wall-clock training time
