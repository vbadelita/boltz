# Flu HA Trimer Single-Ligand Affinity Screen

## Overview

- Conditions: 16
- Completed conditions: 16
- Completed affinity predictions: 16
- Successful configs: `ha_trimer_1lig_aff_s10_n1_p1_r0_as200_an5`, `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5`, `ha_trimer_1lig_aff_s80_n3_p1_r5_as200_an5`
- HA modeled as a homotrimer: chains `A`, `B`, `C`.
- Ligand modeled as a single affinity binder: chain `D`.

## Summary

| Protein | Ligand | Config | Status | Affinity | Binder prob | MSA depth | pLDDT | ipTM | ligand ipTM |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| A/duck/Bavaria/1/77 | alpha2-3 SIA glycan | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 1.363 | 0.472 | 1227 | 0.904 | 0.877 | 0.845 |
| A/duck/Bavaria/1/77 | alpha2-6 SIA glycan | `ha_trimer_1lig_aff_s80_n3_p1_r5_as200_an5` | ok | 1.632 | 0.423 | 1227 | 0.908 | 0.887 | 0.808 |
| A/duck/Bavaria/1/77 | NGC | `ha_trimer_1lig_aff_s10_n1_p1_r0_as200_an5` | ok | 1.271 | 0.414 | 1227 | 0.462 | 0.378 | 0.429 |
| A/duck/Bavaria/1/77 | SIA | `ha_trimer_1lig_aff_s10_n1_p1_r0_as200_an5` | ok | -0.040 | 0.372 | 1227 | 0.454 | 0.413 | 0.440 |
| A/Memphis/7/1980 | alpha2-3 SIA glycan | `ha_trimer_1lig_aff_s80_n3_p1_r5_as200_an5` | ok | 1.154 | 0.499 | 1224 | 0.918 | 0.897 | 0.810 |
| A/Memphis/7/1980 | alpha2-6 SIA glycan | `ha_trimer_1lig_aff_s80_n3_p1_r5_as200_an5` | ok | 1.740 | 0.436 | 1224 | 0.917 | 0.890 | 0.801 |
| A/Memphis/7/1980 | NGC | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 1.578 | 0.536 | 1224 | 0.921 | 0.886 | 0.866 |
| A/Memphis/7/1980 | SIA | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 0.901 | 0.652 | 1224 | 0.919 | 0.881 | 0.887 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | alpha2-3 SIA glycan | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 1.264 | 0.557 | 1235 | 0.917 | 0.858 | 0.869 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | alpha2-6 SIA glycan | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 1.678 | 0.514 | 1235 | 0.915 | 0.858 | 0.793 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | NGC | `ha_trimer_1lig_aff_s80_n3_p1_r5_as200_an5` | ok | 1.124 | 0.579 | 1235 | 0.925 | 0.884 | 0.900 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | SIA | `ha_trimer_1lig_aff_s80_n3_p1_r5_as200_an5` | ok | 0.342 | 0.752 | 1235 | 0.920 | 0.873 | 0.947 |
| A/dairy_cow/Texas/24_009108-001/2024 | alpha2-3 SIA glycan | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 1.298 | 0.526 | 1237 | 0.919 | 0.871 | 0.853 |
| A/dairy_cow/Texas/24_009108-001/2024 | alpha2-6 SIA glycan | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 1.698 | 0.560 | 1237 | 0.919 | 0.879 | 0.845 |
| A/dairy_cow/Texas/24_009108-001/2024 | NGC | `ha_trimer_1lig_aff_s80_n3_p1_r5_as200_an5` | ok | 0.994 | 0.580 | 1237 | 0.921 | 0.877 | 0.905 |
| A/dairy_cow/Texas/24_009108-001/2024 | SIA | `ha_trimer_1lig_aff_s40_n1_p1_r3_as200_an5` | ok | 0.304 | 0.763 | 1237 | 0.922 | 0.883 | 0.941 |

## Files

- Metadata CSV: `flu_trimer_single_ligand_affinity/metadata.csv`
- Structures: `flu_trimer_single_ligand_affinity/structures/`
- Raw outputs: `flu_trimer_single_ligand_affinity/runs/`
- Logs: `flu_trimer_single_ligand_affinity/logs/`
