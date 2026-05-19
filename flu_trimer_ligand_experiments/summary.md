# Flu HA Trimer Ligand Screen

## Overview

- Conditions: 16
- Completed conditions: 16
- Config: `ha_trimer_3lig_s40_n1_p1_r3`
- HA modeled as a homotrimer: chains `A`, `B`, `C`.
- Ligands modeled as three identical copies: chains `D`, `E`, `F`.
- Affinity JSON is intentionally absent: Boltz currently rejects affinity calculation for multi-copy ligands.

## Summary

| Protein | Ligand | Status | MSA depth | Raw MSA depth | complex_plddt | ptm | iptm | ligand_iptm | confidence_score |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| A/duck/Bavaria/1/77 | alpha2-3 SIA glycan | ok | 1227 | 1234 | 0.909 | 0.914 | 0.904 | 0.876 | 0.908 |
| A/duck/Bavaria/1/77 | alpha2-6 SIA glycan | ok | 1227 | 1234 | 0.904 | 0.889 | 0.877 | 0.839 | 0.899 |
| A/duck/Bavaria/1/77 | NGC | ok | 1227 | 1234 | 0.920 | 0.903 | 0.894 | 0.874 | 0.915 |
| A/duck/Bavaria/1/77 | SIA | ok | 1227 | 1234 | 0.916 | 0.891 | 0.883 | 0.898 | 0.909 |
| A/Memphis/7/1980 | alpha2-3 SIA glycan | ok | 1224 | 1231 | 0.915 | 0.885 | 0.873 | 0.829 | 0.906 |
| A/Memphis/7/1980 | alpha2-6 SIA glycan | ok | 1224 | 1231 | 0.909 | 0.889 | 0.879 | 0.819 | 0.903 |
| A/Memphis/7/1980 | NGC | ok | 1224 | 1231 | 0.923 | 0.897 | 0.888 | 0.835 | 0.916 |
| A/Memphis/7/1980 | SIA | ok | 1224 | 1231 | 0.925 | 0.903 | 0.896 | 0.863 | 0.919 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | alpha2-3 SIA glycan | ok | 1235 | 1242 | 0.905 | 0.869 | 0.854 | 0.814 | 0.895 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | alpha2-6 SIA glycan | ok | 1235 | 1242 | 0.910 | 0.877 | 0.863 | 0.852 | 0.901 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | NGC | ok | 1235 | 1242 | 0.921 | 0.870 | 0.859 | 0.896 | 0.909 |
| A/cattle/CA/24-036570-001-original/2024(H5N1) | SIA | ok | 1235 | 1242 | 0.919 | 0.880 | 0.870 | 0.925 | 0.909 |
| A/dairy_cow/Texas/24_009108-001/2024 | alpha2-3 SIA glycan | ok | 1237 | 1244 | 0.910 | 0.882 | 0.869 | 0.824 | 0.901 |
| A/dairy_cow/Texas/24_009108-001/2024 | alpha2-6 SIA glycan | ok | 1237 | 1244 | 0.911 | 0.879 | 0.866 | 0.848 | 0.902 |
| A/dairy_cow/Texas/24_009108-001/2024 | NGC | ok | 1237 | 1244 | 0.922 | 0.899 | 0.890 | 0.906 | 0.915 |
| A/dairy_cow/Texas/24_009108-001/2024 | SIA | ok | 1237 | 1244 | 0.921 | 0.901 | 0.893 | 0.926 | 0.915 |

## Files

- Metadata CSV: `flu_trimer_ligand_experiments/metadata.csv`
- Structures: `flu_trimer_ligand_experiments/structures/`
- Raw outputs: `flu_trimer_ligand_experiments/runs/`
- Logs: `flu_trimer_ligand_experiments/logs/`
