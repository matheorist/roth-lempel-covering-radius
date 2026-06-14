# Magma programs

This directory contains the Magma scripts used for the exact finite-field
computations in the paper:

**A complete covering-radius threshold for length-eight Roth--Lempel codes via trisecant-plane covers**

The programs verify affine classifications, covering radii, deep-hole
certificates, and boundary-field computations for the length-eight
Roth--Lempel codes `RL_{4,delta}(S)` with `|S| = 6`.

## Directory role

- `magma/`: Magma source files.
- `logs/`: output logs from exact Magma runs.
- `certificates/`: explicit boundary certificates.
- `data/`: extracted tables and CSV files.
- `paper/`: manuscript and supplementary material.

## Reproducibility

All computations use exact finite-field arithmetic in Magma. No randomized
search is used for theorem-level classifications or certificates.

The main proof-critical computations are:

1. Complete exceptional-field classification for `q = 7,8,9,11,13`.
2. Boundary certification for `q = 16,17,19`.
3. Verification of the explicit deep-hole product criterion.
4. Verification of the syndrome distance stratification and shell formulas.

Additional larger-field scans are included as independent checks.
