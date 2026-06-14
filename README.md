# Roth--Lempel covering radius

This repository contains computational material for the paper:

**A complete covering-radius threshold for length-eight Roth--Lempel codes via trisecant-plane covers**

## Main result

For every prime power `q >= 16`, every length-eight Roth--Lempel code
`RL_{4,delta}(S)` with `|S| = 6` has covering radius `4`.

For `q = 7,8,9,11,13`, the repository contains the complete affine
classification of the exceptional radius-three cases.

## Repository structure

- `paper/`: manuscript, PDF, supplements, and reproducibility manifest.
- `magma/`: Magma scripts.
- `logs/`: exact output logs from Magma runs.
- `certificates/`: explicit boundary certificates for `q = 16,17,19`.
- `data/`: CSV tables and extracted affine-class data.

## Reproducibility

All theorem-level computations use exact finite-field arithmetic in Magma.
No randomized search is used for classifications or certificates.

The proof-critical computations include:

1. Complete exceptional-field classification for `q = 7,8,9,11,13`.
2. Boundary certification for `q = 16,17,19`.
3. Verification of the explicit deep-hole product criterion.
4. Verification of the syndrome distance stratification and shell formulas.

## Author

Yong Zhang  
School of Mathematics and Statistics, Yancheng Normal University  
Yancheng, Jiangsu 224002, China
