# Deep-hole-extended ESGRS enumeration data

This folder contains reproducibility material for the manuscript
“Global Monomial Equivalence of Deep-Hole Extensions of ESGRS Codes to
Roth--Lempel Codes and Exact Enumeration in the Normalized `[7,4]_q` Family”.

The scripts use exact finite-field arithmetic.  The CSV file records the
small-prime-field values reported in the manuscript, including the MDS
parameter count, tangent-pair count, globally Roth--Lempel-equivalent count,
double-overlap count, and the quotient by the explicitly proved subgroup
`G_q ~= F_q^* x C_2`.

## Contents

- `exact_counts.csv`: exact values for `q = 7, 11, 13, 17, 19, 23`.
- `supplementary_computations.py`: compact regression and counting routines.
- `esgrs_projective_equiv.py`: full finite-field and projective-equivalence
  implementation used for the computations.
- `supplementary_readme.txt`: brief run notes from the submission package.

The code requires Python 3.10 or later and has no third-party dependencies.
All computations are deterministic; the reported values can be regenerated
with the assertions in the supplementary scripts.

## Citation

Yong Zhang, “Global Monomial Equivalence of Deep-Hole Extensions of ESGRS
Codes to Roth--Lempel Codes and Exact Enumeration in the Normalized `[7,4]_q`
Family,” manuscript, 2026.

Author: Yong Zhang, School of Mathematics and Statistics, Yancheng Teachers
University, Yancheng 224002, P. R. China. ORCID:
https://orcid.org/0000-0002-5917-3753.
