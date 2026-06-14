# Codex supplement: complete low-field deep-hole census

Date: 2026-05-25  
Family: `C = RL_{4,delta}(S)`, `|S|=6`, `q in {7,8,9,11,13}`

## Conversion principle

Let `h` be the number of projective syndrome points not covered by the
trisecant planes. For every class with `rho=4`,

```text
number of deep-hole cosets = (q-1) h,
number of deep-hole words  = q^4 (q-1) h.
```

For a class with `rho=3`, the log reports the number `m_3` of projective
syndromes in the third distance shell; then

```text
number of deep-hole cosets = (q-1) m_3,
number of deep-hole words  = q^4 (q-1) m_3.
```

## Complete radius-four defect spectrum

Rows are grouped by finite field, code type, and projective-hole count. The
rightmost two columns give the deep-hole count for one code in the row.

| q | type | holes | affine classes | raw cases | DH cosets per code | DH words per code |
|---:|---|---:|---:|---:|---:|---:|
| 8 | NMDS | 1 | 1 | 56 | 7 | 28,672 |
| 11 | MDS | 1 | 1 | 22 | 10 | 146,410 |
| 11 | NMDS | 1 | 13 | 1,430 | 10 | 146,410 |
| 11 | NMDS | 2 | 2 | 220 | 20 | 292,820 |
| 13 | MDS | 1 | 4 | 624 | 12 | 342,732 |
| 13 | MDS | 2 | 1 | 156 | 24 | 685,464 |
| 13 | MDS | 3 | 2 | 312 | 36 | 1,028,196 |
| 13 | NMDS | 1 | 11 | 1,716 | 12 | 342,732 |
| 13 | NMDS | 2 | 34 | 5,304 | 24 | 685,464 |
| 13 | NMDS | 3 | 22 | 3,328 | 36 | 1,028,196 |
| 13 | NMDS | 4 | 27 | 4,082 | 48 | 1,370,928 |
| 13 | NMDS | 5 | 15 | 2,340 | 60 | 1,713,660 |
| 13 | NMDS | 6 | 10 | 1,560 | 72 | 2,056,392 |
| 13 | NMDS | 7 | 8 | 1,248 | 84 | 2,399,124 |
| 13 | NMDS | 8 | 2 | 234 | 96 | 2,741,856 |
| 13 | NMDS | 9 | 2 | 312 | 108 | 3,084,588 |
| 13 | NMDS | 12 | 1 | 156 | 144 | 4,112,784 |

There are no radius-four classes at `q=7` or `q=9`.

## Complete radius-three deep-hole spectrum

| q | type | DH cosets per code | DH words per code | affine classes | raw cases |
|---:|---|---:|---:|---:|---:|
| 7 | NMDS | 1,380 | 3,313,380 | 1 | 7 |
| 7 | NMDS | 1,398 | 3,356,598 | 1 | 42 |
| 8 | NMDS | 2,709 | 11,096,064 | 3 | 168 |
| 9 | MDS | 4,704 | 30,862,944 | 1 | 24 |
| 9 | NMDS | 4,752 | 31,177,872 | 8 | 516 |
| 9 | NMDS | 4,776 | 31,335,336 | 5 | 216 |
| 11 | MDS | 11,760 | 172,178,160 | 1 | 55 |
| 11 | NMDS | 11,790 | 172,617,390 | 14 | 1,540 |
| 11 | NMDS | 11,820 | 173,056,620 | 16 | 1,705 |
| 11 | NMDS | 11,850 | 173,495,850 | 1 | 110 |
| 13 | MDS | 24,432 | 697,802,352 | 3 | 312 |
| 13 | NMDS | 24,468 | 698,830,548 | 4 | 624 |

## Reproducibility

- Script: `magma/codex_phase7_d3_complete_small_primepowers_q7_q9.m`
- Wrapper: `magma/codex_phase7_d3_defect_spectrum_q7_q13_codex.m`
- Complete shell log:
  `magma/codex_phase7_d3_deephole_shell_spectrum_q7_q13_codex.log`

The log contains each individual affine representative, its orbit size,
MDS/NMDS type, covering radius, projective distance-shell vector, deep-hole
coset count, and deep-hole word count.
