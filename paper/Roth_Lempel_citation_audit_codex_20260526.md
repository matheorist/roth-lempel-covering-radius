# Citation audit report (codex version)

Manuscript checked: `Roth_Lempel_trisecant_covers_DM_codex.tex`

## Audit outcome

The previous draft had eight bibliography entries but used only four of them
in the text. It cited the original Roth--Lempel construction and one
Roth--Lempel NMDS paper, but did not position the manuscript against the two
closest 2025--2026 developments in covering radii/deep holes and published
Roth--Lempel analysis.

After revision, the manuscript contains 12 distinct cited keys and 12
bibliography entries. There are no missing cited entries and no uncited
bibliography entries. Two-pass LaTeX compilation of the citation-revised PDF
reported no unresolved references or box warnings.

## Corrections made

| Location | Revision | Reason |
| --- | --- | --- |
| Introduction | Cited standard covering-code and coding-theory sources (`Cohen1997`, `Huffman2003`, `MWS1977`) | Supports the statement that covering radius and deep holes are classical coding-theoretic questions. |
| Introduction | Cited finite-geometric and NMDS foundational sources (`Ball2015`, `DodunekovLandgev1995`) | Supports projective/arcs formulation and the NMDS terminology. |
| Introduction | Added related-work positioning using `XuZhou2026`, `LiangLiao2026`, and `LiZhuSun2025` | Distinguishes the manuscript's covering-radius threshold result from recent NMDS/weight-distribution and deep-hole construction literature. |
| Roth--Lempel codes section | Added `DodunekovLandgev1995` at the NMDS convention | Attributes the standard NMDS definition. |
| Trisecant-plane proof | Added `Cohen1997` and `Huffman2003` at the redundancy bound | Attributes the standard covering-radius bound used in the proof. |
| Bibliography | Added four directly relevant entries and DOI metadata for journal papers already cited | Improves traceability and removes the central literature gap. |

## Newly added directly relevant references

1. S. M. Dodunekov and I. N. Landgev, "On near-MDS codes," *Journal of Geometry*, 54(1--2), 30--43, 1995. DOI: `10.1007/BF01222850`.
2. Y. Li, S. Zhu, and Z. Sun, "Covering radii and deep holes of two classes of extended twisted GRS codes and their applications," *IEEE Transactions on Information Theory*, 71(5), 3516--3530, 2025. DOI: `10.1109/TIT.2025.3541799`.
3. H. Xu and H. Zhou, "Analysis of Roth--Lempel codes," *IEEE Transactions on Information Theory*, 72(1), 246--252, 2026. DOI: `10.1109/TIT.2025.3626824`.
4. Z. Liang and Q. Liao, "Two classes of NMDS codes from Roth--Lempel codes," *Finite Fields and Their Applications*, 111, 102779, 2026. DOI: `10.1016/j.ffa.2025.102779`.

## Scope decision

Works on quantum-code applications, decoding algorithms, and cryptanalysis
were not added to the short paper because the manuscript does not make claims
in those directions. They would enlarge the bibliography without supporting
the present covering-radius and deep-hole classification results.
