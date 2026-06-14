# Study of `Roth_Lempel_trisecant_covers_DM_codex修改意见.tex`

## Main usable recommendations

1. Replace hard-coded proposition numbers by dynamic LaTeX references.
2. Add a direct reference to `cor:q13-seven` when discussing the seven
   exceptional classes over `F_13`.
3. Strengthen the finite-geometric interpretation of uncovered syndrome
   points as projection centers producing plane `8`-arcs.
4. Add a cautious matroid/Tutte-polynomial bridge for the second syndrome
   shell.
5. Clarify the code-availability path for review and final publication.
6. Place the `Ball2015` citation where affine normalization and projective
   configuration equivalence are actually used.

## Verification against the current manuscript

- The alleged hard-coded `Proposition 1` and `Proposition 2` problems were
  already fixed in the current manuscript: the source uses
  `Proposition~\ref{prop:exhaustion}` and
  `Proposition~\ref{prop:duality}`.
- `cor:q13-seven` was indeed defined but not referenced; this has now been
  fixed in the Discussion.
- The broader geometry and matroid comments were useful, but the wording in
  the suggestion file was too promotional for a journal article. The revised
  manuscript uses lower-key, mathematically checkable statements.
- A fake public repository URL was not inserted. Instead, the reproducibility
  section now states that the anonymous review package is supplied
  electronically and that the final version will use an archival repository
  or DOI link.

## Changes made

- Added a short finite-geometric paragraph after Proposition
  `prop:projection`.
- Added a `Ball2015` citation in the complete-affine-classification setup.
- Added a dynamic reference to Corollary `cor:q13-seven`.
- Added a matroid/Tutte-polynomial sentence in the Discussion and cited
  Greene's foundational paper:
  Curtis Greene, "Weight enumeration and the geometry of linear codes,"
  *Studies in Applied Mathematics* 55(2), 119--128, 1976.
- Added a code-availability sentence to Section 6.

## Checks

- Distinct cited keys: 13.
- Bibliography entries: 13.
- Missing bibliography entries: none.
- Uncited bibliography entries: none.
- Two-pass LaTeX compilation completed with no unresolved references and no
  box warnings.
