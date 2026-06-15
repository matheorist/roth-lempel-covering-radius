# Reproducibility manifest (updated 2026-06-15)

This manifest accompanies `Roth_Lempel_covering_radius.tex`.
It separates computations required by the mathematical argument from larger
independent scans retained for checking and future extension. It makes explicit
that the `q=16` two-slice scan, the extension-field scans, and the thirty
boundary certificates are included in the reproducibility record.

## Proof-critical computations

| Role in the paper | Field range | Program(s) | Output log(s) |
| --- | --- | --- | --- |
| Complete exceptional classification and projective shells | `q=7,8,9,11,13` | `magma/phase7_d3_complete_small_primepowers_q7_q9.m`; `magma/phase7_d3_defect_spectrum_q7_q13.m` | `magma/phase7_d3_complete_small_primepowers_q7_q9.log`; `magma/phase7_d3_deephole_shell_spectrum_q7_q13.log` |
| Explicit product criterion | `q=7,...,16` | `magma/phase8_explicit_deephole_product_criterion.m` | `magma/phase8_explicit_deephole_product_criterion.log` |
| Distance-stratification and shell formula check | `q=7,...,13` | `magma/phase9_explicit_distance_stratification.m`; `magma/phase9_lowfield_triplesum_hole_spectrum.m` | `magma/phase9_explicit_distance_lowfield_complete.log`; `magma/phase9_lowfield_triplesum_hole_spectrum.log` |
| Residual boundary-count determination | `q=16,17,19` | `magma/round42_stratified_boundary_need.m` | `magma/round42_stratified_boundary_need.log` |
| Explicit residual boundary certificates | `q=16,17,19` | `magma/rounds11_18_boundary_certificate_analysis.m`; `magma/round24_infinite_plane_centers_boundary.m` | `magma/rounds11_18_boundary_certificate_analysis.log`; `magma/round24_infinite_plane_centers_boundary.log` |

The thirty displayed residual certificates are supplied in `supplement.tex`.

## Formula and covariance verification

| Verified statement | Program | Output log |
| --- | --- | --- |
| Closed second-shell formula on boundary instances | `magma/round16_m2_formula_boundary.m` | `magma/round16_m2_formula_boundary.log` |
| Low-field explicit shell tabulation by triple-sum multiplicity | `magma/round16_lowfield_explicit_shell_by_r.m` | `magma/round16_lowfield_explicit_shell_by_r.log` |
| Affine covariance of the classification | `magma/round17_affine_covariance_check.m` | `magma/round17_affine_covariance_check.log` |
| Refined theoretical lower bound by `r_S(delta)` | `magma/round23_r_refined_bound_check.m` | `magma/round23_r_refined_bound_check.log` |
| Dual-forbidden-set pairing calculation | `magma/round30_dual_forbidden_pairing.m` | `magma/round30_dual_forbidden_pairing.log` |

## Supplementary independent scans

These computations corroborate the threshold theorem but are not used in its
proof after the theoretical `q >= 23` argument.

| Field range | Program family | Representative logs |
| --- | --- | --- |
| `q=17,19,23,29,31,37,41` | `magma/phase7_d3_restricted_centers_q17_q41.m` and structured-then-full wrappers | `magma/phase7_d3_structured_then_full_q17.log`; `magma/phase7_d3_structured_then_full_q19_q23.log`; `magma/phase7_d3_structured_then_full_q29_q37.log`; `magma/phase7_d3_structured_then_full_q41.log` |
| `q=16,25,27,32,49` | `magma/phase7_d3_multislice_extension_scan.m` | `magma/phase7_d3_extension_q16_only.log`; `magma/phase7_d3_extension_q25_only.log`; `magma/phase7_d3_extension_q27_only.log`; `magma/phase7_d3_extension_q32_only.log`; `magma/phase7_d3_extension_q49_only.log` |
| `q=43,47,53` | `magma/phase7_d3_multislice_prime_scan_q43_q47.m`; `magma/phase7_d3_multislice_q53_full.m` | `magma/phase7_d3_multislice_q43_only.log`; `magma/phase7_d3_multislice_q47_only.log`; `magma/phase7_d3_multislice_q53_full.log` |

The `q=16` log above is the two-slice scan cited in the proof of the complete
small-field classification. The same boundary field is also covered by the
theoretical threshold proof after the twelve residual representatives are
certified by the boundary-certificate supplement.

### Historical one-slice scan summary

Before the two-slice scan was run, a one-slice scan searched first in
`P=[1,u,v,v]` and invoked a full-space fallback only when that slice contained
no hole. It is retained as a cross-check, not as an independent proof input in
the final article.

| `q` | affine classes | `P=[1,u,v,v]` certificates | fallback certificates | uncertified |
| --- | ---: | ---: | ---: | ---: |
| 17 | 777 | 695 | 82 | 0 |
| 19 | 1514 | 1415 | 99 | 0 |
| 23 | 4596 | 4396 | 200 | 0 |
| 29 | 16978 | 16164 | 814 | 0 |
| 31 | 24562 | 23434 | 1128 | 0 |
| 37 | 64604 | 61632 | 2972 | 0 |
| 41 | 112439 | 107526 | 4913 | 0 |

## Data extracts

- Complete low-field census: `plans/low_field_deephole_census_supplement_20260525.md`.
- Full representative data for `q=11,13,17`:
  `Roth_Lempel_complete_affine_classes.tex`.
- Machine-readable scan tables for `q=11,13,17,19,23,29,31,37` are in
  `plans/` with prefix `phase7_d3_`.

All programs use exact Magma finite-field arithmetic. The classification
scripts enumerate affine-orbit representatives; no randomized search is used
for any theorem-level certificate.
