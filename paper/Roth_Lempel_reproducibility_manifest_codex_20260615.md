# Reproducibility manifest (codex version, updated 2026-06-15)

This manifest accompanies `Roth_Lempel_trisecant_covers_DM_codex.tex`.
It separates computations required by the mathematical argument from larger
independent scans retained for checking and future extension. It supersedes
`Roth_Lempel_reproducibility_manifest_codex_20260526.md` by making explicit
that the `q=16` two-slice scan, the extension-field scans, and the thirty
boundary certificates are included in the reproducibility record.

## Proof-critical computations

| Role in the paper | Field range | Program(s) | Output log(s) |
| --- | --- | --- | --- |
| Complete exceptional classification and projective shells | `q=7,8,9,11,13` | `magma/codex_phase7_d3_complete_small_primepowers_q7_q9.m`; `magma/codex_phase7_d3_defect_spectrum_q7_q13_codex.m` | `magma/codex_phase7_d3_complete_small_primepowers_q7_q9_codex.log`; `magma/codex_phase7_d3_deephole_shell_spectrum_q7_q13_codex.log` |
| Explicit product criterion | `q=7,...,16` | `magma/codex_phase8_explicit_deephole_product_criterion.m` | `magma/codex_phase8_explicit_deephole_product_criterion_codex.log` |
| Distance-stratification and shell formula check | `q=7,...,13` | `magma/codex_phase9_explicit_distance_stratification_codex.m`; `magma/codex_phase9_lowfield_triplesum_hole_spectrum_codex.m` | `magma/codex_phase9_explicit_distance_lowfield_complete_codex.log`; `magma/codex_phase9_lowfield_triplesum_hole_spectrum_codex.log` |
| Residual boundary-count determination | `q=16,17,19` | `magma/codex_round42_stratified_boundary_need.m` | `magma/codex_round42_stratified_boundary_need_codex.log` |
| Explicit residual boundary certificates | `q=16,17,19` | `magma/codex_rounds11_18_boundary_certificate_analysis.m`; `magma/codex_round24_infinite_plane_centers_boundary.m` | `magma/codex_rounds11_18_boundary_certificate_analysis_codex.log`; `magma/codex_round24_infinite_plane_centers_boundary_codex.log` |

The thirty displayed residual certificates are supplied in
`Roth_Lempel_boundary_certificates_codex_supplement.tex`.

## Formula and covariance verification

| Verified statement | Program | Output log |
| --- | --- | --- |
| Closed second-shell formula on boundary instances | `magma/codex_round16_m2_formula_boundary_codex.m` | `magma/codex_round16_m2_formula_boundary_codex.log` |
| Low-field explicit shell tabulation by triple-sum multiplicity | `magma/codex_round16_lowfield_explicit_shell_by_r_codex.m` | `magma/codex_round16_lowfield_explicit_shell_by_r_codex.log` |
| Affine covariance of the classification | `magma/codex_round17_affine_covariance_check_codex.m` | `magma/codex_round17_affine_covariance_check_codex.log` |
| Refined theoretical lower bound by `r_S(delta)` | `magma/codex_round23_r_refined_bound_check.m` | `magma/codex_round23_r_refined_bound_check_codex.log` |
| Dual-forbidden-set pairing calculation | `magma/codex_round30_dual_forbidden_pairing.m` | `magma/codex_round30_dual_forbidden_pairing_codex.log` |

## Supplementary independent scans

These computations corroborate the threshold theorem but are not used in its
proof after the theoretical `q >= 23` argument.

| Field range | Program family | Representative logs |
| --- | --- | --- |
| `q=17,19,23,29,31,37,41` | `magma/codex_phase7_d3_restricted_centers_q17_q41.m` and structured-then-full wrappers | `magma/codex_phase7_d3_structured_then_full_q17_codex.log`; `magma/codex_phase7_d3_structured_then_full_q19_q23_codex.log`; `magma/codex_phase7_d3_structured_then_full_q29_q37_codex.log`; `magma/codex_phase7_d3_structured_then_full_q41_codex.log` |
| `q=16,25,27,32,49` | `magma/codex_phase7_d3_multislice_extension_scan.m` | `magma/codex_phase7_d3_extension_q16_only_codex.log`; `magma/codex_phase7_d3_extension_q25_only_codex.log`; `magma/codex_phase7_d3_extension_q27_only_codex.log`; `magma/codex_phase7_d3_extension_q32_only_codex.log`; `magma/codex_phase7_d3_extension_q49_only_codex.log` |
| `q=43,47,53` | `magma/codex_phase7_d3_multislice_prime_scan_q43_q47.m`; `magma/codex_phase7_d3_multislice_q53_full_codex.m` | `magma/codex_phase7_d3_multislice_q43_only_codex.log`; `magma/codex_phase7_d3_multislice_q47_only_codex.log`; `magma/codex_phase7_d3_multislice_q53_full_codex.log` |

The `q=16` log above is the two-slice scan cited in the proof of the complete
small-field classification. The same boundary field is also covered by the
theoretical threshold proof after the twelve residual representatives are
certified by the boundary-certificate supplement.

## Data extracts

- Complete low-field census: `plans/codex_low_field_deephole_census_supplement_20260525.md`.
- Full representative data for `q=11,13,17`:
  `Roth_Lempel_trisecant_covers_DM_codex_full_classes.tex`.
- Machine-readable scan tables for `q=11,13,17,19,23,29,31,37` are in
  `plans/` with prefix `codex_phase7_d3_`.

All programs use exact Magma finite-field arithmetic. The classification
scripts enumerate affine-orbit representatives; no randomized search is used
for any theorem-level certificate.
