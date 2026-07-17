# Environmental linear/nonlinear form assessment

## Purpose

This workflow classifies the functional form of each of the 14 correlation/VIF-screened environmental predictors before final GAMM fitting. It compares a linear focal term with a thin-plate regression spline while holding the analysis rows and all other adjustment terms constant. The result is a conservative, auditable choice of `linear` or `nonlinear`; it is not a test of whether an environmental predictor should survive the separate downstream parsimony step.

The ML fits here are used only to compare candidate functional forms by AIC and supporting diagnostics. Final coefficient estimation and inference use REML in workflow 06.

## Input files

- `input/environmental_variables_selected_by_corr0.8_VIF10.csv` is the modelling table produced by workflow 03 at `../03_environmental-correlation-vif-screening/output/tables/environmental_variables_selected_by_corr0.8_VIF10.csv`.
- The copy under this module's `input/` directory is deliberately included so this workflow can be run locally without writing into or depending at run time on the upstream workflow directory.
- `r_requirements.txt` records the required R and `mgcv` versions.

The script requires the Shannon response, log10 pool size, season, year, terrain, host genus, longitude, latitude, the 13 original land-use dummy columns, and all 14 screened environmental variables. Numeric conversion, finite values, required columns, and binary land-use encoding are checked explicitly.

## Script

Run from any working directory with:

```powershell
Rscript 05_environmental_association_and_gamm/04_environmental-linear-nonlinear-form-assessment/scripts/assess_environmental_model_forms.R
```

`scripts/assess_environmental_model_forms.R` resolves paths relative to its own location, reconstructs the broad land-use factor, assesses every environmental variable, and writes the tables and diagnostic PDF below.

All comparisons use one common complete-case dataset across the response, every screened environmental predictor, and every adjustment variable. Both members of each linear-versus-smooth comparison have the same adjustment structure: log10 pool size, season, year, broad land use, random-effect smooths for terrain and host genus, one two-dimensional longitude/latitude smooth, and the other 13 environmental predictors as linear terms. Only the focal predictor changes between a linear term and a candidate smooth. The random-effect smooths and spatial smooth are adjustment structures; they are not classifications of environmental-variable nonlinearity.

Models used for form screening are fitted with maximum likelihood (`method = "ML"`), which permits AIC comparison of candidate fixed-effect forms. The selected form is passed downstream; workflow 06 performs the final model fit with REML.

## Prespecified seven-criterion rule

The primary nonlinear candidate is a thin-plate spline with `k = 6`. Basis sensitivity normally compares `k = 4/6/8`; a predictor with fewer than eight unique values uses `k = 4/5/6` instead (currently FVC). A predictor is classified as `nonlinear` only when all seven criteria pass:

1. **C1 — AIC improvement:** `smooth_AIC - linear_AIC <= -2`.
2. **C2 — material curvature:** the primary smooth has `edf > 1.5`.
3. **C3 — smooth support:** the primary smooth has `p < 0.05`.
4. **C4 — adequate basis dimension:** the focal `k.check` result has `p >= 0.05`.
5. **C5 — distance from the EDF ceiling:** `(k - 1) - edf >= 0.5` for the primary `k = 6` smooth.
6. **C6 — stable curve across basis sizes:** the minimum pairwise curve correlation is `> 0.98` and the maximum curve difference relative to the primary curve's range is `< 0.20`.
7. **C7 — robustness to boundary observations:** after trimming both tails at 2.5%, the smooth still has delta AIC `<= -2`, `edf > 1.5`, and `p < 0.05`, and its partial-effect curve correlates with the full-data primary curve at `> 0.98`.

Any failed criterion retains the linear form. This deliberately makes isolated evidence of curvature insufficient for a nonlinear designation.

Core computations fail loudly: an unavailable primary linear fit, primary `k = 6` smooth, primary AIC, EDF, or p-value stops the workflow rather than producing an incomplete decision table. Failures in secondary diagnostics or sensitivity fits are recorded in `decision_reason`; the affected criterion evaluates as failed and the conservative result is `linear`.

## Outputs

1. `output/tables/environmental_predictor_model_form_assessment.csv` — the 14-row primary assessment table. Its 31 columns are:
   - identity and sample structure: `variable`, `n`, `unique_n`, `k_values_tested`, `primary_k`;
   - primary model comparison: `linear_AIC`, `smooth_AIC`, `delta_AIC`, `smooth_edf`, `smooth_p_value`;
   - basis diagnostics: `k_index`, `k_check_p_value`, `maximum_available_edf`, `edf_ceiling_margin`, `minimum_curve_correlation`, `maximum_relative_curve_difference`;
   - boundary-trimming diagnostics: `trimmed_n`, `trimmed_delta_AIC`, `trimmed_edf`, `trimmed_p_value`, `trimmed_curve_correlation`;
   - criterion flags: `criterion_1_aic`, `criterion_2_edf`, `criterion_3_p`, `criterion_4_k_check`, `criterion_5_edf_ceiling`, `criterion_6_curve_stability`, `criterion_7_boundary_robustness`;
   - decision fields: `selected_form`, `selected_k`, `decision_reason`.
2. `output/tables/environmental_predictor_model_form_criteria_summary.csv` — a fixed-format nine-row summary: seven criterion pass counts followed by the `linear` and `nonlinear` selection counts, each with the common total.
3. `output/figures/environmental_predictor_model_form_diagnostics.pdf` — a 14-page PDF, one page per predictor, showing the linear and candidate smooth partial-effect curves, data rug, 2.5% tail boundaries, criterion status, and final form decision.

## Current-data result

The current run uses 225 common complete cases. All 14 predictors are classified as `linear`; the fixed summary records 14 linear and 0 nonlinear selections.

Precipitation of Warmest Quarter (Bio18) provides the only partial indication of curvature: it passes C1 (AIC improvement), C2 (EDF), and C5 (EDF-ceiling margin), but fails C3 (smooth p-value), C4 (`k.check`), C6 (curve stability across basis sizes), and C7 (2.5%-tail-trim robustness). It therefore cannot be interpreted as a stable nonlinear effect, and the prespecified all-seven rule retains its linear form.

## Downstream interface

- Workflow 05 reads the computed model-form assessment and then applies its own independent parsimony rule. Functional-form selection and predictor removal are separate decisions.
- Workflow 06 reads `selected_form` (and `selected_k` when applicable) dynamically when constructing the final model, then fits that final specification with REML. It must not hard-code the incidental current-data result that all predictors are linear.

This interface allows a future rerun to propagate a justified nonlinear decision without silently changing the screening rule or conflating it with parsimony.
