# 06 Environmental Driver Parsimonious GAMM

## Purpose

Fit the final parsimonious environmental-diversity GAMM with each retained environmental predictor represented in the model form independently assessed upstream.

## Input Files

- `input/environmental_variables_selected_by_corr0.8_VIF10.csv`: screened, standardized pool-level modelling data from workflow 03.
- `input/parsimonious_environmental_predictors.csv`: the independently reduced predictor set from workflow 05, including `selected_form`, `selected_k` and `model_form_source`.
- `input/environmental_predictor_model_form_assessment.csv`: the complete workflow 04 seven-criterion model-form assessment.

The three files are copied into this module so it can run as a self-contained workflow. Workflow 04 compares linear and smooth candidates with ML and records the seven selection criteria; workflow 05 performs the separate parsimony analysis and propagates those decisions for retained variables. Workflow 06 validates both sources, constructs the final formula dynamically and fits it with REML.

## Script

- `scripts/fit_parsimonious_gamm.R`: validates the three input contracts, builds an explicit environmental-term mapping, fits the final model with `mgcv::gam(method = "REML")`, performs drop-one refits and exports tables and figures.

For a retained variable with `selected_form = linear`, the formula contains its internal `env_i` term. For a validated `selected_form = nonlinear`, the formula contains `s(env_i, bs='tp', k=<selected_k>)`. The current 13 retained predictors are all linear, but this is an input result rather than a hard-coded restriction. A nonlinear decision is accepted only when all seven upstream criteria are true and `selected_k` is a finite integer of at least 4.

The input contract is literal rather than permissive: the complete workflow 04 assessment must contain exactly one row for each of the predefined 14 screened environmental predictors, with no missing, duplicate or unknown variables. Only after validating all 14 assessment rows does workflow 06 select and compare the 13 retained predictors from workflow 05. Model forms must be exactly `linear` or `nonlinear`; C1-C7 must be exactly `TRUE` or `FALSE`; a linear `selected_k` must be genuinely empty/NA; and a nonlinear `selected_k` must be a non-empty, strictly numeric finite integer of at least 4. Across the complete workflow 04 assessment, `selected_form = nonlinear` if and only if C1-C7 are all `TRUE`; therefore an all-`TRUE` row labelled linear and a nonlinear row containing any `FALSE` are both rejected. Workflow 04 and workflow 05 are validated independently before their retained-variable decisions are compared.

The standalone modelling table is also validated strictly. The response, pool size, coordinates and all retained environmental columns must convert losslessly to finite numeric values wherever the source value is non-missing. All 13 named land-use dummy columns are required, may not be missing or non-finite, and must contain only 0/1; missing columns are never synthesized and missing dummy values are never changed to zero. Non-missing factor values may not be empty. Complete cases are then selected from the actual model fields, with original, retained and removed row counts reported. Each categorical model field must retain at least two levels, coordinates must provide at least 10 unique pairs, and complete sample size must be at least both 50 and 11 greater than the parametric design rank.

## Expected Outputs

### Figures

1. `output/figures/environmental_observed_associations_key_variables.pdf`
   - Descriptive, unadjusted observed Shannon-diversity associations for significant retained environmental terms, plus any validated nonlinear retained term so its assessed shape is visible; a deterministic environmental fallback is used if neither set is available.
   - Linear terms use an `lm` line and confidence band. Nonlinear terms use `mgcv::gam(y ~ s(x, bs='tp', k=selected_k), method='REML')` with a smooth confidence band. These panels describe the observed bivariate association and are not plots of the adjusted final-model partial effect.
2. `output/figures/ML_drop_one_contribution.pdf`
   - Environmental drop-one changes in deviance explained from fixed-structure comparisons in which both the full and reduced models are fitted with ML.

### Tables

1. `output/tables/environmental_predictor_model_form_check.csv`
   - Pre-fit and final audit with exactly `variable`, `C1`-`C7`, `upstream_selected_form`, `upstream_selected_k`, `final_formula_term`, `final_model_form` and `action`.
   - `action` records that each upstream decision was validated and implemented.
2. `output/tables/model_fit_summary.csv`
   - Term-level estimates or smooth statistics, including `term`, `role`, `model_form`, `term_type`, `estimate`, `edf`, `statistic` and `p_value`.
   - Smooth internal names such as `s(env_i)` are mapped back to the original environmental label.
3. `output/tables/gamm_drop_one_contribution_summary.csv`
   - Drop-one results with `term_group`, `term_name`, `model_comparison_method = ML`, and ML-based changes in deviance explained, adjusted R-squared and AIC. The corresponding full-model p-value, where available, is taken from the final REML summary and is not an ML model-comparison p-value.

## Model Structure

Environmental terms alone inherit the linear/nonlinear classification from workflow 04. Sampling effort is a linear adjustment. Season, year and broad land-use are categorical adjustments. Terrain and host genus use random-effect smooths (`s(terrain, bs='re')` and `s(host_genus, bs='re')`), while broad residual spatial structure uses `s(longitude, latitude, bs='tp', k=10)`. These random-effect and spatial terms are adjustment structures, not environmental model-form classifications.

The reported final model is fitted with REML. Because drop-one candidates have different fixed-effect structures, the comparison full model and every reduced model are fitted separately with ML; their AIC, adjusted R-squared and deviance-explained differences are therefore ML comparisons. Term-level estimates and `full_model_term_p_value_if_available` remain descriptive results from the final REML summary. The explicit internal-to-label mapping is used consistently for summary extraction, term roles, model forms, p-values, drop-one removal and key-variable selection, so future validated nonlinear predictors do not become confused with random-effect or spatial smooths.
