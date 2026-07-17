# 05 Environmental Association and Generalized Additive Mixed Model (GAMM)

This block corresponds to the environmental association and adjusted environmental modelling analyses. The folder order follows manuscript presentation: correlation/Mantel summaries first, followed by the technical generalized additive mixed model (GAMM) preprocessing and modelling chain.

## Workflow Folders

| Order | Workflow | Purpose |
|---|---|---|
| 1 | `01_environmental-correlation-analysis/` | Produces environmental correlation, pairwise association and Mantel-style summaries. |
| 2 | `02_environmental-variable-standardization/` | Standardizes candidate environmental predictors for downstream GAMM screening and modelling. |
| 3 | `03_environmental-correlation-vif-screening/` | Validates the retained predictor set against correlation and VIF criteria and audits maximum-size feasible predictor sets. |
| 4 | `04_environmental-linear-nonlinear-form-assessment/` | Compares linear and smooth forms for each of the 14 screened predictors under a prespecified seven-criterion rule. |
| 5 | `05_environmental-final-variable-selection/` | Applies an independent contribution-based parsimony rule and propagates the assessed forms for retained predictors. |
| 6 | `06_environmental-driver-gamm-final-model/` | Constructs environmental terms dynamically from the retained form metadata, fits the final adjusted GAMM with REML and exports model summaries and PDF figures. |

Technical workflow order: run workflows 2-6 in order for the GAMM chain. Workflow 1 is placed first because it appears earlier in the manuscript narrative but is not required as a direct input to the GAMM chain.

## Method Rationale

The environmental block separates exploratory association summaries from the parsimonious GAMM chain. Correlation and Mantel summaries document broad environmental structure, while standardization, correlation/VIF validation, model-form assessment, parsimony-based variable selection and GAMM fitting form the reproducible modelling chain. This separation makes it clear which outputs are descriptive checks and which outputs are direct modelling inputs.

The retained 14-predictor environmental set is documented as a maximum-size set satisfying the combined absolute Spearman correlation < 0.8 and VIF < 10 criteria. The correlation/VIF workflow audits candidate sets of sizes 16, 15 and 14, records that no 16- or 15-predictor set passes both criteria, and marks the retained 14-predictor set as the feasible 14-predictor set with the lowest max VIF.

Workflow 04 assesses functional form on one common complete-case dataset while keeping the adjustment structure identical between each focal linear-versus-smooth comparison. It uses maximum likelihood for valid AIC comparison and requires all seven prespecified conditions for a nonlinear decision: AIC improvement, material EDF, smooth-term support, adequate basis dimension, distance from the EDF ceiling, stability across basis sizes and robustness to boundary trimming. The current assessment classifies all 14 predictors as linear. `Precipitation of Warmest Quarter (Bio18)` shows partial evidence of curvature but does not satisfy all seven conditions, so it remains linear at this stage.

Workflow 05 does not reassess functional form. It independently applies the predefined contribution-based parsimony rule, removes `Precipitation of Warmest Quarter (Bio18)`, and carries the upstream form metadata for the 13 retained predictors into workflow 06. Workflow 06 then constructs each retained environmental term dynamically from `selected_form` and `selected_k` and fits the final specification with REML. The current final model therefore contains 13 linear environmental predictors, but linearity is a computed input result rather than a hard-coded restriction.

Environmental-variable smooths, random-effect smooths and the spatial smooth have distinct roles. Only retained environmental predictors receive the upstream linear/nonlinear classification. Terrain and host genus are random-effect adjustments, while the two-dimensional longitude/latitude smooth adjusts for residual spatial structure; neither type is evidence that an environmental predictor is nonlinear.

## GAMM Chain Interfaces

- Workflow 04 reads workflow 03's screened standardized modelling table and writes the 14-row model-form assessment, a criterion summary and diagnostic PDF.
- Workflow 05 reads that assessment, the preliminary contribution table and the screened modelling table, then writes the 13-row parsimonious predictor table and a removed-term audit.
- Workflow 06 reads the screened modelling table, the parsimonious predictor table and the complete model-form assessment, then writes the form implementation audit, final model and drop-one summaries, and PDF figures.
