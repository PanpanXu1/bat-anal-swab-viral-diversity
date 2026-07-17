# 01 Adjusted Shannon Model

## Purpose

Fit the adjusted Shannon diversity model for terrain-level comparison while accounting for sampling and host-composition covariates.

## Input Files

- `input/adjusted_shannon_model_input.csv`

## Scripts

- `scripts/fit_adjusted_shannon_terrain_model.R`: fits the terrain model, estimates adjusted marginal means under the configured weighting approaches, exports pairwise contrasts and writes PDF summaries.

## Expected Outputs

### Figures

1. `output/figures/primary_adjusted_predicted_Shannon_equal_weights.pdf` - adjusted predicted Shannon diversity by terrain using equal weighting of the covariate distribution.
2. `output/figures/primary_adjusted_predicted_Shannon_proportional_weights.pdf` - adjusted predicted Shannon diversity by terrain using proportional weighting of the observed covariate distribution.

### Tables

1. `output/tables/primary_model_summary.txt`
   - Text summary of the fitted primary model, including model terms and diagnostics reported by R.
2. `output/tables/primary_adjusted_predicted_Shannon_by_terrain.csv`
   - Adjusted marginal means by terrain and weighting approach.
   - Key columns: `weighting`, `terrain`, `emmean`, `SE`, confidence limits and model test statistics.
3. `output/tables/mountain_vs_other_terrain_equal_weights_Tukey.csv`
   - Tukey-adjusted contrasts comparing mountain terrain with other terrain categories under the equal-weighting setting.

## Notes

Run the script from this workflow folder or pass the script path directly to R. Figure outputs are standardized to PDF.

## Method Rationale

Adjusted marginal means are used because terrain comparisons can be affected by sampling size, season, year and host-composition covariates. Reporting both equal-weight and proportional-weight summaries makes the weighting assumption explicit: equal weighting compares terrains under a balanced covariate distribution, whereas proportional weighting preserves the observed covariate distribution. The contrast table is exported separately so that the plotted adjusted means and the pairwise comparison values can be audited independently.
