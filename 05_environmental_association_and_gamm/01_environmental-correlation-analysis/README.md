# 01 Environmental Correlation Analysis

## Purpose

Summarize environmental correlations and Mantel-style association checks.

## Input Files

- `input/shannon_environmental_analysis_input.csv`

## Script Workflow

1. `scripts/01_calculate_environmental_spearman_correlations.py`
   - Reads the environmental analysis input table.
   - Calculates the environmental-variable Spearman correlation matrix.
   - Writes both matrix-format and long-format pairwise correlation tables.

2. `scripts/02_run_shannon_environment_mantel_tests.py`
   - Uses the same input table to evaluate Shannon-environment associations with permutation-based Mantel-style tests.
   - Writes one row per environmental variable, including the test statistic, empirical P value, permutation count and random seed.

3. `scripts/03_plot_environmental_correlation_mantel.py`
   - Reads the outputs from Steps 1 and 2.
   - Draws the combined correlation/Mantel-link PDF figure.
   - This script depends on the two upstream output folders and should be run after Steps 1 and 2.

## Expected Outputs

### Figures

1. `output/figures/environmental_correlation_mantel_plot/environmental_correlation_mantel_plot.pdf`
   - Combined environmental correlation matrix and Shannon-environment Mantel-link plot.
   - The matrix portion is derived from `output/tables/environmental_spearman/environmental_spearman_correlation_matrix.csv`; the Mantel links are derived from `output/tables/shannon_environment_mantel/shannon_environment_mantel_tests.csv`.

### Tables

1. `output/tables/environmental_spearman/environmental_spearman_correlation_matrix.csv`
   - Square matrix of pairwise Spearman rho values among environmental variables.
   - Rows and columns are environmental predictors.
2. `output/tables/environmental_spearman/environmental_spearman_pairwise_tests.csv`
   - Long-format correlation table.
   - Key columns: `Factor1`, `Factor2`, `Spearman_rho` and `p_value`.
3. `output/tables/shannon_environment_mantel/shannon_environment_mantel_tests.csv`
   - Mantel-style association checks between Shannon diversity and each environmental variable.
   - Key columns: `Environmental_variable`, `Mantel_r`, `p_value`, `permutations` and `random_seed`.

## Notes

Run the scripts in numeric order. Figure outputs are standardized to PDF.

## Method Rationale

Spearman correlations are used for environmental predictors because they summarize monotonic relationships and are less sensitive to non-normal predictor distributions than Pearson correlations. Mantel-style tests are used for Shannon-environment association checks because the workflow evaluates distance/similarity structure rather than only independent scalar correlations. Permutation-derived P values are recorded with the permutation count and random seed so that empirical significance calculations are reproducible.

Only Mantel links at or below the configured plotting threshold are drawn in the figure to keep the visual layer interpretable. The complete numerical test table is still exported, so omitted links are a display decision rather than missing analysis output.
