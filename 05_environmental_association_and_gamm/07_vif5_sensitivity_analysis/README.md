# VIF < 5 sensitivity GAMM

This workflow applies a VIF < 5 sensitivity analysis to the 13 environmental predictors in the primary GAMM (workflow 06). At each selection step the highest-VIF predictor was removed and VIFs were recalculated. Bio17, HFT and GMR were removed, leaving ten predictors with a maximum VIF of 4.85. The primary GAMM is provided by workflow 06.

## Files

- `input/environmental_variables_VIF5.csv`: the same 225 pool-level observations and non-environmental fields used by workflow 06, with only the ten retained environmental columns.
- `input/parsimonious_environmental_predictors.csv`: the ten predictors and their inherited model-form decisions.
- `input/environmental_predictor_model_form_assessment.csv`: workflow 04's original 14-row form assessment, required by the model-form validation in the R script. Only the ten predictors above enter this fit.
- `scripts/fit_parsimonious_gamm.R`: workflow 06's GAMM code with only the modelling-input filename changed to `environmental_variables_VIF5.csv`.
- `scripts/run_vif5_gamm.py`: reconstructs the VIF < 5 deletion sequence from workflow 06's saved 13 predictors, verifies the ten-predictor input, runs the R fit and formats its term table.
- `output/tables/model_fit_summary.csv`: unrounded term-level results from the refitted GAMM.
- `output/tables/gamm_vif5_term_estimates.csv`: the nine-column term table, aligned with the primary GAMM table. Estimates and test statistics are displayed to two decimals; P values near 0.05 or below 0.01 retain additional precision.
- `output/tables/environmental_predictor_vif.csv`: individual VIF values, rounded to two decimals, for the ten environmental predictors in the VIF < 5 GAMM. Pool size, season, year and land-use category are adjustment covariates outside this environmental VIF screen.

## Run

From `05_environmental_association_and_gamm/`:

```bash
python 07_vif5_sensitivity_analysis/scripts/run_vif5_gamm.py
```

Python requires NumPy and pandas (listed in the repository-root `requirements.txt`); R requires `mgcv` (`r_requirements.txt`). If `Rscript` is not on `PATH`, pass `--rscript` followed by its executable path.

The refitted model uses 225 pools, explains 42.56% of deviance and has adjusted R² = 0.33. Its environmental terms are Bio2, Bio4, Bio5, Bio8, CHEQ, NDVI, FVC, PSD, GR and GLH. The analysis starts from the primary model's 13 environmental predictors and keeps its adjustment structure unchanged.
