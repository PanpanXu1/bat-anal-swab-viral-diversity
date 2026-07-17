# 02 Environmental Variable Standardization

## Purpose

Standardize environmental predictors for downstream screening and modelling.

## Input Files

- `input/environmental_preprocessing_input.csv`

## Scripts

- `scripts/standardize_environmental_variables.py`: reads the pool-level environmental input table, keeps modelling metadata fields unchanged, and standardizes numeric environmental predictors for downstream correlation/VIF screening and GAMM fitting.

## Expected Outputs

### Tables

1. `output/tables/standardized_environmental_variables_by_pool.csv`
   - Pool-level modelling table with standardized environmental predictors.
   - Metadata columns such as `Sample group`, `Shannon index`, `Season`, `Year`, `SampleSize (log10)`, `latitude`, `longitude`, `Terrain` and `Host genus` are retained to support downstream models.
   - Environmental predictor columns are standardized so that later correlation, VIF and GAMM steps use comparable numeric scales.

## Notes

This output is copied into the next module as an input table so the VIF-screening workflow can be run independently.

## Method Rationale

Environmental predictors are standardized before screening and modelling because they are measured on different scales. Standardization makes correlation screening, VIF evaluation and model fitting less sensitive to arbitrary measurement units. Metadata and response fields are kept unstandardized where they are identifiers, grouping variables or model covariates that should retain their original interpretation.
