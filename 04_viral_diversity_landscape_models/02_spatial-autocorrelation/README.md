# 02 Spatial Autocorrelation

## Purpose

Test whether pool-level Shannon diversity has global spatial autocorrelation under a k-nearest-neighbor spatial weights definition.

## Input Files

- `input/viral_shannon_coordinates.csv`: pool-level Shannon diversity and rounded coordinate fields.

## Scripts

- `scripts/calculate_global_moran_spatial_autocorrelation.py`: builds k-nearest-neighbor weights, calculates global Moran's I and estimates an empirical P value using permutation testing.

## Expected Outputs

### Tables

1. `output/tables/global_moran_summary.csv`
   - One-row summary of the global Moran's I test.
   - Key columns: `n_samples`, `k_neighbors`, `permutations`, `random_seed`, `global_moran_i`, `expected_moran_i`, `permutation_p_value` and `permutation_p_value_display`.

## Method Notes

Global Moran's I is used as a screening test for broad spatial autocorrelation in pool-level Shannon diversity. A k-nearest-neighbor weight definition is used because it gives each point the same number of spatial neighbors, which is useful when sampling locations are unevenly distributed. Coordinates in the released input are rounded for locality protection, so this check should be interpreted as a reproducible landscape-scale diagnostic rather than a precise cave-location analysis.

When no permuted statistic is at least as extreme as the observed statistic, the empirical P value is reported at the minimum resolvable value for the configured permutation count. With 9999 permutations, this lower bound is `0.0001`, and the display field records `<0.0001`. Reporting a lower bound avoids writing a literal zero P value, which would overstate the numerical precision of a permutation test.
