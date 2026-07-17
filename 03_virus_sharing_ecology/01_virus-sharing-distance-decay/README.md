# 01 Virus Sharing Distance Decay

## Purpose

Evaluate the relationship between geographic distance and shared viral-cluster counts across sampling groups.

## Input Files

- `input/distance_decay_viral_cluster_membership.csv`: viral-cluster membership table used to count shared viral clusters between sampling groups.
- `input/sampling_site_pair_distances.csv`: pairwise sampling-site distance table used to attach geographic distances to the shared-cluster comparisons.

## Scripts

- `scripts/analyze_virus_sharing_distance_decay.py`: validates the pairwise input table, keeps site pairs with non-zero shared viral clusters, calculates log10(shared clusters), runs the Spearman association test and writes the distance-decay PDF.

## Expected Outputs

### Figures

1. `output/figures/distance_decay_log10.pdf` - scatterplot of geographic distance against log10(shared clusters) for site pairs with `SharedClusters > 0`, with the fitted visual trend used for checking distance-decay structure.

### Tables

1. `output/tables/shared_clusters_distance.csv`
   - Pairwise analysis table used for plotting and testing.
   - Key columns: `SampleGroup1`, `SampleGroup2`, `SharedClusters`, `Distance_km` and `Log10SharedClusters`.
2. `output/tables/spearman_correlation.csv`
   - Spearman test summary for distance and viral-cluster sharing.
   - Key columns: `N_site_pairs`, `Spearman_rho`, `P_value` and `Included_pairs`.

## Method Rationale

Only site pairs with `SharedClusters > 0` are included in this distance-decay analysis, and those positive counts are transformed as log10(shared clusters). Zero-sharing site pairs are excluded before transformation because log10(0) is undefined and the workflow tests whether the magnitude of observed sharing declines with distance among pairs that share at least one viral cluster. Spearman correlation is used because the distance-decay question is rank-based and does not require a linear relationship or normally distributed residuals.
