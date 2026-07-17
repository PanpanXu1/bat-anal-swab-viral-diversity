# 02 Terrain Virus Sharing

## Purpose

Analyze whether viral-cluster sharing differs across terrain categories.

## Input Files

- `input/terrain_viral_cluster_membership.csv`: viral-cluster membership table with sample group and terrain assignments.

## Script Workflow

1. `scripts/01_plot_terrain_shared_cluster_heatmap.py`
   - Converts sample-level viral-cluster membership into a sample-pool-by-sample-pool shared-cluster matrix ordered by terrain.
   - Writes the heatmap used to visually inspect sample-pool sharing patterns with terrain annotations.
   - This script provides the descriptive matrix view; it does not run the formal terrain association tests.

2. `scripts/02_test_terrain_association_with_virus_sharing.py`
   - Converts the same membership table into unique non-self sample-pair records.
   - Writes the pairwise table used as the statistical testing input.
   - Runs ANOVA-style and Kruskal-Wallis-style permutation tests to evaluate terrain-associated differences in shared-cluster counts.

## Expected Outputs

### Figures

1. `output/figures/terrain_shared_cluster_heatmap.pdf` - sample-pool-by-sample-pool heatmap of log10(shared viral clusters + 1). Rows and columns are sample pools ordered by terrain; terrain color bars annotate each sample pool.

### Tables

1. `output/tables/terrain_pairwise_shared_clusters.csv`
   - Pairwise sample table used for statistical testing.
   - Key columns: `SampleGroup1`, `SampleGroup2`, `Terrain1`, `Terrain2`, `TerrainPair` and `SharedClusters`.
2. `output/tables/terrain_anova_results.csv`
   - Permutation-based ANOVA-style test summary.
   - Key columns: `F_statistic`, `P_value`, `Permutations` and `Random_seed`.
3. `output/tables/terrain_kruskal_wallis_results.csv`
   - Permutation-based Kruskal-Wallis-style test summary.
   - Key columns: `H_statistic`, `P_value`, `Permutations` and `Random_seed`.

## Method Notes

The heatmap and the formal tests are separated because they serve different purposes: the heatmap gives a transparent terrain-annotated sample-pool matrix, whereas the test script evaluates whether pairwise sharing differs by terrain grouping. Both scripts use the same membership input, which keeps the descriptive and testing outputs traceable to one source table.

Permutation tests are used because shared-cluster counts are discrete, pairwise and not guaranteed to satisfy parametric distribution assumptions. The terrain-pair test uses 999 permutations by default as a lightweight reproducibility setting. Increase `PERMUTATIONS` in the script for a higher-resolution empirical P value.
