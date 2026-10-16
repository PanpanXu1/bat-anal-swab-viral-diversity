# 01 Virus Sharing Distance Decay

## Purpose

Evaluate whether shared viral-cluster counts decline with geographic distance among sampling-group pairs sharing at least one viral cluster. Compare the complete positive-sharing subset with subsets separated by >0, >25, >50 and >100 km.

## Input Files

- `input/distance_decay_viral_cluster_membership.csv`: sequence identifiers, viral-cluster identifiers and sampling-group identifiers used to count distinct shared clusters.
- `input/sampling_site_pair_distances.csv`: distances in kilometres for all 26,335 unordered pairs of the 230 sampling groups. Columns are `SampleGroup1`, `SampleGroup2` and `Distance_km`; self-pairs are excluded. Distances were calculated from the original coordinates using GeoPy's WGS84 geodesic distance and retained to three decimal places. The 7,871 positive-sharing pair distances used in the observed analysis are included in this table.

The complete distance table supports label permutations. Observed correlations and scatterplots include only pairs with at least one shared viral cluster, regardless of the number of rows in the distance input.

## Coordinate Protection and Reproducibility

Public coordinate fields elsewhere in this repository are rounded to two decimal places to protect sensitive cave roost locations. Precise coordinates are not included in this workflow. The supplied pairwise distances were calculated before coordinate rounding; the scripts read this distance table directly. Recalculating distances from publicly rounded coordinates will not reproduce the supplied distances and may change which pairs have zero distance or fall within a distance threshold.

For exact reproduction of this workflow, use the included membership and distance tables. No private coordinate file or local source path is required to run either script.

## Installation and Execution

Use Python 3.10 or later. From this folder, run:

```bash
python -m pip install -r requirements.txt
python scripts/analyze_virus_sharing_distance_decay.py
python scripts/analyze_distance_threshold_sensitivity.py
```

The first script generates the analysis without a distance cutoff. The second generates the five-row sensitivity summary and all five distance-decay figures. Both use the same permutation implementation, seed and ordering of sampling groups.

## Method

Each sampling group's viral clusters are represented as a set. The number of shared clusters is the size of the intersection between two sets. Spearman's rank correlation coefficient is calculated between geographic distance and shared-cluster counts, restricted to positive-sharing pairs. Average ranks are used for ties. Sampling groups can share a geographic location; the analysis without a distance cutoff includes these 0 km pairs.

Two-sided P values use 9,999 unrestricted sampling-group label permutations, with NumPy random seed 42. The complete shared-cluster matrix includes zero-sharing entries. For each permutation, its rows and columns are reordered using the same random sampling-group index vector, while the geographic distance matrix remains fixed. Only unique off-diagonal pairs are considered. The positive-sharing filter and the applicable geographic distance threshold are applied anew, and Spearman ranks are recalculated for each resulting subset. Thus, threshold-specific subset sizes can vary across permutations; the reported pair counts describe the observed data.

For each distance subset, the P value is `(b + 1) / (9999 + 1)`, where `b` is the number of permutations with an absolute correlation at least as large as the observed absolute correlation. The same sequence of sampling-group permutations is used across all five subsets.

### Permutation Hypothesis

The null hypothesis is that the correspondence between sampling-group identities in the viral-sharing matrix and their geographic labels is exchangeable. Simultaneous row and column permutation preserves the complete sharing matrix's internal pair structure while changing its alignment with geography. This is a sampling-group label randomization test of positive-sharing correlations, rather than a standard Mantel test using every pair in its correlation statistic. It does not estimate a geographic effect adjusted for host identity, season or other sampling-group attributes.

### Distance Thresholds and Figures

The observed analysis is repeated for All positive-sharing pairs and for distances >0, >25, >50 and >100 km. The >50 km subset evaluates whether the association extends beyond nearby sampling groups.

Plots show geographic distance against log10(shared-cluster counts), with a fitted linear trend for visualization. The correlation statistic uses the untransformed positive counts; log10 preserves their rank order. Each figure displays Spearman's rank correlation coefficient, its permutation P value and the distance filter. The >0 km figure uses the minimum included positive distance as its first x-axis tick.

## Expected Outputs

### Figures

PDF figures are written to `output/figures/`:

1. `distance_decay_log10.pdf`: all positive-sharing pairs, including 0 km pairs.
2. `distance_decay_over_0km.pdf`: distance >0 km.
3. `distance_decay_over_25km.pdf`: distance >25 km.
4. `distance_decay_over_50km.pdf`: distance >50 km.
5. `distance_decay_over_100km.pdf`: distance >100 km.

### Tables

1. `output/tables/shared_clusters_distance.csv`: observed positive-sharing pair table, including the two sampling groups, shared-cluster count, distance and log10(shared-cluster count).
2. `output/tables/spearman_correlation.csv`: correlation and label-permutation summary without a distance cutoff.
3. `output/tables/distance_threshold_sensitivity.csv`: corresponding summary for All, >0, >25, >50 and >100 km.

Summary tables use full column names, including `Number of sampling-group pairs`, `Spearman's rank correlation coefficient`, `P value`, `Number of extreme permutations`, `Number of permutations` and `Permutation unit`.

Displayed correlations, axis ticks and numeric pair-table outputs use two decimal places; counts and distance cutoffs are integers. Calculations use the input precision. P values below 0.001 are displayed as `<0.001`; values from 0.001 to below 0.01 are displayed as `<0.01` to avoid a misleading `0.00`. Other P values use two decimal places. The extreme-permutation count and total number of permutations retain the information needed to reconstruct each numerical P value. With 9,999 permutations, zero extreme permutations gives a numerical P value of 0.0001, not zero.

## Repository Contents

The `input/` directory contains the two analysis tables, `scripts/` contains the two analysis scripts, and `requirements.txt` records the runtime dependencies. Generated outputs, temporary files and Python caches are excluded by `.gitignore`.
