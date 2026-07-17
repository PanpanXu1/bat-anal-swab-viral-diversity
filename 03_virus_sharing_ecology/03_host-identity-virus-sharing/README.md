# 03 Host Identity Virus Sharing

## Purpose

Compare host sequence identity with viral-cluster sharing across host species.

## Input Files

- `input/host_identity_viral_cluster_membership.csv`: viral-cluster membership by host species.
- `input/host_sequence_identity_percent.csv`: symmetric host sequence identity matrix. Host species names use spaces between genus and species.

## Scripts

- `scripts/analyze_host_identity_virus_sharing.py`: calculates shared-cluster matrices, transforms sharing as log10(shared clusters + 1), runs a two-sided Spearman Mantel test and writes heatmaps.

## Expected Outputs

### Figures

1. `output/figures/host_sequence_identity_heatmap.pdf` - heatmap of host sequence identity percentages. The same host order is used in the sharing heatmap to support direct visual comparison.
2. `output/figures/host_shared_clusters_heatmap.pdf` - heatmap of log10(shared viral clusters + 1) values across host species.

### Tables

1. `output/tables/host_shared_cluster_matrix.csv`
   - Symmetric matrix of raw shared viral-cluster counts by host pair.
   - Host names are stored with spaces between genus and species for readable table display.
2. `output/tables/host_shared_cluster_log10_plus1_matrix.csv`
   - Symmetric matrix after log10(shared clusters + 1) transformation.
   - This transformed matrix is used for the sharing heatmap and Mantel-style comparison.
3. `output/tables/host_pair_identity_and_sharing.csv`
   - Long-format host-pair table linking `Host_identity_percent`, `Shared_clusters` and `Log10_shared_clusters_plus1`.
   - This table is the most direct audit table for the relationship between sequence identity and sharing.
4. `output/tables/mantel_test_results.csv`
   - Test summary recording the Mantel statistic, P value, permutation count, random seed and sharing transformation.

## Method Rationale

The analysis uses pairwise host matrices because both host sequence identity and viral sharing are naturally defined between host pairs. Viral sharing is transformed as log10(shared clusters + 1) to keep zero-sharing pairs while reducing the influence of very large counts. A Mantel-style permutation test is used because the observations are pairwise distances/similarities rather than independent single-host records; permutation testing provides an empirical significance calculation under the matrix structure.
