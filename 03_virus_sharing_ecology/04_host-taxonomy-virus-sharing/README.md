# 04 Host Taxonomy Virus Sharing

## Purpose

Test whether viral-cluster sharing is higher within the same host family or genus than between different families or genera.

## Input Files

- `input/host_taxonomy.csv`: host species with family and genus assignments. Species names use spaces between genus and species.
- `input/host_taxonomy_viral_cluster_membership.csv`: viral-cluster membership by host species.

## Scripts

- `scripts/analyze_host_taxonomy_virus_sharing.py`: constructs host-pair sharing records and runs one-sided Mann-Whitney U tests for same-family and same-genus comparisons.

## Expected Outputs

### Tables

1. `output/tables/host_pair_taxonomy_and_sharing.csv`
   - Long-format host-pair table used as the statistical-test input.
   - Key columns include `Family1`, `Family2`, `Genus1`, `Genus2`, `Same_family`, `Same_genus`, `Shared_clusters` and `Log10_shared_clusters_plus1`.
2. `output/tables/mann_whitney_u_results.csv`
   - Mann-Whitney U test summary for same-family and same-genus comparisons.
   - Key columns include the taxonomic level, alternative hypothesis, group sizes, U statistic, P value and mean/median log-transformed sharing values.

## Output Relationship

The pairwise table is the transparent intermediate table. The Mann-Whitney result table is calculated directly from the `Same_family`, `Same_genus` and `Log10_shared_clusters_plus1` fields in that pairwise table.

## Method Rationale

Family and genus comparisons are performed on host-pair records because the biological question is whether sharing differs for within-taxon versus between-taxon host pairs. The Mann-Whitney U test is used because the transformed sharing values are not assumed to be normally distributed and the comparison is between two pair categories. The log10(shared clusters + 1) transformation keeps zero-sharing pairs and reduces count-scale skew before group comparison.
