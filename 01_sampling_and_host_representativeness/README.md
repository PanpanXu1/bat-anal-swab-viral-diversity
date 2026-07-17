# 01 Sampling and Host Representativeness

This block corresponds to the manuscript's introductory sampling representativeness analyses. It groups workflows that describe sampling scale, host occurrence coverage and richness representation before the manuscript moves into viral detections.

## Workflow Folders

| Order | Workflow | Purpose |
|---|---|---|
| 1 | `01_pool-size-distribution/` | Summarizes individual swabs per sequencing pool. |
| 2 | `02_gbif-host-diversity-comparison/` | Compares in-house host occurrence coverage with curated Global Biodiversity Information Facility (GBIF) records. |
| 3 | `03_seasonal-species-accumulation-extrapolation/` | Plots seasonal species accumulation and iNEXT-based extrapolation curves used to support sampling representativeness. |

## Method Rationale

This block is placed before viral analyses because the sampling design and host representativeness checks define the population context for downstream viral summaries. The workflows intentionally separate pool-size summaries, GBIF comparison and species-accumulation analyses because each addresses a different source of sampling concern: sequencing-pool composition, independent host-occurrence coverage and richness representation under observed sampling effort. Each workflow writes generated products to its own `output/` directory so the provenance of figures and tables remains local to the input and script that generated them.
