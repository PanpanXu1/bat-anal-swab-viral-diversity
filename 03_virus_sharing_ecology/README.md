# 03 Virus Sharing Ecology

This block corresponds to the manuscript section on ecological and host-associated patterns in viral sharing.

## Workflow Folders

| Order | Workflow | Purpose |
|---|---|---|
| 1 | `01_virus-sharing-distance-decay/` | Evaluates geographic distance decay among site pairs with non-zero shared viral clusters. |
| 2 | `02_terrain-virus-sharing/` | Plots a terrain-annotated sample-pool sharing matrix and tests sharing across terrain-pair categories. |
| 3 | `03_host-identity-virus-sharing/` | Relates host sequence identity to shared viral clusters. |
| 4 | `04_host-taxonomy-virus-sharing/` | Tests host-taxonomic grouping and viral-cluster sharing. |

Coordinate-bearing inputs used for spatial distance context are rounded to two decimal places for sensitive locality protection.

## Method Rationale

The workflows are ordered from geography, to landscape category, to host genetic identity, and then host taxonomy. This mirrors a progression from spatial separation to ecological grouping and then biological host relatedness. Pairwise tables are retained in the submodules because viral sharing is naturally defined between sites or hosts, and those pairwise records provide the audit trail for each statistical summary. Transformations are documented within each workflow because the distance-decay analysis uses positive shared-cluster pairs, whereas the matrix-style host-sharing analyses retain zero-sharing pairs with log10(x + 1).
