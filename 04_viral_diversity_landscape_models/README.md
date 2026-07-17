# 04 Viral Diversity and Landscape Models

This block corresponds to the manuscript analyses of pool-level viral diversity after the earlier sampling and viral-spectrum sections.

## Workflow Folders

| Order | Workflow | Purpose |
|---|---|---|
| 1 | `01_adjusted-shannon-model/` | Fits the adjusted Shannon diversity model using terrain, pool size and temporal covariates. |
| 2 | `02_spatial-autocorrelation/` | Calculates global Moran's I for pool-level Shannon diversity. |

The adjusted model controls for pool size and temporal sampling variables to reduce sampling-effort bias. Coordinate-bearing inputs are rounded to two decimal places before public release.

## Method Rationale

The adjusted Shannon model and spatial autocorrelation check are kept as separate workflows because they answer different reliability questions. The adjusted model evaluates terrain-associated diversity after accounting for sampling and temporal covariates, while the Moran's I workflow checks whether residual spatial structure may be relevant at the released coordinate resolution. Keeping them separate avoids mixing a model-estimation step with a spatial diagnostic.
