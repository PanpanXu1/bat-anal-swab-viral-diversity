# 02 Viral Detection and Spectrum

## Purpose

Reproduce the code-supported viral detection accumulation and viral spectrum visualizations.

## Workflows

| Order | Workflow | Purpose |
|---|---|---|
| 1 | `01_read-based-viral-family-detection-accumulation/` | Plots cumulative viral-family detection against assigned reads. |
| 2 | `02_province-host-viral-spectrum_in-this-study/` | Plots province-host viral spectrum for samples in this study. |
| 3 | `03_host-genus-viral-spectrum_in-this-study/` | Plots host-genus viral spectrum for samples in this study. |
| 4 | `04_host-genus-viral-spectrum_previous-datasets/` | Plots host-genus viral spectrum for previous datasets using the same visualization logic. |

## Relationship Between In-Study and Previous-Dataset Spectrum Analyses

The host-genus spectrum analyses for samples in this study and for previous datasets use comparable matrix-based visualization to examine whether the same host-genus signal is observed across independent data sources. This pairing is included to support reproducible comparison of overall viral loads across viral families, with particular attention to the *Rhinolophus* host signal relevant to coronavirus-related analyses. The previous-dataset workflow is treated as an independent comparison layer using the same plotting contract, rather than as a duplicate of the in-study samples.

## Method Rationale

The three viral-spectrum workflows use parallel input structures and plotting rules so that differences between panels reflect the data source and grouping level rather than ad hoc figure construction. RPM and log10 RPM matrices are used to reduce sensitivity to sequencing-depth differences and right-skewed abundance values. Row-standardized heatmaps are used where the aim is to compare relative host or host-genus patterns within each viral family, rather than to let a few high-abundance families dominate the color scale across the full matrix.
