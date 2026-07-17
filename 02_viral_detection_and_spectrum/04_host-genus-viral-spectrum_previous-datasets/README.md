# 04 Host-Genus Viral Spectrum (Previous Datasets)

## Purpose

Plot host-genus viral spectrum for previous datasets using the same visualization logic as the in-study host-genus workflow.

## Input Files

- `input/viral_spectrum_previous_datasets_by_host_genus_log10_rpm.csv`: viral-family-by-host-genus log10 RPM matrix from previous datasets.
- `input/viral_spectrum_previous_datasets_by_host_genus_row_metadata.csv`: viral family type metadata.
- `input/viral_spectrum_previous_datasets_by_host_genus_pool_counts.csv`: host-genus pool counts used as supporting metadata.

## Scripts

- `scripts/plot_viral_spectrum.py`: clusters host genera and plots the row-standardized viral spectrum heatmap.

## Expected Outputs

- `output/figures/viral_spectrum_previous_datasets_by_host_genus.pdf`: PDF heatmap for previous datasets.

## Notes

This workflow is an independent comparison layer. It uses a plotting structure matched to the in-study host-genus workflow so that the *Rhinolophus* host-genus signal can be assessed across independent datasets without mixing previous records into the in-study sample matrix.

## Method Rationale

The previous-datasets matrix is kept in a separate workflow to avoid combining independent data sources into the in-study sample matrix. Using the same host-genus grouping, log10 RPM transformation and row-standardized heatmap logic makes the comparison method consistent while preserving dataset provenance. This design supports a reproducible cross-dataset check without treating previous records as additional samples from this study.
