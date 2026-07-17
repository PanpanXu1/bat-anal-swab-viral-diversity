# 03 Host-Genus Viral Spectrum (In This Study)

## Purpose

Plot host-genus viral spectrum for samples in this study using row-standardized viral-family RPM values.

## Input Files

- `input/viral_spectrum_by_host_genus_log10_rpm.csv`: viral-family-by-host-genus log10 RPM matrix.
- `input/viral_spectrum_by_host_genus_row_metadata.csv`: viral family type metadata.
- `input/viral_spectrum_by_host_genus_pool_counts.csv`: host-genus pool counts used as supporting metadata.

## Scripts

- `scripts/plot_viral_spectrum.py`: clusters host genera and plots the row-standardized viral spectrum heatmap.

## Expected Outputs

- `output/figures/viral_spectrum_by_host_genus.pdf`: PDF heatmap comparing viral spectrum across host genera in this study.

## Notes

This workflow is paired with `04_host-genus-viral-spectrum_previous-datasets/` to compare whether host-genus viral-load patterns are consistent across independent data sources.

## Method Rationale

Host genus is used as the grouping level to make the matrix comparable across samples with uneven species coverage while retaining a biologically interpretable host-taxonomic level. The input values are log10 RPM to reduce the influence of highly skewed viral-family abundance values. Row standardization is applied before plotting because the visualization compares host-genus patterns within each viral family; it is not intended to rank viral families by absolute abundance across the whole matrix.
