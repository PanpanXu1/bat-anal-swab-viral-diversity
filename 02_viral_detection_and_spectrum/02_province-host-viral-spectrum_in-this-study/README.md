# 02 Province-Host Viral Spectrum (In This Study)

## Purpose

Plot viral-family RPM across province-host combinations for samples in this study.

## Input Files

- `input/viral_spectrum_by_province_host_rpm.csv`: viral-family-by-province-host RPM matrix.
- `input/viral_spectrum_by_province_host_row_metadata.csv`: viral family order, type and display-order metadata.
- `input/viral_spectrum_by_province_host_column_metadata.csv`: province and host-species metadata for each matrix column.

## Scripts

- `scripts/plot_viral_spectrum.py`: generates the province-host viral spectrum bubble heatmap.

## Expected Outputs

- `output/figures/viral_spectrum_by_province_host.pdf`: PDF heatmap; host species and viral family/order labels are rendered in italic where applicable.

## Method Rationale

Province-host combinations are used here because this panel needs to preserve both geographic and host context from the in-study samples. RPM values are used to make viral-family signal less dependent on raw sequencing depth. The metadata files are kept separate from the abundance matrix so that label order, taxonomic grouping and display formatting can be audited without changing the numeric input.
