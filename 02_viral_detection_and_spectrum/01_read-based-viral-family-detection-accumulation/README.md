# 01 Read-Based Viral Family Detection Accumulation

## Purpose

Plot cumulative viral-family detection as assigned reads increase within each sample.

## Input Files

- `input/viral_family_read_counts_by_sample.csv`: sample-level viral family read-count table with `sample_id`, `viral_family_tax_id` and `assigned_reads`.

## Scripts

- `scripts/plot_viral_family_detection_accumulation_curve.py`: calculates cumulative detected viral families per sample and writes a PDF curve plot.

## Expected Outputs

- `output/figures/viral_family_detection_accumulation_curve_0_20000_reads.pdf`: read-depth accumulation plot limited to 0-20,000 assigned reads per viral family.

## Method Rationale

The plot uses cumulative viral-family detection across assigned reads to evaluate whether family-level detection changes strongly with read depth. Assigned-read counts are highly right-skewed in the input table: the median is approximately 4,491 reads, the 75th percentile is approximately 24,054 reads, and the maximum exceeds 12 million reads. A fixed 0-20,000 read window therefore focuses the figure on the main low-to-moderate read-depth accumulation region where most curve changes occur, while preventing a small number of extreme high-depth observations from stretching the x-axis and visually compressing the informative part of the curves.

The 20,000-read limit is a display window, not a data-filtering rule for the released input table. The full input table is retained, and the plotting window can be changed in the script when a separate inspection of the long high-depth tail is needed.
