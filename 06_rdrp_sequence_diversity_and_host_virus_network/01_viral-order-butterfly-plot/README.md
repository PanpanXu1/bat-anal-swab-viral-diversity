# 01 Viral Order Butterfly Plot

## Purpose

Visualize viral-order abundance patterns between the NIPB and ZOVER groups using a butterfly-style bar plot.

## Input Files

- `input/viral_order_log10_abundance_by_pool.csv`: viral-order table with log10(counts + 1) values for NIPB and ZOVER.

## Scripts

- `scripts/plot_viral_order_butterfly.py`: validates the viral-order abundance table and writes the butterfly plot as a PDF.

## Expected Outputs

### Figures

1. `output/figures/viral_order_butterfly_plot.pdf` - horizontal butterfly plot comparing NIPB and ZOVER viral-order abundance values. The inset summarizes sequence proportions recovered from the supplied log10(counts + 1) values.

## Method Rationale

The butterfly layout is used because the two groups share the same viral-order axis and can be compared symmetrically without creating separate panels. Input values are log10(counts + 1) to keep zero counts valid and reduce the influence of highly abundant orders. The inset sequence proportion is recalculated from the supplied log-transformed values so that it remains tied to the released input table.
