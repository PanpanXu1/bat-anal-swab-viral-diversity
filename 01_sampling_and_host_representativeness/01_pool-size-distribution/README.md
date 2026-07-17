# 01 Pool Size Distribution

## Purpose

Summarize the number of individual swabs included in each sequencing pool.

## Input Files

- `input/pool_sample_metadata.xlsx`: minimal pool membership table with `Number` and `pool_id`. `Number` is a sequential record identifier from `S1` to `S5498`; `pool_id` identifies the sequencing pool.

## Scripts

- `scripts/plot_pool_sample_count_distribution.py`: counts records within each `pool_id` and plots the pool-size distribution.

## Expected Outputs

- `output/figures/pool_sample_count_distribution.pdf`: PDF bar plot showing the number of pools within each individual-swab count bin.

## Method Rationale

The input table keeps only `Number` and `pool_id` because this workflow only needs the record order and the pool grouping field. The `Number` field is a serial label rather than a biological sample identifier. Pool-size counts are calculated directly from `pool_id`.
