# Reproducible Analyses of Bat Anal-Swab Viral Diversity, Host-Virus Sharing, Environmental Associations, and Host-Level Coronavirus Surveillance Prioritization

This repository provides reproducible code and analysis-ready data for bat anal-swab viral diversity, host-virus sharing, environmental associations, and host-level prioritization for bat coronavirus surveillance. Each numbered top-level folder represents one analysis block, and each subfolder is an individual workflow with its own `scripts/`, `input/` and `output/` directories.

Scripts write generated figures, tables and diagnostics to workflow-level `output/` directories. When a small generated table is required by a downstream workflow, it is copied into that workflow's `input/` directory and its provenance is described in the workflow README.

## Installation

Clone or download the repository, then create the software environment from the repository root.

```bash
conda env create -f environment.yml
conda activate manuscript_analysis_reproducibility
```

On the tested Ubuntu system described below, a fresh installation from `environment.yml` with an empty package cache took **approximately 11.5 min** wall-clock time. This duration is based on the `real` output from `/usr/bin/time -p`; download and environment-solving time will vary with network speed and hardware.

Python-only workflows can also be run with Python 3.10 or later by installing the root dependencies together with any workflow-specific dependencies listed in the corresponding workflow README or `requirements.txt`.

```bash
python -m pip install -r requirements.txt
```

R workflows require R 4.3 or later and the packages listed in `environment.yml` and workflow-level R requirement files. With a standard R installation, missing R packages can be installed from CRAN, for example:

```r
install.packages(c("iNEXT", "mgcv", "emmeans", "ggplot2", "readr", "dplyr", "tidyr", "tibble"))
```

## Tested Environment

All released analysis scripts were tested successfully in one fresh Conda environment on 22 August 2026.

| Component | Tested value |
|---|---|
| Operating system | Ubuntu 22.04.2 LTS (Jammy Jellyfish), running under WSL2 |
| Architecture | x86_64 |
| Python | 3.10.21 |
| R | 4.5.3 |
| Conda | 23.1.0 |
| CPU | 12th Gen Intel(R) Core(TM) i7-12700; 20 logical CPUs available to WSL2 |
| RAM | 47 GiB available to Ubuntu |

The tested direct Python dependencies were Matplotlib 3.10.9, NumPy 2.2.6, openpyxl 3.1.5, pandas 2.3.3 and ReportLab 5.0.1. The tested direct R dependencies were dplyr 1.2.1, emmeans 2.0.4, ggplot2 4.0.3, iNEXT 3.0.2, mgcv 1.9-4, readr 2.2.0, tibble 3.3.1 and tidyr 1.3.2. These versions are the resolved versions from the successful fresh `environment.yml` installation; the environment file retains compatible constraints rather than an exact platform lock.

No GPU or other non-standard hardware is required. All 30 released analysis scripts completed successfully using CPU execution. Although the tested system exposed a CUDA-capable host to WSL2, none of the released scripts imports, configures or uses GPU software.

## Running the Workflows

Each workflow is self-contained and includes analysis-ready input data, executable scripts, workflow-specific instructions, and expected outputs. The included analysis-ready input data can be used to run and test the corresponding workflow.

Each workflow can be run independently by following its workflow-level README.

To rerun all released Python workflows from the repository root:

```bash
python reproducibility/run_all_python_workflows.py
```

To rerun all released R workflows from the repository root, make sure `Rscript` is on `PATH`, then run:

```bash
python reproducibility/run_all_r_workflows.py
```

To run one workflow independently from the repository root, use the exact script path listed in its workflow-level README:

```bash
python path/to/workflow/scripts/script.py
Rscript path/to/workflow/scripts/script.R
```

Scripts read the analysis-ready files in their workflow-level `input/` directory and create or refresh files under the corresponding `output/` directory. Each workflow README lists the expected output filenames, formats and interpretation.

## Repository Structure

```text
repository_root/
|-- 01_sampling_and_host_representativeness/
|   |-- 01_pool-size-distribution/
|   |-- 02_gbif-host-diversity-comparison/
|   `-- 03_seasonal-species-accumulation-extrapolation/
|-- 02_viral_detection_and_spectrum/
|   |-- 01_read-based-viral-family-detection-accumulation/
|   |-- 02_province-host-viral-spectrum_in-this-study/
|   |-- 03_host-genus-viral-spectrum_in-this-study/
|   `-- 04_host-genus-viral-spectrum_previous-datasets/
|-- 03_virus_sharing_ecology/
|   |-- 01_virus-sharing-distance-decay/
|   |-- 02_terrain-virus-sharing/
|   |-- 03_host-identity-virus-sharing/
|   `-- 04_host-taxonomy-virus-sharing/
|-- 04_viral_diversity_landscape_models/
|   |-- 01_adjusted-shannon-model/
|   `-- 02_spatial-autocorrelation/
|-- 05_environmental_association_and_gamm/
|   |-- 01_environmental-correlation-analysis/
|   |-- 02_environmental-variable-standardization/
|   |-- 03_environmental-correlation-vif-screening/
|   |-- 04_environmental-linear-nonlinear-form-assessment/
|   |-- 05_environmental-final-variable-selection/
|   `-- 06_environmental-driver-gamm-final-model/
|-- 06_rdrp_sequence_diversity_and_host_virus_network/
|   |-- 01_viral-order-butterfly-plot/
|   |-- 02_rdrp-amino-acid-identity-barplot/
|   `-- 03_rdrp-contig-host-virus-network-centrality/
|-- 07_coronavirus_surveillance_priority/
|   |-- 01_bat-host-surveillance-priority-assessment/
|   `-- 02_bat-crs-scatterplot/
|-- reproducibility/
|   |-- run_all_python_workflows.py
|   `-- run_all_r_workflows.py
|-- CITATION.cff
|-- CODE_OF_CONDUCT.md
|-- CONTRIBUTING.md
|-- DATA_AVAILABILITY.md
|-- DATA_LICENSE.md
|-- environment.yml
|-- LICENSE
|-- NOTICE.md
|-- README.md
|-- SECURITY.md
`-- requirements.txt
```

## Analysis Blocks

| Order | Repository block | Manuscript role |
|---|---|---|
| 01 | `01_sampling_and_host_representativeness/` | Sampling scale, in-house host coverage, Global Biodiversity Information Facility (GBIF) comparison and seasonal species accumulation/extrapolation checks. |
| 02 | `02_viral_detection_and_spectrum/` | Viral-family detection accumulation, in-study viral-spectrum heatmaps and previous-dataset host-genus spectrum comparison. |
| 03 | `03_virus_sharing_ecology/` | Geographic distance decay, terrain-associated sharing, host identity sharing and host-taxonomy sharing analyses. |
| 04 | `04_viral_diversity_landscape_models/` | Adjusted Shannon diversity comparison and spatial autocorrelation checks. |
| 05 | `05_environmental_association_and_gamm/` | Environmental correlation/Mantel summaries followed by an explicit seven-criterion functional-form assessment, independent parsimony selection and dynamically specified final environmental generalized additive mixed model (GAMM). |
| 06 | `06_rdrp_sequence_diversity_and_host_virus_network/` | Viral-order abundance, RdRp amino-acid identity summaries and RdRp host-virus network centrality. |
| 07 | `07_coronavirus_surveillance_priority/` | Coronavirus Ranking Score (CRS) calculation, host surveillance-priority assessment, selection of prioritized bat hosts and coronavirus surveillance priority visualization. |

## Workflow Convention

Each workflow folder uses the same local structure:

```text
workflow-name/
|-- scripts/      # executable R or Python scripts
|-- input/        # minimal public analysis-ready input tables
|-- output/       # created or refreshed locally when the workflow is run
|-- README.md     # workflow-specific notes
`-- requirements.txt or r_requirements.txt when needed
```

## Software Environment

Python scripts require Python 3.10 or later and the packages listed in the root and workflow-level `requirements.txt` files. Python requirements use compatible lower-bound version constraints to improve cross-platform environment creation. R scripts require R 4.3 or later and the packages listed in workflow-level R requirement files or `environment.yml`. The root environment files summarize the common dependencies used across the included modules: `pandas`, `numpy`, `matplotlib`, `openpyxl`, `reportlab`, `mgcv`, `emmeans`, `ggplot2`, `readr`, `dplyr`, `tidyr`, `tibble` and `iNEXT`.

Most plotting scripts try to use Times New Roman when it is available. On Linux or macOS, Matplotlib and R may fall back to another installed serif font unless Times-compatible fonts are installed. This can cause small visual differences in text metrics without changing the underlying numerical outputs.

## Measured Runtime

The following wall-clock runtimes were measured with `/usr/bin/time -p` in the tested environment above. Each script was run once against its included analysis-ready input data in a temporary copy of the repository. The table reports each measured `real` value rounded to the nearest second; runtime can vary with CPU load, storage performance and system configuration.

| Workflow | Script | Approximate measured wall-clock runtime |
|---|---|---:|
| `01_sampling_and_host_representativeness/01_pool-size-distribution` | `plot_pool_sample_count_distribution.py` | ~2 s |
| `01_sampling_and_host_representativeness/02_gbif-host-diversity-comparison` | `compare_inhouse_host_diversity_with_gbif.py` | ~1 s |
| `01_sampling_and_host_representativeness/03_seasonal-species-accumulation-extrapolation` | `01_plot_seasonal_species_accumulation_extrapolation.R` | ~25 s |
| `01_sampling_and_host_representativeness/03_seasonal-species-accumulation-extrapolation` | `02_calculate_terrain_inext_extrapolation_metrics.R` | ~22 s |
| `01_sampling_and_host_representativeness/03_seasonal-species-accumulation-extrapolation` | `03_plot_terrain_species_accumulation_resampling.py` | ~4 s |
| `02_viral_detection_and_spectrum/01_read-based-viral-family-detection-accumulation` | `plot_viral_family_detection_accumulation_curve.py` | ~2 s |
| `02_viral_detection_and_spectrum/02_province-host-viral-spectrum_in-this-study` | `plot_viral_spectrum.py` | ~3 s |
| `02_viral_detection_and_spectrum/03_host-genus-viral-spectrum_in-this-study` | `plot_viral_spectrum.py` | ~3 s |
| `02_viral_detection_and_spectrum/04_host-genus-viral-spectrum_previous-datasets` | `plot_viral_spectrum.py` | ~4 s |
| `03_virus_sharing_ecology/01_virus-sharing-distance-decay` | `analyze_virus_sharing_distance_decay.py` | ~32 s |
| `03_virus_sharing_ecology/02_terrain-virus-sharing` | `01_plot_terrain_shared_cluster_heatmap.py` | ~8 s |
| `03_virus_sharing_ecology/02_terrain-virus-sharing` | `02_test_terrain_association_with_virus_sharing.py` | ~40 s |
| `03_virus_sharing_ecology/03_host-identity-virus-sharing` | `analyze_host_identity_virus_sharing.py` | ~9 s |
| `03_virus_sharing_ecology/04_host-taxonomy-virus-sharing` | `analyze_host_taxonomy_virus_sharing.py` | ~4 s |
| `04_viral_diversity_landscape_models/01_adjusted-shannon-model` | `fit_adjusted_shannon_terrain_model.R` | ~2 s |
| `04_viral_diversity_landscape_models/02_spatial-autocorrelation` | `calculate_global_moran_spatial_autocorrelation.py` | ~2 s |
| `05_environmental_association_and_gamm/01_environmental-correlation-analysis` | `01_calculate_environmental_spearman_correlations.py` | ~63 s |
| `05_environmental_association_and_gamm/01_environmental-correlation-analysis` | `02_run_shannon_environment_mantel_tests.py` | ~9 s |
| `05_environmental_association_and_gamm/01_environmental-correlation-analysis` | `03_plot_environmental_correlation_mantel.py` | ~6 s |
| `05_environmental_association_and_gamm/02_environmental-variable-standardization` | `standardize_environmental_variables.py` | ~1 s |
| `05_environmental_association_and_gamm/03_environmental-correlation-vif-screening` | `screen_environmental_correlation_vif.py` | ~43 s |
| `05_environmental_association_and_gamm/04_environmental-linear-nonlinear-form-assessment` | `assess_environmental_model_forms.R` | ~41 s |
| `05_environmental_association_and_gamm/05_environmental-final-variable-selection` | `define_parsimonious_environmental_predictors.py` | ~3 s |
| `05_environmental_association_and_gamm/06_environmental-driver-gamm-final-model` | `fit_parsimonious_gamm.R` | ~8 s |
| `06_rdrp_sequence_diversity_and_host_virus_network/01_viral-order-butterfly-plot` | `plot_viral_order_butterfly.py` | ~3 s |
| `06_rdrp_sequence_diversity_and_host_virus_network/02_rdrp-amino-acid-identity-barplot` | `plot_rdrp_amino_acid_identity_barplot.py` | ~5 s |
| `06_rdrp_sequence_diversity_and_host_virus_network/03_rdrp-contig-host-virus-network-centrality` | `plot_rdrp_contig_host_virus_network_centrality.py` | ~4 s |
| `07_coronavirus_surveillance_priority/01_bat-host-surveillance-priority-assessment` | `01_calculate_entropy_weighted_crs.py` | ~4 s |
| `07_coronavirus_surveillance_priority/01_bat-host-surveillance-priority-assessment` | `02_select_prioritized_bat_hosts.py` | ~3 s |
| `07_coronavirus_surveillance_priority/02_bat-crs-scatterplot` | `plot_bat_crs_scatterplot.py` | ~8 s |

## Using Your Own Data

To use a workflow with user-supplied data, first make a separate working copy of that workflow so the released analysis-ready inputs and expected outputs remain available for comparison. Prepare replacement files using the corresponding files under `input/` as structural templates and follow the workflow-level README. Preserve the documented filenames, column names, data types, units, transformations, category labels and missing-value conventions because the scripts validate or directly reference these fields.

Place the replacement files in the copied workflow's `input/` directory, install the listed workflow dependencies, and run the corresponding script or scripts in the documented order. Generated figures, tables and diagnostics will be written to that workflow's `output/` directory. Compare the generated filenames and formats with the Expected Outputs section of the workflow README. Workflows that consume an upstream-derived table include a local analysis-ready copy of that table; users replacing it should reproduce the same structure and preserve the provenance relationship described in the README.

## Data and Coordinate Protection

Only small analysis-ready input tables required by the released scripts are included. Input provenance is described in the relevant workflow README files, especially when a downstream workflow uses a small upstream result table as its local input.

To protect sensitive cave roost locations, geographic coordinates are rounded to two decimal places. Detailed cave locations and site-use information are not publicly released. Access to precise locality information requires a reasonable request to the corresponding author, Dr. Zhiqiang Wu ([zqwu_lab@163.com](mailto:zqwu_lab@163.com)).

## License and Citation

Code and scripts in this repository are released under the MIT License. See `LICENSE`.

Repository documentation, README files, method notes and public analysis-ready input tables are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0), unless otherwise noted. See `DATA_LICENSE.md`.

Precise cave roost coordinates, sensitive locality information, controlled-access raw data and third-party materials are not relicensed by this repository. External datasets and software packages remain subject to their own licenses, citation requirements and terms of use.

If you use this repository, please cite the repository using `CITATION.cff` and cite the associated manuscript when the final citation becomes available.

## Output Policy

Workflow-level `output/` directories are created or refreshed when scripts are run. Expected output names, formats and roles are documented in each workflow README, so readers can verify generated products against the workflow contract even when pre-generated files are not included. Rerunning a workflow overwrites or refreshes the corresponding local products under that workflow's `output/` directory. Figure outputs are standardized to PDF; PNG, JPG, JPEG, TIFF and SVG files are not target outputs for the released workflows.

## Method-Reliability Notes

Workflow-specific README files document the required inputs, key parameters, transformations and analysis rationale. Analysis-ready inputs and generated outputs are organized at the workflow level to facilitate tracing of each released analysis.
