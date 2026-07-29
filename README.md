# Bat Anal-Swab Viral Discovery, Diversity, and Spatial Prioritization for Coronavirus Surveillance

Reproducible analysis code and public analysis-ready data accompanying the submitted manuscript.

## Study Overview

This repository contains reproducible workflows and analysis-ready data supporting the metatranscriptomic discovery of bat-associated viruses from 5,498 georeferenced bat anal swabs collected across 250 sites. The released workflows cover viral discovery and diversity, host-virus sharing, environmental associations, landscape-level analyses, and integrated spatial prioritization for bat coronavirus surveillance across China and mainland Southeast Asia.

## Manuscript

**Title:** Discovery of Novel Bat RNA Viruses and Integrated Spatial Prioritization Highlight Mountain Landscapes for Coronavirus Surveillance

**Authors:** Panpan Xu, Yelin Han, Kun Zhao, Wenliang Zhao, Shixuan Dong, Bo Liu, Xingyu Zhang, Yuyang Wang, Lamei Zhao, Xiujuan Yu, Qing Tang, Junpeng Zhang, Guangjian Zhu, Shuyi Zhang, Jian Yang, Qi Jin, Edward C. Holmes, and Zhiqiang Wu.

**Equal contribution:** Panpan Xu, Yelin Han, Kun Zhao, Wenliang Zhao, and Shixuan Dong.

**Correspondence:** Edward C. Holmes and Zhiqiang Wu.

**Status:** Manuscript submitted.

## Repository Organization

This repository organizes the main analysis scripts and minimal analysis-ready input tables used for method-level reproduction of the manuscript's code-supported analyses. Each numbered top-level folder represents one analysis block, and each subfolder is an individual workflow with its own `scripts/`, `input/` and `output/` directories.

Scripts write generated figures, tables and diagnostics to workflow-level `output/` directories. When a small generated table is required by a downstream workflow, it is copied into that workflow's `input/` directory and its provenance is described in the workflow README.

## Quick Start

Clone or download the repository, then create the software environment from the repository root.

```bash
conda env create -f environment.yml
conda activate manuscript_analysis_reproducibility
```

If conda is not available, Python-only workflows can be run with Python 3.10 or later after installing the root dependency summary:

```bash
python -m pip install -r requirements.txt
```

R workflows require R 4.3 or later and the packages listed in `environment.yml` and workflow-level R requirement files. With a standard R installation, missing R packages can be installed from CRAN, for example:

```r
install.packages(c("iNEXT", "mgcv", "emmeans", "ggplot2", "readr", "dplyr", "tidyr", "tibble"))
```

To rerun all released Python workflows from the repository root:

```bash
python reproducibility/run_all_python_workflows.py
```

To rerun all released R workflows from the repository root, make sure `Rscript` is on `PATH`, then run:

```bash
python reproducibility/run_all_r_workflows.py
```

The Python workflows include permutation-based analyses and may take several minutes on a typical desktop. Individual workflows can also be run directly from their own `scripts/` directory; each workflow README documents the expected inputs.

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
|-- 07_coronavirus_host_risk_prioritization/
|   |-- 01_bat-host-prioritization/
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
| 07 | `07_coronavirus_host_risk_prioritization/` | Coronavirus Risk Score (CRS) calculation, high-priority host selection and CRS visualization. |

## Workflow Convention

Each workflow folder uses the same local structure:

```text
workflow-name/
|-- scripts/      # executable R or Python scripts
|-- input/        # minimal public analysis-ready input tables
|-- output/       # created locally when workflows are run; not included in the repository
|-- README.md     # workflow-specific notes
`-- requirements.txt or r_requirements.txt when needed
```

## Software Environment

Python scripts require Python 3.10 or later and the packages listed in the root and workflow-level `requirements.txt` files. Python requirements use compatible lower-bound version constraints to improve cross-platform environment creation. R scripts require R 4.3 or later and the packages listed in workflow-level R requirement files or `environment.yml`. The root environment files summarize the common dependencies used across the included modules: `pandas`, `numpy`, `matplotlib`, `openpyxl`, `reportlab`, `mgcv`, `emmeans`, `ggplot2`, `readr`, `dplyr`, `tidyr`, `tibble` and `iNEXT`.

Most plotting scripts try to use Times New Roman when it is available. On Linux or macOS, Matplotlib and R may fall back to another installed serif font unless Times-compatible fonts are installed. This can cause small visual differences in text metrics without changing the underlying numerical outputs.

## Data and Coordinate Protection

Only small analysis-ready input tables required by the released scripts are included. Input provenance is described in the relevant workflow README files, especially when a downstream workflow uses a small upstream result table as its local input.

To protect sensitive cave roost locations, geographic coordinates are rounded to two decimal places. Detailed cave locations and site-use information are not publicly released. Access to precise locality information requires a reasonable request to the corresponding author, Dr. Zhiqiang Wu ([zqwu_lab@163.com](mailto:zqwu_lab@163.com)).

## License and Citation

Code and scripts in this repository are released under the MIT License. See `LICENSE`.

Repository documentation, README files, method notes and public analysis-ready input tables are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0), unless otherwise noted. See `DATA_LICENSE.md`.

Precise cave roost coordinates, sensitive locality information, controlled-access raw data and third-party materials are not relicensed by this repository. External datasets and software packages remain subject to their own licenses, citation requirements and terms of use.

If you use this repository, please cite the repository using `CITATION.cff` and cite the associated manuscript when the final citation becomes available.

## Output Policy

Workflow-level output/ directories are generated locally when scripts are run and are not included in the released repository. Expected output names and formats are documented in each workflow README.

## Method-Reliability Notes

Workflow-specific README files document the main input fields, parameters and analysis rationale. The repository uses standardized transformations and explicit input tables for downstream workflows so that each figure can be traced to the input fields required by its script. Count-like variables are transformed according to each workflow's analysis contract: zero-retaining pairwise matrices use log10(x + 1), while positive-only shared-cluster analyses use log10(x). Rank-based or permutation-based tests are used where distributional assumptions are not central to the method. Environmental predictors are standardized before screening and modelling so that effect scales, correlation checks and Variance Inflation Factor (VIF) calculations are comparable. Coordinate rounding is applied before public release to protect sensitive cave-roost localities while preserving province- and landscape-scale reproducibility.
