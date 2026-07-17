# Data and Code Availability

This repository provides the analysis scripts and small analysis-ready input tables required for the included analysis modules. The repository is intended for method-level reproduction of the code-supported analyses, including input handling, parameter choices, statistical procedures and figure/table generation.

## Included Analysis Inputs

Input files are stored in the `input/` directory of each workflow. These files are curated, analysis-ready tables containing the fields required by the corresponding scripts. When a workflow uses a small table produced by an upstream step, that table is copied into the downstream `input/` directory and its provenance is described in the workflow README. This keeps each workflow locally runnable while preserving the upstream relationship.

## Generated Outputs

Scripts write generated figures, model summaries, diagnostic tables and intermediate products to workflow-level `output/` directories. These files are generated products rather than primary input data. If a generated table is required by a later workflow, it is provided in that later workflow's `input/` directory with a provenance note.

## Public Data and Documentation License

Public analysis-ready input tables, README files and method documentation included in this repository are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0), unless otherwise noted. See `DATA_LICENSE.md` for scope and attribution details.

This public license does not apply to precise cave roost coordinates, sensitive locality information, controlled-access raw data, unreleased primary data or third-party materials that are not included in this repository.

## Sensitive Locality Data

To protect sensitive cave roost locations, geographic coordinates are rounded to two decimal places. Detailed cave locations and site-use information are not publicly released. Access to precise locality information requires a reasonable request to the corresponding author, Dr. Zhiqiang Wu ([zqwu_lab@163.com](mailto:zqwu_lab@163.com)).

## Software Environment

Python and R dependencies are summarized in the root `requirements.txt`, workflow-level requirement files and `environment.yml`. Workflow-specific README files describe the expected input files, scripts, generated outputs and method notes.

## Contact

Questions about data access, sensitive locality information or repository use should be directed to the corresponding author, Dr. Zhiqiang Wu ([zqwu_lab@163.com](mailto:zqwu_lab@163.com)).
