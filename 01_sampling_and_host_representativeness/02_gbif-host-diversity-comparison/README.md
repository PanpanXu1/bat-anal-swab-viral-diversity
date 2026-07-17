# 02 Global Biodiversity Information Facility (GBIF) Host Diversity Comparison

## Purpose

Compare host taxonomic coverage in the in-house sampling data with curated Global Biodiversity Information Facility (GBIF) occurrence records from the manuscript-relevant provinces. This workflow documents how the host taxa represented by field sampling relate to independent occurrence records.

## Input Files

1. `input/inhouse_host_species_by_province.csv` - in-house host species, genus and family records by standardized province.
2. `input/gbif_curated_occurrence_records.csv` - curated GBIF records with taxonomic fields, rounded coordinates, province assignment and basic record provenance. The `Province_standard` field is already curated before this repository step. Raw GBIF downloading, nine-province filtering and coordinate-to-province boundary assignment were completed before this repository step.

Coordinates are rounded to two decimal places in accordance with the repository-wide locality-protection policy.

## Script Workflow

1. `scripts/compare_inhouse_host_diversity_with_gbif.py`
   - Reads the in-house host table and curated GBIF occurrence table.
   - Standardizes taxonomic comparison fields at species and genus levels.
   - Counts shared and source-specific taxa, summarizes GBIF record support, and writes a multi-sheet workbook.
   - Generates four PDF figures from the same summary objects so that workbook and figures are traceable to the same input records.

## Expected Outputs

### Figures

1. `output/figures/inhouse_data_vs_GBIF_species_genus_richness.pdf` - species- and genus-level richness comparison between the in-house data and curated GBIF records.
2. `output/figures/inhouse_data_only_species_absent_from_GBIF.pdf` - species present in the in-house data but not represented in the curated GBIF subset, summarized by genus.
3. `output/figures/species_coverage_by_GBIF_record_threshold.pdf` - sensitivity curve showing how in-house species coverage changes as the minimum GBIF record-support threshold increases.
4. `output/figures/genus_coverage_by_GBIF_record_threshold.pdf` - genus-level version of the GBIF record-support threshold check.

### Tables

1. `output/tables/Inhouse_data_vs_GBIF.xlsx`
   - `Summary`: overall counts for in-house, GBIF, shared and source-specific taxa.
   - `Taxon_comparison`: species-level comparison table with taxon status, GBIF record counts, province counts and support classes.
   - `Inhouse_data_only_species`: subset of species present in the in-house data but absent from curated GBIF records.
   - `GBIF_species_frequency_support`: GBIF record-support details by species, including province coverage.
   - `Coverage_by_GBIF_frequency`: species-level coverage summarized across GBIF record-count thresholds.
   - `GBIF_genus_frequency_support`: GBIF record-support details by genus.
   - `Genus_coverage_by_frequency`: genus-level coverage summarized across GBIF record-count thresholds.
   - `GBIF_assignment_QC` and `GBIF_province_QC`: provenance checks for province assignment and GBIF metadata fields.

## Provenance Notes

The GBIF table is a curated analysis input. It keeps `decimalLongitude`, `decimalLatitude`, `Province_assignment_method`, `basisOfRecord` and `coordinateUncertaintyInMeters` so readers can identify how the curated province field was supported without exposing precise locality information.

## Method Rationale

GBIF is used as the independent host-occurrence reference because it provides record-level, taxonomically curated occurrence data that can be filtered to the same province scope as the in-house sampling data. This is more appropriate for this comparison than broad IUCN range polygons, which are estimated area layers and can be too coarse for evaluating host representation within the nine manuscript-relevant provinces.

The curated GBIF table was restricted to the nine provinces relevant to this study before this repository step. For records with a usable province label, the standardized province field was taken from the supplied locality information. For records lacking a usable province label but containing coordinates, rounded coordinates were used to assign the record to a province. The provenance fields in `input/gbif_curated_occurrence_records.csv` retain this distinction through `Province_assignment_method`, while `decimalLongitude`, `decimalLatitude`, `basisOfRecord` and `coordinateUncertaintyInMeters` provide record-level QC context.

GBIF contains many sparse and occasionally rare occurrence records. A taxon represented by a single record may be less stable as external evidence than a taxon documented repeatedly across records or provinces. Therefore, coverage is summarized across GBIF record-count thresholds rather than relying only on an all-records comparison. These threshold plots show how in-house coverage changes when increasingly well-supported GBIF taxa are considered, while leaving the underlying curated input records unchanged.
