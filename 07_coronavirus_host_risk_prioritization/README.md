# 07 Coronavirus Host Risk Prioritization

This block organizes the coronavirus host-risk scoring, representative high-priority host selection and CRS scatterplot workflows.

## Workflow Folders

| Order | Workflow | Purpose |
|---|---|---|
| 1 | `01_bat-host-prioritization/` | Calculates entropy-weighted CRS values and selects representative high-priority bat hosts using documented data-driven thresholds. |
| 2 | `02_bat-crs-scatterplot/` | Generates the CRS scatterplot from the selected-host plotting input produced by the prioritization workflow. |

Run the CRS calculation script before the host-selection script in `01_bat-host-prioritization/`. The host-selection script also updates the reduced input table used by `02_bat-crs-scatterplot/`.

## Method Rationale

The scoring workflow retains the full host-level CRS table, threshold summary and per-host audit trail. Selection is based on CRS percentile ranking and predefined evidence-structure filters. The scatterplot workflow then visualizes only the selected representative high-priority hosts, keeping calculation, selection audit and plotting steps traceable but separated.
