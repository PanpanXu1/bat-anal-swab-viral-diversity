from pathlib import Path

import numpy as np
import pandas as pd


ENV_VARS = [
    "Annual Mean Temperature (Bio1)",
    "Mean Diurnal Range (Mean of monthly max temp - min temp) (Bio2)",
    "Isothermality (BIO2/BIO7) (x100) (Bio3)",
    "Temperature Seasonality (standard deviation x100) (Bio4)",
    "Max Temperature of Warmest Month (Bio5)",
    "Min Temperature of Coldest Month (Bio6)",
    "Temperature Annual Range (BIO5-BIO6) (Bio7)",
    "Mean Temperature of Wettest Quarter (Bio8)",
    "Mean Temperature of Driest Quarter (Bio9)",
    "Mean Temperature of Warmest Quarter (Bio10)",
    "Mean Temperature of Coldest Quarter (Bio11)",
    "Annual Precipitation (Bio12)",
    "Precipitation of Wettest Month (Bio13)",
    "Precipitation of Driest Month (Bio14)",
    "Precipitation Seasonality (Coefficient of Variation) (Bio15)",
    "Precipitation of Wettest Quarter (Bio16)",
    "Precipitation of Driest Quarter (Bio17)",
    "Precipitation of Warmest Quarter (Bio18)",
    "Precipitation of Coldest Quarter (Bio19)",
    "Global Mammal Richness (GMR)",
    "Global Railway (GR)",
    "Global Linear Hydrography (GLH)",
    "Human Footprint (HFT)",
    "Fractional Vegetation Cover (FVC)",
    "Normalized Difference Vegetation Index (NDVI)",
    "China GDP Spatial Distribution (GDP)",
    "China High-Resolution Ecological Environment Quality (CHEQ)",
    "China Population Spatial Distribution (PSD)",
    "Digital Elevation Model (DEM)",
]
RETAINED_PREDICTORS = [
    "Mean Diurnal Range (Mean of monthly max temp - min temp) (Bio2)",
    "Temperature Seasonality (standard deviation x100) (Bio4)",
    "Max Temperature of Warmest Month (Bio5)",
    "Mean Temperature of Wettest Quarter (Bio8)",
    "Precipitation of Driest Quarter (Bio17)",
    "Precipitation of Warmest Quarter (Bio18)",
    "Human Footprint (HFT)",
    "Global Mammal Richness (GMR)",
    "China High-Resolution Ecological Environment Quality (CHEQ)",
    "Normalized Difference Vegetation Index (NDVI)",
    "Fractional Vegetation Cover (FVC)",
    "China Population Spatial Distribution (PSD)",
    "Global Railway (GR)",
    "Global Linear Hydrography (GLH)",
]
CORR_THRESHOLD = 0.8
VIF_THRESHOLD = 10.0
WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
MODULE_ROOT = WORKFLOW_ROOT.parent
MAXIMUM_SET_SIZES_TO_AUDIT = [16, 15, 14]


def calculate_vif(data, variables):
    x_matrix = data[variables].astype(float).to_numpy()
    n_rows, n_vars = x_matrix.shape
    records = []
    for idx, variable in enumerate(variables):
        y = x_matrix[:, idx]
        other_idx = [j for j in range(n_vars) if j != idx]
        if not other_idx:
            vif = 1.0
        else:
            x_other = x_matrix[:, other_idx]
            x_design = np.column_stack([np.ones(n_rows), x_other])
            coefficients, *_ = np.linalg.lstsq(x_design, y, rcond=None)
            y_hat = x_design @ coefficients
            ss_res = np.sum((y - y_hat) ** 2)
            ss_tot = np.sum((y - y.mean()) ** 2)
            r_squared = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan
            vif = np.inf if np.isclose(1 - r_squared, 0) else 1 / (1 - r_squared)
        records.append({"variable": variable, "VIF": vif})
    return pd.DataFrame(records).sort_values("VIF", ascending=False)


def max_vif(data, variables):
    return float(calculate_vif(data, list(variables))["VIF"].max())


def build_compatibility_graph(corr_matrix, variables):
    return {
        variable: {
            other
            for other in variables
            if other != variable and abs(corr_matrix.loc[variable, other]) < CORR_THRESHOLD
        }
        for variable in variables
    }


def enumerate_compatible_sets(compatibility, variables, target_size):
    ordered = sorted(variables, key=lambda variable: len(compatibility[variable]))
    compatible_sets = []

    def search(chosen, candidates):
        if len(chosen) + len(candidates) < target_size:
            return
        if len(chosen) == target_size:
            compatible_sets.append(tuple(chosen))
            return

        remaining = list(candidates)
        while remaining:
            variable = remaining.pop(0)
            next_candidates = set(remaining) & compatibility[variable]
            search(chosen + [variable], next_candidates)

    search([], set(ordered))
    return compatible_sets


def summarize_predictor_set(variables):
    return "; ".join(variables)


def audit_maximal_corr_vif_sets(data, corr_matrix, output_dir):
    compatibility = build_compatibility_graph(corr_matrix, ENV_VARS)
    retained_set = frozenset(RETAINED_PREDICTORS)
    audit_rows = []
    valid_14_rows = []
    maximum_vif_qualified_size = 0

    for target_size in MAXIMUM_SET_SIZES_TO_AUDIT:
        compatible_sets = enumerate_compatible_sets(compatibility, ENV_VARS, target_size)
        vif_summaries = []
        for predictor_set in compatible_sets:
            set_max_vif = max_vif(data, predictor_set)
            vif_summaries.append((predictor_set, set_max_vif))

        qualified_sets = [
            (predictor_set, set_max_vif)
            for predictor_set, set_max_vif in vif_summaries
            if set_max_vif < VIF_THRESHOLD
        ]
        if qualified_sets:
            maximum_vif_qualified_size = max(maximum_vif_qualified_size, target_size)

        qualified_sets_sorted = sorted(qualified_sets, key=lambda item: item[1])
        current_rank = ""
        current_max_vif = ""
        for rank, (predictor_set, set_max_vif) in enumerate(qualified_sets_sorted, start=1):
            if frozenset(predictor_set) == retained_set:
                current_rank = rank
                current_max_vif = set_max_vif
                break

        audit_rows.append(
            {
                "predictor_set_size": target_size,
                "pairwise_corr_compatible_sets_checked": len(compatible_sets),
                "vif_qualified_sets": len(qualified_sets),
                "minimum_max_vif": min((item[1] for item in vif_summaries), default=np.nan),
                "retained_set_is_this_size": len(RETAINED_PREDICTORS) == target_size,
                "retained_set_vif_qualified_at_this_size": current_rank != "",
                "retained_set_max_vif": current_max_vif,
                "retained_set_rank_by_lowest_max_vif": current_rank,
                "correlation_threshold_rule": f"absolute Spearman rho < {CORR_THRESHOLD}",
                "vif_threshold_rule": f"max VIF < {VIF_THRESHOLD}",
            }
        )

        if target_size == len(RETAINED_PREDICTORS):
            for rank, (predictor_set, set_max_vif) in enumerate(qualified_sets_sorted, start=1):
                row = {
                    "rank_by_lowest_max_vif": rank,
                    "predictor_count": target_size,
                    "max_vif": set_max_vif,
                    "is_retained_predictor_set": frozenset(predictor_set) == retained_set,
                    "predictors": summarize_predictor_set(predictor_set),
                }
                for idx, predictor in enumerate(predictor_set, start=1):
                    row[f"predictor_{idx}"] = predictor
                valid_14_rows.append(row)

    for row in audit_rows:
        row["maximum_vif_qualified_predictor_set_size"] = maximum_vif_qualified_size
        row["retained_set_is_maximum_size_vif_qualified_set"] = (
            len(RETAINED_PREDICTORS) == maximum_vif_qualified_size
        )

    pd.DataFrame(audit_rows).to_csv(
        output_dir / "maximal_corr_vif_predictor_set_audit.csv",
        index=False,
    )
    pd.DataFrame(valid_14_rows).to_csv(
        output_dir / "valid_14_predictor_sets_corr0.8_VIF10.csv",
        index=False,
    )


def main():
    input_path = WORKFLOW_ROOT / "input" / "standardized_environmental_variables_by_pool.csv"
    output_dir = WORKFLOW_ROOT / "output" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(input_path)

    missing = [var for var in ENV_VARS if var not in data.columns]
    if missing:
        raise ValueError(f"Missing environmental variables: {missing}")
    missing_retained = [var for var in RETAINED_PREDICTORS if var not in data.columns]
    if missing_retained:
        raise ValueError(f"Missing retained predictors: {missing_retained}")

    corr_matrix = data[ENV_VARS].corr(method="spearman")
    corr_matrix.to_csv(output_dir / "all_candidate_spearman_correlation_matrix.csv")

    records = []
    for i, factor1 in enumerate(ENV_VARS):
        for factor2 in ENV_VARS[i + 1:]:
            rho = corr_matrix.loc[factor1, factor2]
            records.append({"factor1": factor1, "factor2": factor2, "spearman_rho": rho, "abs_rho": abs(rho)})
    long_corr = pd.DataFrame(records).sort_values("abs_rho", ascending=False)
    high_corr_n = int((long_corr["abs_rho"] >= CORR_THRESHOLD).sum())
    if high_corr_n:
        print(f"High-correlation pairs with absolute Spearman rho >= {CORR_THRESHOLD}: {high_corr_n}")
    else:
        print(f"No variable pairs had absolute Spearman rho >= {CORR_THRESHOLD}.")

    selected = RETAINED_PREDICTORS
    final_vif = calculate_vif(data, selected)
    if (final_vif["VIF"] > VIF_THRESHOLD).any():
        raise ValueError("Retained predictor set contains VIF values above threshold.")
    audit_maximal_corr_vif_sets(data, corr_matrix, output_dir)
    final_vif.to_csv(output_dir / "selected_predictor_VIF_corr0.8_VIF10.csv", index=False)
    data[[col for col in data.columns if col not in ENV_VARS] + selected].to_csv(
        output_dir / "environmental_variables_selected_by_corr0.8_VIF10.csv", index=False
    )
    pd.DataFrame({"variable": selected}).to_csv(output_dir / "selected_environmental_predictors.csv", index=False)


if __name__ == "__main__":
    main()
