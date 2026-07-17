import math
from pathlib import Path

import numpy as np
import pandas as pd


WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
MODULE_ROOT = WORKFLOW_ROOT.parent
P_VALUE_THRESHOLD = 0.05
MODEL_FORM_SOURCE = "workflow_04_seven_criterion_assessment"
CRITERION_COLUMNS = [
    "criterion_1_aic",
    "criterion_2_edf",
    "criterion_3_p",
    "criterion_4_k_check",
    "criterion_5_edf_ceiling",
    "criterion_6_curve_stability",
    "criterion_7_boundary_robustness",
]


def parse_boolean(value, column, variable):
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().upper() in {"TRUE", "FALSE"}:
        return value.strip().upper() == "TRUE"
    raise ValueError(
        f"Invalid boolean in {column!r} for {variable!r}: {value!r}; "
        "expected TRUE or FALSE."
    )


def parse_selected_k(value, selected_form, variable):
    if pd.isna(value) or (isinstance(value, str) and not value.strip()):
        if selected_form == "nonlinear":
            raise ValueError(f"Nonlinear predictor {variable!r} requires selected_k.")
        return pd.NA
    if selected_form == "linear":
        raise ValueError(
            f"Linear predictor {variable!r} must have an empty selected_k; got {value!r}."
        )
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid selected_k for {variable!r}: {value!r}.") from exc
    if not math.isfinite(numeric) or numeric != int(numeric) or numeric < 4:
        raise ValueError(
            f"Invalid selected_k for {variable!r}: {value!r}; expected a finite integer >= 4."
        )
    return int(numeric)


def validate_model_form_assessment(assessment, environmental_terms):
    required = {"variable", "selected_form", "selected_k", *CRITERION_COLUMNS}
    missing = required - set(assessment.columns)
    if missing:
        raise ValueError(f"Model-form assessment is missing required columns: {sorted(missing)}")
    duplicates = assessment.loc[assessment["variable"].duplicated(keep=False), "variable"].tolist()
    if duplicates:
        raise ValueError(f"Model-form assessment contains duplicate variables: {sorted(set(duplicates))}")

    selected = assessment[assessment["variable"].isin(environmental_terms)].copy()
    counts = selected["variable"].value_counts()
    invalid = [variable for variable in environmental_terms if counts.get(variable, 0) != 1]
    if invalid:
        raise ValueError(
            "Model-form assessment must contain exactly one row for every screened predictor; "
            f"invalid variables: {invalid}"
        )
    selected = selected.set_index("variable").loc[environmental_terms].reset_index()
    normalized_k = []
    for _, row in selected.iterrows():
        variable = row["variable"]
        form = row["selected_form"]
        if form not in {"linear", "nonlinear"}:
            raise ValueError(f"Invalid selected_form for {variable!r}: {form!r}.")
        criteria = [parse_boolean(row[column], column, variable) for column in CRITERION_COLUMNS]
        if form == "nonlinear" and not all(criteria):
            failed = [column for column, passed in zip(CRITERION_COLUMNS, criteria) if not passed]
            raise ValueError(
                f"Nonlinear predictor {variable!r} does not pass all seven criteria; failed: {failed}"
            )
        if form == "linear" and all(criteria):
            raise ValueError(
                f"Linear predictor {variable!r} passes all seven criteria and must be nonlinear."
            )
        normalized_k.append(parse_selected_k(row["selected_k"], form, variable))
    selected["selected_k"] = pd.array(normalized_k, dtype="Int64")
    return selected[["variable", "selected_form", "selected_k"]]


def validate_contribution_rows(contribution, environmental_terms):
    selected = contribution[contribution["variable"].isin(environmental_terms)].copy()
    counts = selected["variable"].value_counts()
    missing = [variable for variable in environmental_terms if counts.get(variable, 0) == 0]
    if missing:
        raise ValueError(
            "Preliminary contribution table is missing screened predictors: "
            f"{missing}"
        )
    duplicates = [variable for variable in environmental_terms if counts.get(variable, 0) > 1]
    if duplicates:
        raise ValueError(
            "Preliminary contribution table contains duplicate screened predictors: "
            f"{duplicates}"
        )
    return selected.set_index("variable").loc[environmental_terms].reset_index()


def main():
    input_path = (
        WORKFLOW_ROOT
        / "input"
        / "environmental_variables_selected_by_corr0.8_VIF10.csv"
    )
    contribution_path = WORKFLOW_ROOT / "input" / "preliminary_environmental_contribution.csv"
    assessment_path = WORKFLOW_ROOT / "input" / "environmental_predictor_model_form_assessment.csv"
    output_dir = WORKFLOW_ROOT / "output" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(input_path)
    contribution = pd.read_csv(contribution_path)
    assessment = pd.read_csv(assessment_path)
    required_contribution_columns = {
        "variable",
        "delta_deviance_explained",
        "delta_AIC",
        "full_model_term_p_value",
    }
    missing_contribution_columns = required_contribution_columns - set(contribution.columns)
    if missing_contribution_columns:
        raise ValueError(
            "Preliminary contribution table is missing required columns: "
            f"{sorted(missing_contribution_columns)}"
        )
    environmental_terms = [
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
    selected = [var for var in environmental_terms if var in data.columns]
    if selected != environmental_terms:
        missing_terms = [var for var in environmental_terms if var not in data.columns]
        raise ValueError(f"Screened modelling table is missing predictors: {missing_terms}")
    selected_forms = validate_model_form_assessment(assessment, environmental_terms)
    selected_contribution = validate_contribution_rows(contribution, environmental_terms)
    for column in ("delta_deviance_explained", "delta_AIC", "full_model_term_p_value"):
        selected_contribution[column] = pd.to_numeric(
            selected_contribution[column], errors="raise"
        )
        non_finite = ~np.isfinite(selected_contribution[column].to_numpy(dtype=float))
        if non_finite.any():
            variables = selected_contribution.loc[non_finite, "variable"].tolist()
            raise ValueError(
                f"Preliminary contribution column {column!r} contains non-finite "
                f"values for variables: {variables}"
            )
    invalid_p = ~selected_contribution["full_model_term_p_value"].between(0, 1)
    if invalid_p.any():
        variables = selected_contribution.loc[invalid_p, "variable"].tolist()
        raise ValueError(
            "Preliminary contribution column 'full_model_term_p_value' must be "
            f"between 0 and 1 for variables: {variables}"
        )

    parsimony_removed = selected_contribution[
        (selected_contribution["delta_deviance_explained"] < 0)
        & (selected_contribution["full_model_term_p_value"] >= P_VALUE_THRESHOLD)
    ].copy()
    removed_terms = set(parsimony_removed["variable"])
    final_predictors = [var for var in selected if var not in removed_terms]
    final_table = selected_forms[selected_forms["variable"].isin(final_predictors)].copy()
    final_table["model_form_source"] = MODEL_FORM_SOURCE

    removed_records = []
    for _, row in parsimony_removed.iterrows():
        removed_records.append(
            {
                "variable": row["variable"],
                "reason": (
                    "Removed by the predefined parsimony rule: negative conditional "
                    "deviance contribution and no statistical support in the preliminary "
                    "adjusted model."
                ),
                "preliminary_delta_deviance_explained": row["delta_deviance_explained"],
                "preliminary_delta_AIC": row["delta_AIC"],
                "preliminary_p_value": row["full_model_term_p_value"],
                "p_value_threshold": P_VALUE_THRESHOLD,
            }
        )

    final_table.to_csv(
        output_dir / "parsimonious_environmental_predictors.csv",
        index=False,
    )
    pd.DataFrame(removed_records).to_csv(output_dir / "parsimony_removed_terms.csv", index=False)


if __name__ == "__main__":
    main()
