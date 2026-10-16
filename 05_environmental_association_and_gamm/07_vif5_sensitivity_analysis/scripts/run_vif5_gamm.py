"""Fit the ten-predictor VIF < 5 GAMM and export its term table."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MODEL_INPUT = "environmental_variables_VIF5.csv"
PREDICTOR_INPUT = "parsimonious_environmental_predictors.csv"
FORM_INPUT = "environmental_predictor_model_form_assessment.csv"
COLUMNS = ("term", "role", "model_form", "term_type", "estimate", "std_error", "edf", "statistic", "p_value")


def calculate_vif(data: pd.DataFrame, predictors: list[str]) -> pd.DataFrame:
    matrix = data[predictors].astype(float).to_numpy()
    records = []
    for index in range(matrix.shape[1]):
        response = matrix[:, index]
        other = matrix[:, [j for j in range(matrix.shape[1]) if j != index]]
        design = np.column_stack([np.ones(len(response)), other])
        coefficients, *_ = np.linalg.lstsq(design, response, rcond=None)
        residual = response - design @ coefficients
        r_squared = 1 - np.sum(residual ** 2) / np.sum((response - response.mean()) ** 2)
        vif = np.inf if np.isclose(1 - r_squared, 0) else 1 / (1 - r_squared)
        records.append({"variable": predictors[index], "VIF": vif})
    return pd.DataFrame(records).sort_values("VIF", ascending=False)


def number(value: object) -> str:
    return "" if pd.isna(value) else f"{float(value):.2f}"


def p_value(value: object) -> str:
    p = float(value)
    if p < 0.001:
        return "P < 0.001"
    if p < 0.01 or 0.045 <= p < 0.055:
        return f"P = {p:.3f}"
    return f"{p:.2f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default=shutil.which("Rscript"), help="Rscript executable")
    args = parser.parse_args()
    if not args.rscript:
        parser.error("Rscript was not found; pass --rscript PATH")
    input_dir = ROOT / "input"
    output_tables = ROOT / "output" / "tables"
    data = pd.read_csv(input_dir / MODEL_INPUT)
    predictors = pd.read_csv(input_dir / PREDICTOR_INPUT)["variable"].tolist()
    assessments = pd.read_csv(input_dir / FORM_INPUT)["variable"].tolist()
    if len(predictors) != 10 or set(data.columns) & set(assessments) != set(predictors):
        raise ValueError("The modelling input must contain exactly the ten VIF < 5 predictors")
    final_vif = calculate_vif(data, predictors)
    max_vif = float(final_vif["VIF"].max())
    if max_vif >= 5:
        raise ValueError(f"The modelling input exceeds VIF < 5: {max_vif:.2f}")

    output_tables.mkdir(parents=True, exist_ok=True)
    reference = ROOT.parent / "06_environmental-driver-gamm-final-model" / "input"
    reference_data = pd.read_csv(reference / "environmental_variables_selected_by_corr0.8_VIF10.csv")
    primary_predictors = pd.read_csv(reference / PREDICTOR_INPUT)["variable"].tolist()
    retained = primary_predictors.copy()
    while True:
        stage_vif = calculate_vif(reference_data, retained)
        highest = stage_vif.iloc[0]
        if highest["VIF"] < 5:
            break
        retained.remove(highest["variable"])
    if retained != predictors or not np.allclose(reference_data[predictors], data[predictors], rtol=0, atol=1e-12):
        raise RuntimeError("Saved ten-predictor input does not match iterative VIF < 5 selection from the final 13")

    with tempfile.TemporaryDirectory(prefix="vif5_gamm_") as temporary:
        case = Path(temporary) / "gamm"
        (case / "input").mkdir(parents=True)
        (case / "scripts").mkdir()
        for filename in (MODEL_INPUT, PREDICTOR_INPUT, FORM_INPUT):
            shutil.copy2(input_dir / filename, case / "input" / filename)
        source_r = ROOT / "scripts" / "fit_parsimonious_gamm.R"
        shutil.copy2(source_r, case / "scripts" / source_r.name)
        run = subprocess.run([args.rscript, str(case / "scripts" / source_r.name)], capture_output=True, text=True, encoding="utf-8", errors="replace")
        if run.returncode:
            raise RuntimeError(f"GAMM fit failed:\n{run.stdout}\n{run.stderr}")
        observations = re.search(r"complete model rows:\s*(\d+)", run.stdout)
        deviance = re.search(r"Deviance explained \(%\):\s*([0-9.]+)", run.stdout)
        if not observations or not deviance:
            raise RuntimeError(f"GAMM fit metrics missing:\n{run.stdout}")
        shutil.copy2(case / "output" / "tables" / "model_fit_summary.csv", output_tables / "model_fit_summary.csv")
        shutil.copy2(case / "output" / "tables" / "environmental_predictor_vif.csv", output_tables / "environmental_predictor_vif.csv")

    fitted = pd.read_csv(output_tables / "model_fit_summary.csv")
    if len(fitted) != 26 or fitted.loc[fitted["role"] == "environmental predictor", "term"].tolist() != predictors:
        raise RuntimeError("Fitted terms do not match the ten-predictor GAMM")
    rows = []
    for item in fitted.to_dict(orient="records"):
        rows.append({
            "term": item["term"], "role": item["role"],
            "model_form": "linear fixed effect" if item["model_form"] == "environmental linear fixed effect" else item["model_form"],
            "term_type": item["term_type"], "estimate": number(item["estimate"]),
            "std_error": number(item["std_error"]), "edf": number(item["edf"]),
            "statistic": number(item["statistic"]), "p_value": p_value(item["p_value"]),
        })
    with (output_tables / "gamm_vif5_term_estimates.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Ten predictors; maximum VIF {max_vif:.2f}; {observations.group(1)} pools; deviance explained {float(deviance.group(1)):.2f}%")


if __name__ == "__main__":
    main()
