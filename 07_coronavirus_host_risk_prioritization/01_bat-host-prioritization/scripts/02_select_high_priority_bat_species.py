from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "output" / "tables"
SCORES_PATH = RESULTS_DIR / "bat_risk_score_ewm.csv"
THRESHOLD_PATH = RESULTS_DIR / "selection_threshold_summary.csv"
AUDIT_PATH = RESULTS_DIR / "selection_audit_trail.csv"
FINAL_SELECTION_PATH = RESULTS_DIR / "high_priority_bat_species_for_sdm.csv"
SCATTER_INPUT_PATH = (
    PROJECT_ROOT.parent
    / "02_bat-crs-scatterplot"
    / "input"
    / "bat_crs_scatterplot_input.csv"
)
EXPECTED_FINAL_SELECTION = 12
OUTPUT_DECIMAL_PLACES = 6
MIN_ZOONOTIC_SEQUENCE_NUMBER = 2

REQUIRED_COLUMNS = [
    "Host",
    "Zoonotic Sequence Number",
    "EWM_rank",
    "EWM_CRS",
    "VD_norm",
    "ZSC_norm",
    "ZSN_corrected_norm",
    "VS_ln_norm",
]

NUMERIC_COLUMNS = REQUIRED_COLUMNS[1:]


def read_csv_with_fallback(path: Path) -> pd.DataFrame:
    for encoding in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path)


def warn(message: str) -> None:
    print(f"WARNING: {message}")


def write_rounded_csv(data: pd.DataFrame, path: Path) -> None:
    output = data.copy()
    for column in output.select_dtypes(include=[np.number]).columns:
        if pd.api.types.is_integer_dtype(output[column]):
            continue
        output[column] = output[column].round(OUTPUT_DECIMAL_PLACES)
    output.to_csv(path, index=False, encoding="utf-8-sig")


def assign_selection(row: pd.Series) -> pd.Series:
    if not row["passes_CRS_P90_floor"]:
        return pd.Series(
            {
                "final_selected": False,
                "selection_path": "Not selected",
                "exclusion_reason": "Excluded: EWM_CRS below the CRS P90 high-priority floor.",
            }
        )
    if not row["passes_ZSN_floor"]:
        return pd.Series(
            {
                "final_selected": False,
                "selection_path": "Not selected",
                "exclusion_reason": (
                    "Excluded: zoonotic sequence evidence below the minimum repeated-sequence "
                    f"floor of {MIN_ZOONOTIC_SEQUENCE_NUMBER}."
                ),
            }
        )
    if row["passes_path_A"]:
        return pd.Series(
            {
                "final_selected": True,
                "selection_path": "Path A",
                "exclusion_reason": "Selected via Path A.",
            }
        )
    if row["passes_path_B"]:
        return pd.Series(
            {
                "final_selected": True,
                "selection_path": "Path B",
                "exclusion_reason": "Selected via Path B.",
            }
        )
    if np.isclose(row["ZSC_norm"], 1.00):
        return pd.Series(
            {
                "final_selected": False,
                "selection_path": "Not selected",
                "exclusion_reason": "Excluded: ZSC_norm = 1.00 but VD_norm below P90.",
            }
        )
    if np.isclose(row["ZSC_norm"], 0.50):
        return pd.Series(
            {
                "final_selected": False,
                "selection_path": "Not selected",
                "exclusion_reason": "Excluded: ZSC_norm = 0.50 but VD_norm below P95.",
            }
        )
    return pd.Series(
        {
            "final_selected": False,
            "selection_path": "Not selected",
            "exclusion_reason": "Excluded: did not meet either Path A or Path B.",
        }
    )


def make_threshold_summary(scores: pd.DataFrame, crs_p90: float, vd_p90: float, vd_p95: float) -> pd.DataFrame:
    n_hosts = len(scores)
    return pd.DataFrame(
        [
            {
                "threshold_name": "CRS_P90",
                "source_column": "EWM_CRS",
                "quantile": 0.90,
                "threshold_value": crs_p90,
                "n_hosts": n_hosts,
                "calculation_basis": "Full host-level EWM_CRS distribution after taxonomic harmonization.",
                "rule": "EWM_CRS >= CRS_P90",
            },
            {
                "threshold_name": "ZSN_min_sequence_evidence",
                "source_column": "Zoonotic Sequence Number",
                "quantile": None,
                "threshold_value": MIN_ZOONOTIC_SEQUENCE_NUMBER,
                "n_hosts": n_hosts,
                "calculation_basis": "Minimum repeated sequence-evidence floor applied to the raw zoonotic sequence count.",
                "rule": f"Zoonotic Sequence Number >= {MIN_ZOONOTIC_SEQUENCE_NUMBER}",
            },
            {
                "threshold_name": "VD_P90",
                "source_column": "VD_norm",
                "quantile": 0.90,
                "threshold_value": vd_p90,
                "n_hosts": n_hosts,
                "calculation_basis": "Full host-level normalized viral-diversity distribution.",
                "rule": "Used in Path A with ZSC_norm = 1.00 and the minimum repeated zoonotic-sequence evidence floor.",
            },
            {
                "threshold_name": "VD_P95",
                "source_column": "VD_norm",
                "quantile": 0.95,
                "threshold_value": vd_p95,
                "n_hosts": n_hosts,
                "calculation_basis": "Full host-level normalized viral-diversity distribution.",
                "rule": "Used in Path B with ZSC_norm = 0.50 and the minimum repeated zoonotic-sequence evidence floor.",
            },
        ]
    )


def write_scatterplot_input(final_selection: pd.DataFrame) -> None:
    scatter_input = final_selection[
        [
            "Host",
            "VD_norm",
            "ZSC_norm",
            "ZSN_corrected_norm",
            "VS_ln_norm",
            "EWM_CRS",
        ]
    ].rename(
        columns={
            "VD_norm": "Viral Diversity_norm",
            "ZSC_norm": "Zoonotic Species Count_norm",
            "VS_ln_norm": "VS_log_norm",
            "EWM_CRS": "CRS",
        }
    )
    SCATTER_INPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    scatter_input.to_csv(SCATTER_INPUT_PATH, index=False, encoding="utf-8-sig")


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    scores = read_csv_with_fallback(SCORES_PATH)

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in scores.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    for column in NUMERIC_COLUMNS:
        scores[column] = pd.to_numeric(scores[column], errors="raise")

    scores = scores.sort_values(["EWM_rank", "Host"], ascending=[True, True]).copy()
    crs_p90 = scores["EWM_CRS"].quantile(0.90)
    vd_p90 = scores["VD_norm"].quantile(0.90)
    vd_p95 = scores["VD_norm"].quantile(0.95)

    threshold_summary = make_threshold_summary(scores, crs_p90, vd_p90, vd_p95)
    write_rounded_csv(threshold_summary, THRESHOLD_PATH)

    audit = scores.copy()
    audit["CRS_P90"] = crs_p90
    audit["VD_P90"] = vd_p90
    audit["VD_P95"] = vd_p95
    audit["passes_CRS_P90_floor"] = audit["EWM_CRS"] >= crs_p90
    audit["passes_ZSN_floor"] = audit["Zoonotic Sequence Number"] >= MIN_ZOONOTIC_SEQUENCE_NUMBER
    audit["passes_path_A"] = (
        audit["passes_ZSN_floor"]
        & np.isclose(audit["ZSC_norm"], 1.00)
        & (audit["VD_norm"] >= vd_p90)
    )
    audit["passes_path_B"] = (
        audit["passes_ZSN_floor"]
        & np.isclose(audit["ZSC_norm"], 0.50)
        & (audit["VD_norm"] >= vd_p95)
    )

    selection_columns = audit.apply(assign_selection, axis=1)
    audit = pd.concat([audit, selection_columns], axis=1)

    audit_columns = [
        "Host",
        "Zoonotic Sequence Number",
        "EWM_rank",
        "EWM_CRS",
        "VD_norm",
        "ZSC_norm",
        "ZSN_corrected_norm",
        "VS_ln_norm",
        "CRS_P90",
        "VD_P90",
        "VD_P95",
        "passes_CRS_P90_floor",
        "passes_ZSN_floor",
        "passes_path_A",
        "passes_path_B",
        "final_selected",
        "selection_path",
        "exclusion_reason",
    ]
    audit = audit[audit_columns].sort_values(["EWM_rank", "Host"]).copy()
    write_rounded_csv(audit, AUDIT_PATH)

    final_selection = audit[audit["final_selected"]].copy()
    final_selection["selection_status"] = "Selected_for_downstream_spatial_handoff"
    final_selection["selection_reason"] = final_selection["exclusion_reason"]
    final_columns = [
        "Host",
        "EWM_rank",
        "EWM_CRS",
        "VD_norm",
        "ZSC_norm",
        "ZSN_corrected_norm",
        "VS_ln_norm",
        "selection_path",
        "selection_status",
        "selection_reason",
    ]
    final_selection = final_selection[final_columns].sort_values(["EWM_rank", "Host"])
    write_rounded_csv(final_selection, FINAL_SELECTION_PATH)
    write_scatterplot_input(final_selection)

    if len(final_selection) != EXPECTED_FINAL_SELECTION:
        warn(
            f"Expected {EXPECTED_FINAL_SELECTION} final selected hosts, "
            f"but found {len(final_selection)}."
        )

    print(f"Loaded species records: {len(scores)}")
    print(f"CRS P90: {crs_p90:.{OUTPUT_DECIMAL_PLACES}f}")
    print(f"VD_norm P90: {vd_p90:.{OUTPUT_DECIMAL_PLACES}f}")
    print(f"VD_norm P95: {vd_p95:.{OUTPUT_DECIMAL_PLACES}f}")
    print(f"Selected via Path A: {int((final_selection['selection_path'] == 'Path A').sum())}")
    print(f"Selected via Path B: {int((final_selection['selection_path'] == 'Path B').sum())}")
    print(f"Selected representative high-priority hosts: {len(final_selection)}")
    print(f"Output saved: {FINAL_SELECTION_PATH.relative_to(PROJECT_ROOT)}")
    print(f"Scatterplot input saved: {SCATTER_INPUT_PATH.relative_to(PROJECT_ROOT.parent)}")


if __name__ == "__main__":
    main()
