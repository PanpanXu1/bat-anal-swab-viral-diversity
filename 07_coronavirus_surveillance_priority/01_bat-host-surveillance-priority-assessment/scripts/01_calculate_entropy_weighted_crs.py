from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = PROJECT_ROOT / "input" / "bat_coronavirus_host_surveillance_priority_indicators.csv"
RESULTS_DIR = PROJECT_ROOT / "output" / "tables"
WEIGHTS_PATH = RESULTS_DIR / "entropy_weights.csv"
SCORES_PATH = RESULTS_DIR / "bat_crs_rankings_ewm.csv"
EXPECTED_HOST_RECORDS = 229
OUTPUT_DECIMAL_PLACES = 6

RAW_COLUMNS = [
    "Host",
    "Viral Diversity",
    "Zoonotic Species Count",
    "Zoonotic Sequence Number",
    "Viral Sequence Count",
]

NUMERIC_COLUMNS = [
    "Viral Diversity",
    "Zoonotic Species Count",
    "Zoonotic Sequence Number",
    "Viral Sequence Count",
]

NORMALIZED_COLUMNS = [
    "VD_norm",
    "ZSC_norm",
    "ZSN_corrected_norm",
    "VS_ln_norm",
]

INDICATOR_LABELS = [
    "VD_norm",
    "ZSC_norm",
    "ZSN_corrected_norm",
    "VS_ln_norm",
]


def read_csv_with_fallback(path: Path) -> pd.DataFrame:
    for encoding in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding, low_memory=False)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, low_memory=False)


def warn(message: str) -> None:
    print(f"WARNING: {message}")


def min_max_normalize(series: pd.Series) -> pd.Series:
    min_value = series.min()
    max_value = series.max()
    if max_value == min_value:
        warn(f"{series.name} has a constant value; normalized values were set to 0.")
        return pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    return (series - min_value) / (max_value - min_value)


def write_rounded_csv(data: pd.DataFrame, path: Path) -> None:
    output = data.copy()
    for column in output.select_dtypes(include=[np.number]).columns:
        if pd.api.types.is_integer_dtype(output[column]):
            continue
        output[column] = output[column].round(OUTPUT_DECIMAL_PLACES)
    output.to_csv(path, index=False, encoding="utf-8-sig")


def prepare_input() -> pd.DataFrame:
    data = read_csv_with_fallback(DATA_PATH)
    data = data.dropna(how="all").copy()

    missing_columns = [column for column in RAW_COLUMNS if column not in data.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    data = data[RAW_COLUMNS].copy()
    before_drop = len(data)
    data = data.dropna(subset=RAW_COLUMNS).copy()
    dropped = before_drop - len(data)
    if dropped:
        warn(f"Dropped {dropped} records with missing required fields.")

    for column in NUMERIC_COLUMNS:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    before_numeric_drop = len(data)
    data = data.dropna(subset=NUMERIC_COLUMNS).copy()
    numeric_dropped = before_numeric_drop - len(data)
    if numeric_dropped:
        warn(f"Dropped {numeric_dropped} records with non-numeric indicator values.")

    for column in NUMERIC_COLUMNS:
        if (data[column] < 0).any():
            raise ValueError(f"Negative values detected in {column}.")

    if (data["Viral Sequence Count"] <= 0).any():
        raise ValueError("Viral Sequence Count must be greater than zero.")

    if data["Host"].duplicated().any():
        duplicated_hosts = data.loc[data["Host"].duplicated(), "Host"].tolist()
        raise ValueError(f"Duplicated Host values detected: {duplicated_hosts}")

    for column in NUMERIC_COLUMNS:
        if np.allclose(data[column], np.round(data[column])):
            data[column] = data[column].astype(int)

    if len(data) != EXPECTED_HOST_RECORDS:
        warn(
            f"Expected {EXPECTED_HOST_RECORDS} host records after cleaning, "
            f"but found {len(data)}."
        )

    return data


def calculate_entropy_weights(normalized_data: pd.DataFrame) -> pd.DataFrame:
    matrix = normalized_data[NORMALIZED_COLUMNS].to_numpy(dtype=float)
    n_records = matrix.shape[0]
    proportions = np.zeros_like(matrix)

    column_sums = matrix.sum(axis=0)
    for index, column_sum in enumerate(column_sums):
        if column_sum == 0:
            warn(
                f"{NORMALIZED_COLUMNS[index]} has a zero column sum; "
                "a uniform proportion vector was used for entropy calculation."
            )
            proportions[:, index] = 1 / n_records
        else:
            proportions[:, index] = matrix[:, index] / column_sum

    proportions_for_log = np.where(proportions == 0, 1e-12, proportions)
    entropy_coefficient = 1 / np.log(n_records)
    entropy = -entropy_coefficient * np.sum(
        proportions_for_log * np.log(proportions_for_log), axis=0
    )
    redundancy = 1 - entropy
    redundancy = np.where(redundancy < 0, 0, redundancy)

    redundancy_sum = redundancy.sum()
    if redundancy_sum == 0:
        warn("All redundancy values are zero; equal weights were used.")
        weights = np.full(len(NORMALIZED_COLUMNS), 1 / len(NORMALIZED_COLUMNS))
    else:
        weights = redundancy / redundancy_sum

    return pd.DataFrame(
        {
            "indicator": INDICATOR_LABELS,
            "entropy": entropy,
            "redundancy": redundancy,
            "weight": weights,
        }
    )


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    data = prepare_input()
    data["VD_norm"] = min_max_normalize(data["Viral Diversity"])
    data["ZSC_norm"] = min_max_normalize(data["Zoonotic Species Count"])
    data["ZSN_corrected"] = data["Zoonotic Sequence Number"] / np.log(
        data["Viral Sequence Count"] + 1
    )
    data["ZSN_corrected_norm"] = min_max_normalize(data["ZSN_corrected"])
    data["VS_ln"] = np.log(data["Viral Sequence Count"] + 1)
    data["VS_ln_norm"] = min_max_normalize(data["VS_ln"])

    weights = calculate_entropy_weights(data)
    weight_values = weights["weight"].to_numpy(dtype=float)
    data["EWM_CRS"] = data[NORMALIZED_COLUMNS].to_numpy(dtype=float).dot(weight_values)
    data = data.sort_values(["EWM_CRS", "Host"], ascending=[False, True]).copy()
    data["EWM_rank"] = np.arange(1, len(data) + 1)

    scores_output = data[
        [
            "Host",
            "Viral Diversity",
            "Zoonotic Species Count",
            "Zoonotic Sequence Number",
            "Viral Sequence Count",
            "VD_norm",
            "ZSC_norm",
            "ZSN_corrected",
            "ZSN_corrected_norm",
            "VS_ln",
            "VS_ln_norm",
            "EWM_CRS",
            "EWM_rank",
        ]
    ].copy()

    write_rounded_csv(weights, WEIGHTS_PATH)
    write_rounded_csv(scores_output, SCORES_PATH)

    print(f"Processed host records: {len(data)}")
    print(f"Entropy weight sum: {weights['weight'].sum():.{OUTPUT_DECIMAL_PLACES}f}")
    print(f"Output saved: {WEIGHTS_PATH.relative_to(PROJECT_ROOT)}")
    print(f"Output saved: {SCORES_PATH.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
