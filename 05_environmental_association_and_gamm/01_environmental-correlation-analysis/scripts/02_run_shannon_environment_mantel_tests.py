"""Test associations between Shannon-diversity and environmental-distance matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "input" / "shannon_environmental_analysis_input.csv"
DEFAULT_OUTPUT = ROOT / "output" / "tables" / "shannon_environment_mantel"

# Keep permutation streams stable when input columns are reorganized.
MANTEL_VARIABLE_ORDER = [
    "DEM",
    "HFT",
    "Bio19",
    "Bio1",
    "Bio2",
    "Bio3",
    "Bio4",
    "Bio5",
    "Bio6",
    "Bio7",
    "Bio8",
    "Bio9",
    "Bio10",
    "Bio11",
    "Bio12",
    "Bio13",
    "Bio14",
    "Bio15",
    "Bio16",
    "Bio17",
    "Bio18",
    "GMR",
    "CHEQ",
    "NDVI",
    "FVC",
    "GDP",
    "PSD",
    "GR",
    "GLH",
]



def pairwise_euclidean_1d(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    rows, cols = np.triu_indices(len(values), k=1)
    return np.abs(values[rows] - values[cols])


def pearson_statistic(x: np.ndarray, y: np.ndarray) -> float:
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt(np.sum(x ** 2) * np.sum(y ** 2))
    return float(np.sum(x * y) / denom) if denom else float("nan")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Perform Pearson Mantel tests between pairwise Euclidean distances "
            "in Shannon diversity and each environmental variable."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def load_analysis_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    if "Sample" not in data.columns:
        raise ValueError("Input must contain a 'Sample' column.")
    if "Shannon" not in data.columns:
        raise ValueError("Input must contain a 'Shannon' column.")
    if data["Sample"].isna().any() or data["Sample"].duplicated().any():
        raise ValueError("Sample identifiers must be complete and unique.")

    numeric = data.drop(columns="Sample").apply(pd.to_numeric, errors="raise")
    if numeric.isna().any().any():
        raise ValueError("Shannon and environmental variables must not contain missing values.")
    return data


def mantel_test(
    response: np.ndarray,
    predictor: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    response_distance = pairwise_euclidean_1d(response)
    predictor_distance = pairwise_euclidean_1d(predictor)
    observed = pearson_statistic(response_distance, predictor_distance)

    exceedances = 0
    for _ in range(permutations):
        permuted_distance = pairwise_euclidean_1d(rng.permutation(predictor))
        statistic = pearson_statistic(response_distance, permuted_distance)
        exceedances += abs(statistic) >= abs(observed)

    p_value = (exceedances + 1) / (permutations + 1)
    return observed, p_value


def main() -> None:
    args = parse_args()
    if args.permutations < 1:
        raise ValueError("--permutations must be at least 1.")

    data = load_analysis_data(args.input)
    response = data["Shannon"].to_numpy(dtype=float)
    variables = [column for column in data.columns if column not in {"Sample", "Shannon"}]
    if set(variables) != set(MANTEL_VARIABLE_ORDER):
        raise ValueError(
            "Environmental columns must match the documented 29-variable input schema."
        )
    variables = MANTEL_VARIABLE_ORDER
    rng = np.random.default_rng(args.seed)

    records = []
    for variable in variables:
        statistic, p_value = mantel_test(
            response,
            data[variable].to_numpy(dtype=float),
            args.permutations,
            rng,
        )
        records.append(
            {
                "Response": "Shannon",
                "Environmental_variable": variable,
                "Mantel_r": statistic,
                "p_value": p_value,
                "permutations": args.permutations,
                "random_seed": args.seed,
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(
        args.output_dir / "shannon_environment_mantel_tests.csv",
        index=False,
    )

    print(f"Samples analyzed: {len(data)}")
    print(f"Environmental variables analyzed: {len(variables)}")
    print(f"Permutations per test: {args.permutations}")


if __name__ == "__main__":
    main()
