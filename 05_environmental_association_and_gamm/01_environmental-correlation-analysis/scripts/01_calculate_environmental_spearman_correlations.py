"""Calculate pairwise Spearman correlations among environmental variables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "input" / "shannon_environmental_analysis_input.csv"
DEFAULT_OUTPUT = ROOT / "output" / "tables" / "environmental_spearman"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate pairwise Spearman correlations among environmental variables."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()



def rank_average(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2 + 1
        start = end
    return ranks


def spearman_statistic(x: pd.Series, y: pd.Series) -> float:
    rx = rank_average(x.to_numpy(dtype=float))
    ry = rank_average(y.to_numpy(dtype=float))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt(np.sum(rx ** 2) * np.sum(ry ** 2))
    return float(np.sum(rx * ry) / denom) if denom else float("nan")


def permutation_p_value(x: pd.Series, y: pd.Series, observed: float, permutations: int = 999, seed: int = 42) -> float:
    rng = np.random.default_rng(seed)
    y_values = y.to_numpy(dtype=float)
    extreme = 0
    for _ in range(permutations):
        statistic = spearman_statistic(x, pd.Series(rng.permutation(y_values)))
        if abs(statistic) >= abs(observed):
            extreme += 1
    return (extreme + 1) / (permutations + 1)

def load_environmental_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    if "Sample" not in data.columns:
        raise ValueError("Input must contain a 'Sample' column.")
    if data["Sample"].isna().any() or data["Sample"].duplicated().any():
        raise ValueError("Sample identifiers must be complete and unique.")

    if "Shannon" not in data.columns:
        raise ValueError("Input must contain a 'Shannon' column.")

    variables = data.drop(columns=["Sample", "Shannon"]).apply(
        pd.to_numeric, errors="raise"
    )
    if variables.isna().any().any():
        raise ValueError("Environmental variables must not contain missing values.")
    return variables


def correlation_long_table(data: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, float | str]] = []
    columns = list(data.columns)
    for index, factor_1 in enumerate(columns):
        for factor_2 in columns[index + 1 :]:
            rho = spearman_statistic(data[factor_1], data[factor_2])
            p_value = permutation_p_value(data[factor_1], data[factor_2], rho, seed=42 + index)
            records.append(
                {
                    "Factor1": factor_1,
                    "Factor2": factor_2,
                    "Spearman_rho": rho,
                    "p_value": p_value,
                }
            )
    return pd.DataFrame(records)


def main() -> None:
    args = parse_args()
    data = load_environmental_data(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    correlation_matrix = data.corr(method="spearman")
    correlation_matrix.to_csv(
        args.output_dir / "environmental_spearman_correlation_matrix.csv"
    )

    long_table = correlation_long_table(data)
    long_table.to_csv(
        args.output_dir / "environmental_spearman_pairwise_tests.csv",
        index=False,
    )

    print(f"Samples analyzed: {len(data)}")
    print(f"Environmental variables analyzed: {data.shape[1]}")


if __name__ == "__main__":
    main()
