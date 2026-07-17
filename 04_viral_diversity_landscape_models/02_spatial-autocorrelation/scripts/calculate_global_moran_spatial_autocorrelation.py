"""Run global spatial autocorrelation analysis of viral Shannon diversity."""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "input" / "viral_shannon_coordinates.csv"
TABLE_DIR = ROOT / "output" / "tables"

K_NEIGHBORS = 5
PERMUTATIONS = 9999
RANDOM_SEED = 42
REQUIRED_COLUMNS = {"SampleGroup", "Shannon_Index", "lat", "long"}


def load_samples(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")
    if data[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ValueError("Required input columns contain missing values.")
    if data["SampleGroup"].duplicated().any():
        raise ValueError("SampleGroup values must be unique.")
    return data.copy()


def knn_weights(coordinates: np.ndarray, k: int) -> np.ndarray:
    if k >= len(coordinates):
        raise ValueError("K_NEIGHBORS must be smaller than the number of samples.")
    diff = coordinates[:, None, :] - coordinates[None, :, :]
    distances = np.sqrt(np.sum(diff * diff, axis=2))
    np.fill_diagonal(distances, np.inf)
    neighbors = np.argsort(distances, axis=1)[:, :k]
    weights = np.zeros((len(coordinates), len(coordinates)), dtype=float)
    row_index = np.arange(len(coordinates))[:, None]
    weights[row_index, neighbors] = 1.0 / k
    return weights


def moran_i(values: np.ndarray, weights: np.ndarray) -> float:
    centered = values - values.mean()
    numerator = np.sum(weights * np.outer(centered, centered))
    denominator = np.sum(centered ** 2)
    weight_sum = weights.sum()
    return float(len(values) / weight_sum * numerator / denominator)


def permutation_p_value(values: np.ndarray, weights: np.ndarray, observed: float) -> tuple[float, str]:
    rng = np.random.default_rng(RANDOM_SEED)
    extreme = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(values)
        statistic = moran_i(permuted, weights)
        if abs(statistic) >= abs(observed):
            extreme += 1
    if extreme == 0:
        lower = 1 / (PERMUTATIONS + 1)
        return lower, "<0.0001"
    p_value = extreme / PERMUTATIONS
    return p_value, f"{p_value:.4g}"


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    samples = load_samples(INPUT_PATH)
    coordinates = samples[["long", "lat"]].to_numpy(dtype=float)
    values = samples["Shannon_Index"].to_numpy(dtype=float)
    weights = knn_weights(coordinates, K_NEIGHBORS)
    observed = moran_i(values, weights)
    expected = -1 / (len(values) - 1)
    p_value, p_display = permutation_p_value(values, weights, observed)

    global_summary = pd.DataFrame(
        {
            "n_samples": [len(samples)],
            "k_neighbors": [K_NEIGHBORS],
            "permutations": [PERMUTATIONS],
            "random_seed": [RANDOM_SEED],
            "global_moran_i": [observed],
            "expected_moran_i": [expected],
            "permutation_p_value": [p_value],
            "permutation_p_value_display": [p_display],
        }
    )
    global_summary.to_csv(TABLE_DIR / "global_moran_summary.csv", index=False)
    print(global_summary.to_string(index=False))


if __name__ == "__main__":
    main()
