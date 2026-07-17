"""Test the association between terrain-pair type and shared viral clusters."""

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = ROOT / "input" / "terrain_viral_cluster_membership.csv"
TABLES_DIR = ROOT / "output" / "tables"
PERMUTATIONS = 999
RANDOM_SEED = 42


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


def one_way_f(values: np.ndarray, groups: np.ndarray) -> float:
    overall = values.mean()
    labels = np.unique(groups)
    ss_between = sum(np.sum(groups == label) * (values[groups == label].mean() - overall) ** 2 for label in labels)
    ss_within = sum(np.sum((values[groups == label] - values[groups == label].mean()) ** 2) for label in labels)
    df_between = len(labels) - 1
    df_within = len(values) - len(labels)
    return float((ss_between / df_between) / (ss_within / df_within)) if ss_within else float("inf")


def kruskal_h(values: np.ndarray, groups: np.ndarray) -> float:
    ranks = rank_average(values)
    labels = np.unique(groups)
    n = len(values)
    h = 12 / (n * (n + 1)) * sum((ranks[groups == label].sum() ** 2) / np.sum(groups == label) for label in labels) - 3 * (n + 1)
    return float(h)


def permutation_p(values: np.ndarray, groups: np.ndarray, observed: float, statistic_fn) -> float:
    rng = np.random.default_rng(RANDOM_SEED)
    extreme = 0
    for _ in range(PERMUTATIONS):
        statistic = statistic_fn(values, rng.permutation(groups))
        if statistic >= observed:
            extreme += 1
    return (extreme + 1) / (PERMUTATIONS + 1)


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA_FILE).dropna().drop_duplicates()
    sample_clusters = data.groupby("SampleGroup")["Cluster_ID"].apply(set).to_dict()
    sample_terrain = data.drop_duplicates("SampleGroup").set_index("SampleGroup")["Terrain"].to_dict()

    records = []
    for sample_1, sample_2 in combinations(sorted(sample_clusters), 2):
        terrain_pair = "-".join(sorted((sample_terrain[sample_1], sample_terrain[sample_2])))
        records.append({
            "SampleGroup1": sample_1,
            "SampleGroup2": sample_2,
            "Terrain1": sample_terrain[sample_1],
            "Terrain2": sample_terrain[sample_2],
            "TerrainPair": terrain_pair,
            "SharedClusters": len(sample_clusters[sample_1] & sample_clusters[sample_2]),
        })
    pairs = pd.DataFrame(records)
    pairs.to_csv(TABLES_DIR / "terrain_pairwise_shared_clusters.csv", index=False)

    values = pairs["SharedClusters"].to_numpy(dtype=float)
    groups = pairs["TerrainPair"].to_numpy()
    f_stat = one_way_f(values, groups)
    f_p = permutation_p(values, groups, f_stat, one_way_f)
    h_stat = kruskal_h(values, groups)
    h_p = permutation_p(values, groups, h_stat, kruskal_h)

    pd.DataFrame([{
        "Method": "Permutation one-way ANOVA-style F test",
        "Grouping_variable": "TerrainPair",
        "N_unique_nonself_site_pairs": len(pairs),
        "F_statistic": f_stat,
        "P_value": f_p,
        "Permutations": PERMUTATIONS,
        "Random_seed": RANDOM_SEED,
    }]).to_csv(TABLES_DIR / "terrain_anova_results.csv", index=False, float_format="%.10g")
    pd.DataFrame([{
        "Method": "Permutation Kruskal-Wallis rank test",
        "Grouping_variable": "TerrainPair",
        "N_unique_nonself_site_pairs": len(pairs),
        "H_statistic": h_stat,
        "P_value": h_p,
        "Permutations": PERMUTATIONS,
        "Random_seed": RANDOM_SEED,
    }]).to_csv(TABLES_DIR / "terrain_kruskal_wallis_results.csv", index=False, float_format="%.10g")

    print(f"Unique non-self site pairs: {len(pairs)}")
    print(f"Permutation F-test P = {f_p:.4g}")
    print(f"Permutation Kruskal-Wallis P = {h_p:.4g}")


if __name__ == "__main__":
    main()
