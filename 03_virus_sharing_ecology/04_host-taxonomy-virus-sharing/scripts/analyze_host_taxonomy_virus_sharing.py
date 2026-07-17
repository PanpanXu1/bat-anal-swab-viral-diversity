"""Compare viral-cluster sharing within and between host taxonomic groups."""

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from math import erfc, sqrt


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "input"
RESULTS_DIR = ROOT / "output" / "tables"



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


def mann_whitney_greater(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    x_values = x.to_numpy(dtype=float)
    y_values = y.to_numpy(dtype=float)
    n1, n2 = len(x_values), len(y_values)
    combined = np.concatenate([x_values, y_values])
    ranks = rank_average(combined)
    rank_sum_x = ranks[:n1].sum()
    u_statistic = rank_sum_x - n1 * (n1 + 1) / 2
    mean_u = n1 * n2 / 2
    sd_u = sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (u_statistic - mean_u - 0.5) / sd_u if sd_u else 0.0
    p_value = 0.5 * erfc(z / sqrt(2))
    return float(u_statistic), float(p_value)

def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    membership = pd.read_csv(DATA_DIR / "host_taxonomy_viral_cluster_membership.csv")
    taxonomy = pd.read_csv(DATA_DIR / "host_taxonomy.csv")
    required_membership = {"Cluster_ID", "host"}
    required_taxonomy = {"host", "family", "genus"}
    if not required_membership.issubset(membership.columns):
        raise ValueError(f"Membership table must contain {sorted(required_membership)}")
    if not required_taxonomy.issubset(taxonomy.columns):
        raise ValueError(f"Taxonomy table must contain {sorted(required_taxonomy)}")

    membership = membership.dropna(subset=list(required_membership)).drop_duplicates()
    taxonomy = taxonomy.dropna(subset=list(required_taxonomy)).drop_duplicates()
    if taxonomy["host"].duplicated().any():
        raise ValueError("Each host must have exactly one family and genus assignment")
    membership_hosts = set(membership["host"])
    missing = sorted(membership_hosts - set(taxonomy["host"]))
    if missing:
        raise ValueError(f"Taxonomy is missing analyzed hosts: {missing}")
    taxonomy = taxonomy[taxonomy["host"].isin(membership_hosts)].copy()
    return membership, taxonomy.set_index("host")


def calculate_pairwise_table(
    membership: pd.DataFrame,
    taxonomy: pd.DataFrame,
) -> pd.DataFrame:
    host_clusters = membership.groupby("host")["Cluster_ID"].apply(set).to_dict()
    records = []
    for host_1, host_2 in combinations(sorted(host_clusters), 2):
        shared = len(host_clusters[host_1] & host_clusters[host_2])
        records.append(
            {
                "Host1": host_1,
                "Host2": host_2,
                "Family1": taxonomy.loc[host_1, "family"],
                "Family2": taxonomy.loc[host_2, "family"],
                "Genus1": taxonomy.loc[host_1, "genus"],
                "Genus2": taxonomy.loc[host_2, "genus"],
                "Same_family": taxonomy.loc[host_1, "family"]
                == taxonomy.loc[host_2, "family"],
                "Same_genus": taxonomy.loc[host_1, "genus"]
                == taxonomy.loc[host_2, "genus"],
                "Shared_clusters": shared,
                "Log10_shared_clusters_plus1": np.log10(shared + 1),
            }
        )
    return pd.DataFrame(records)


def test_level(pairs: pd.DataFrame, level: str) -> dict[str, object]:
    flag = f"Same_{level.lower()}"
    same = pairs.loc[pairs[flag], "Log10_shared_clusters_plus1"]
    different = pairs.loc[~pairs[flag], "Log10_shared_clusters_plus1"]
    u_statistic, p_value = mann_whitney_greater(same, different)
    return {
        "Taxonomic_level": level,
        "Alternative_hypothesis": f"Same {level} > Different {level}",
        "N_same_group_pairs": len(same),
        "N_different_group_pairs": len(different),
        "U_statistic": u_statistic,
        "P_value": p_value,
        "Same_group_mean_log10_plus1": same.mean(),
        "Different_group_mean_log10_plus1": different.mean(),
        "Same_group_median_log10_plus1": same.median(),
        "Different_group_median_log10_plus1": different.median(),
    }


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    membership, taxonomy = load_inputs()
    pairs = calculate_pairwise_table(membership, taxonomy)
    expected_pairs = len(taxonomy) * (len(taxonomy) - 1) // 2
    if len(pairs) != expected_pairs:
        raise ValueError(f"Expected {expected_pairs} unique nonself pairs; found {len(pairs)}")

    results = [test_level(pairs, "Family"), test_level(pairs, "Genus")]
    pairs.to_csv(
        RESULTS_DIR / "host_pair_taxonomy_and_sharing.csv",
        index=False,
        float_format="%.10g",
    )
    pd.DataFrame(results).to_csv(
        RESULTS_DIR / "mann_whitney_u_results.csv",
        index=False,
        float_format="%.10g",
    )

    print(f"Hosts: {len(taxonomy)}; unique nonself pairs: {len(pairs)}")
    for result in results:
        print(
            f"{result['Taxonomic_level']}: U = {result['U_statistic']:.1f}; "
            f"P = {result['P_value']:.6g}"
        )


if __name__ == "__main__":
    main()
