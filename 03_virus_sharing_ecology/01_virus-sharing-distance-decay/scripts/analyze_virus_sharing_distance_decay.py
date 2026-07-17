"""Run the virus-sharing distance-decay analysis."""

from itertools import combinations
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "input"
TABLES_DIR = ROOT / "output" / "tables"
FIGURES_DIR = ROOT / "output" / "figures"


def configure_fonts() -> None:
    for font_file in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
        font_path = Path("C:/Windows/Fonts") / font_file
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )



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


def permutation_p_value(x: pd.Series, y: pd.Series, observed: float, permutations: int = 9999, seed: int = 42) -> float:
    rng = np.random.default_rng(seed)
    y_values = y.to_numpy(dtype=float)
    extreme = 0
    for _ in range(permutations):
        statistic = spearman_statistic(x, pd.Series(rng.permutation(y_values)))
        if abs(statistic) >= abs(observed):
            extreme += 1
    return max((extreme + 1) / (permutations + 1), 1 / (permutations + 1))

def canonical_pair(sample_1: str, sample_2: str) -> tuple[str, str]:
    return tuple(sorted((sample_1, sample_2)))


def calculate_shared_clusters(membership: pd.DataFrame) -> pd.DataFrame:
    sample_clusters = (
        membership.groupby("SampleGroup", sort=True)["Cluster_ID"].apply(set).to_dict()
    )
    records = []
    for sample_1, sample_2 in combinations(sorted(sample_clusters), 2):
        shared_count = len(sample_clusters[sample_1] & sample_clusters[sample_2])
        if shared_count > 0:
            records.append((sample_1, sample_2, shared_count))
    return pd.DataFrame(
        records, columns=["SampleGroup1", "SampleGroup2", "SharedClusters"]
    )


def main() -> None:
    configure_fonts()
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    membership = pd.read_csv(DATA_DIR / "distance_decay_viral_cluster_membership.csv")
    distances = pd.read_csv(DATA_DIR / "sampling_site_pair_distances.csv")

    required_membership = {"Sequence_ID", "Cluster_ID", "SampleGroup"}
    required_distances = {"SampleGroup1", "SampleGroup2", "Distance_km"}
    if not required_membership.issubset(membership.columns):
        raise ValueError(f"Membership table must contain {sorted(required_membership)}")
    if not required_distances.issubset(distances.columns):
        raise ValueError(f"Distance table must contain {sorted(required_distances)}")

    membership = membership.dropna(subset=list(required_membership)).drop_duplicates()
    distances[["SampleGroup1", "SampleGroup2"]] = distances.apply(
        lambda row: canonical_pair(row["SampleGroup1"], row["SampleGroup2"]),
        axis=1,
        result_type="expand",
    )
    distances = distances.drop_duplicates(["SampleGroup1", "SampleGroup2"])

    shared = calculate_shared_clusters(membership)
    analysis = shared.merge(
        distances, on=["SampleGroup1", "SampleGroup2"], how="left", validate="one_to_one"
    )
    if analysis["Distance_km"].isna().any():
        missing = analysis.loc[
            analysis["Distance_km"].isna(), ["SampleGroup1", "SampleGroup2"]
        ]
        raise ValueError(f"Missing distances for {len(missing)} site pairs")

    analysis["Log10SharedClusters"] = np.log10(analysis["SharedClusters"])
    rho = spearman_statistic(analysis["Distance_km"], analysis["SharedClusters"])
    p_value = permutation_p_value(analysis["Distance_km"], analysis["SharedClusters"], rho)

    analysis = analysis.sort_values(["SampleGroup1", "SampleGroup2"]).reset_index(drop=True)
    analysis.to_csv(
        TABLES_DIR / "shared_clusters_distance.csv",
        index=False,
        float_format="%.6f",
    )
    pd.DataFrame(
        [
            {
                "Method": "Spearman",
                "N_site_pairs": len(analysis),
                "Spearman_rho": rho,
                "P_value": p_value,
                "Included_pairs": "SharedClusters > 0",
            }
        ]
    ).to_csv(
        TABLES_DIR / "spearman_correlation.csv",
        index=False,
        float_format="%.10g",
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(analysis["Distance_km"], analysis["Log10SharedClusters"], alpha=0.65, s=9, color="#4f79a7")
    slope, intercept = np.polyfit(analysis["Distance_km"], analysis["Log10SharedClusters"], 1)
    x_line = np.linspace(analysis["Distance_km"].min(), analysis["Distance_km"].max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="#1f5bd7", linewidth=2)
    ax.grid(color="#E6E6E6", linewidth=0.8)
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Shared Clusters (Log10)")
    y_tick_step = 0.5
    y_upper = np.ceil(analysis["Log10SharedClusters"].max() / y_tick_step) * y_tick_step
    ax.set_ylim(0, y_upper)
    ax.yaxis.set_major_locator(MultipleLocator(y_tick_step))
    ax.text(
        0.98,
        0.97,
        rf"Spearman's $\rho$ = {rho:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
    )
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "distance_decay_log10.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Analyzed {len(analysis)} site pairs with SharedClusters > 0")
    print(f"Spearman's rho = {rho:.4f}; P = {p_value:.4e}")


if __name__ == "__main__":
    main()
