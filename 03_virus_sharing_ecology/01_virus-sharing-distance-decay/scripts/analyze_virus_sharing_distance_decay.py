"""Analyze positive viral sharing using sampling-group label permutations."""

from itertools import combinations
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "input"
TABLES_DIR = ROOT / "output" / "tables"
FIGURES_DIR = ROOT / "output" / "figures"
PERMUTATIONS = 9999
SEED = 42


def configure_fonts() -> None:
    for font_file in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
        font_path = Path("C:/Windows/Fonts") / font_file
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
    plt.rcParams.update({
        "font.family": "Times New Roman",
        "font.serif": ["Times New Roman"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def rank_average(values: np.ndarray) -> np.ndarray:
    _, inverse, counts = np.unique(values, return_inverse=True, return_counts=True)
    ends = np.cumsum(counts)
    average_ranks = (2 * ends - counts + 1) / 2
    return average_ranks[inverse]


def spearman_statistic(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray) -> float:
    rx = rank_average(np.asarray(x, dtype=float))
    ry = rank_average(np.asarray(y, dtype=float))
    rx -= rx.mean()
    ry -= ry.mean()
    denominator = np.sqrt(np.sum(rx ** 2) * np.sum(ry ** 2))
    return float(np.sum(rx * ry) / denominator) if denominator else float("nan")


def canonical_pair(sample_1: str, sample_2: str) -> tuple[str, str]:
    return tuple(sorted((sample_1, sample_2)))


def calculate_shared_clusters(membership: pd.DataFrame) -> pd.DataFrame:
    sample_clusters = membership.groupby("SampleGroup", sort=True)["Cluster_ID"].apply(set).to_dict()
    records = []
    for first, second in combinations(sorted(sample_clusters), 2):
        count = len(sample_clusters[first] & sample_clusters[second])
        if count > 0:
            records.append((first, second, count))
    return pd.DataFrame(records, columns=["SampleGroup1", "SampleGroup2", "SharedClusters"])


def load_analysis() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Load observed positive pairs and complete, identically ordered matrices."""
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
        axis=1, result_type="expand",
    )
    if distances.duplicated(["SampleGroup1", "SampleGroup2"]).any():
        raise ValueError("Distance table must have one row per unordered sampling-group pair")
    if not np.isfinite(distances["Distance_km"]).all() or (distances["Distance_km"] < 0).any():
        raise ValueError("Distances must be finite, non-negative values in kilometres")

    groups = sorted(membership["SampleGroup"].unique())
    all_pairs = pd.DataFrame(combinations(groups, 2), columns=["SampleGroup1", "SampleGroup2"])
    complete = all_pairs.merge(
        distances, on=["SampleGroup1", "SampleGroup2"], how="left", validate="one_to_one",
    )
    missing = int(complete["Distance_km"].isna().sum())
    if missing or len(distances) != len(all_pairs):
        raise ValueError(
            f"Label permutations require distances for all {len(all_pairs)} sampling-group pairs; "
            f"{missing} required distances are missing"
        )

    group_index = {group: index for index, group in enumerate(groups)}
    upper = np.triu_indices(len(groups), 1)
    distance_matrix = np.zeros((len(groups), len(groups)), dtype=float)
    distance_matrix[upper] = complete["Distance_km"].to_numpy()
    distance_matrix[(upper[1], upper[0])] = distance_matrix[upper]

    shared = calculate_shared_clusters(membership)
    rows = shared["SampleGroup1"].map(group_index).to_numpy()
    columns = shared["SampleGroup2"].map(group_index).to_numpy()
    sharing_matrix = np.zeros((len(groups), len(groups)), dtype=int)
    sharing_matrix[rows, columns] = shared["SharedClusters"].to_numpy()
    sharing_matrix[columns, rows] = shared["SharedClusters"].to_numpy()
    analysis = shared.merge(
        distances, on=["SampleGroup1", "SampleGroup2"], how="left", validate="one_to_one",
    ).sort_values(["SampleGroup1", "SampleGroup2"]).reset_index(drop=True)
    analysis["Log10SharedClusters"] = np.log10(analysis["SharedClusters"])
    return analysis, sharing_matrix, distance_matrix


def site_label_permutation(
    sharing_matrix: np.ndarray,
    distance_matrix: np.ndarray,
    thresholds: tuple[int | None, ...] = (None,),
    permutations: int = PERMUTATIONS,
    seed: int = SEED,
) -> pd.DataFrame:
    """Test absolute Spearman correlations under exchangeable group labels.

    One permutation reorders both axes of the sharing matrix with the same
    index vector. Geography stays fixed. Positive sharing and each geographic
    threshold are applied anew, with average ranks recomputed for every subset.
    Only the upper triangle is used; self-pairs are excluded.
    """
    rows, columns = np.triu_indices(len(sharing_matrix), 1)
    distances = distance_matrix[rows, columns]
    counts = sharing_matrix[rows, columns]
    distance_masks = [
        np.ones(len(distances), dtype=bool) if threshold is None else distances > threshold
        for threshold in thresholds
    ]
    observed_masks = [mask & (counts > 0) for mask in distance_masks]
    observed = np.array([
        spearman_statistic(distances[mask], counts[mask]) for mask in observed_masks
    ])
    if not np.isfinite(observed).all():
        raise ValueError("Each observed subset needs variation in distances and shared counts")
    extreme = np.zeros(len(thresholds), dtype=int)
    rng = np.random.default_rng(seed)
    for iteration in range(permutations):
        order = rng.permutation(len(sharing_matrix))
        # These entries equal sharing_matrix[np.ix_(order, order)][rows, columns].
        permuted_counts = sharing_matrix[order[rows], order[columns]]
        positive = permuted_counts > 0
        for index, distance_mask in enumerate(distance_masks):
            mask = distance_mask & positive
            statistic = spearman_statistic(distances[mask], permuted_counts[mask])
            if not np.isfinite(statistic):
                raise ValueError("A permuted subset has an undefined Spearman correlation")
            if abs(statistic) >= abs(observed[index]) - 1e-12:
                extreme[index] += 1
        if (iteration + 1) % 2000 == 0:
            print(f"Completed {iteration + 1}/{permutations} sampling-group label permutations", flush=True)

    return pd.DataFrame({
        "Distance_filter": ["All" if threshold is None else f">{threshold} km" for threshold in thresholds],
        "N_site_pairs": [int(mask.sum()) for mask in observed_masks],
        "Spearman_rho": observed,
        "P_value": (extreme + 1) / (permutations + 1),
        "Extreme_permutations": extreme,
        "Permutations": permutations,
    })


def format_p_value(p_value: float) -> str:
    if p_value < 0.001:
        return "<0.001"
    if p_value < 0.01:
        return "<0.01"
    return f"{p_value:.2f}"


def format_summary(summary: pd.DataFrame) -> pd.DataFrame:
    display = summary.copy()
    display["Spearman_rho"] = display["Spearman_rho"].map(lambda value: f"{value:.2f}")
    display["P_value"] = display["P_value"].map(format_p_value)
    display["Permutation unit"] = "Sampling-group labels"
    return display.rename(columns={
        "Distance_filter": "Distance filter",
        "N_site_pairs": "Number of sampling-group pairs",
        "Spearman_rho": "Spearman's rank correlation coefficient",
        "P_value": "P value",
        "Extreme_permutations": "Number of extreme permutations",
        "Permutations": "Number of permutations",
    })


def plot_distance_decay(
    analysis: pd.DataFrame, rho: float, p_value: float, threshold: int | None,
) -> None:
    configure_fonts()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(analysis["Distance_km"], analysis["Log10SharedClusters"], alpha=0.65, s=9, color="#4f79a7")
    slope, intercept = np.polyfit(analysis["Distance_km"], analysis["Log10SharedClusters"], 1)
    x_line = np.linspace(analysis["Distance_km"].min(), analysis["Distance_km"].max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="#1f5bd7", linewidth=2)
    ax.grid(color="#E6E6E6", linewidth=0.8)
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Shared Clusters (Log10)")
    ax.set_ylim(0, np.ceil(analysis["Log10SharedClusters"].max() / 0.5) * 0.5)
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    if threshold == 0:
        minimum = analysis["Distance_km"].min()
        right_limit = ax.get_xlim()[1]
        ax.set_xticks([minimum, *[tick for tick in ax.get_xticks() if minimum < tick <= right_limit]])
    p_text = format_p_value(p_value)
    relation = "" if p_text.startswith("<") else "= "
    ax.text(
        0.98, 0.97, rf"Spearman's $\rho$ = {rho:.2f}" + "\n" + rf"$P$ {relation}{p_text}",
        transform=ax.transAxes, ha="right", va="top",
    )
    label = "No distance cutoff" if threshold is None else f"Distance > {threshold} km"
    ax.text(0.98, 0.035, label, transform=ax.transAxes, ha="right", va="bottom")
    fig.tight_layout()
    name = "distance_decay_log10" if threshold is None else f"distance_decay_over_{threshold}km"
    fig.savefig(FIGURES_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    analysis, sharing_matrix, distance_matrix = load_analysis()
    summary = site_label_permutation(sharing_matrix, distance_matrix)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    analysis.to_csv(TABLES_DIR / "shared_clusters_distance.csv", index=False, float_format="%.2f")
    format_summary(summary).to_csv(TABLES_DIR / "spearman_correlation.csv", index=False)
    row = summary.iloc[0]
    plot_distance_decay(analysis, row["Spearman_rho"], row["P_value"], None)
    print(f"All: n={len(analysis)}, rho={row['Spearman_rho']:.2f}, P {format_p_value(row['P_value'])}")


if __name__ == "__main__":
    main()
