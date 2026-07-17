"""Analyze viral-cluster sharing in relation to host sequence identity."""

from itertools import combinations
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "input"
TABLES_DIR = ROOT / "output" / "tables"
FIGURES_DIR = ROOT / "output" / "figures"
N_PERMUTATIONS = 9_999
RANDOM_SEED = 42
SPECIES_TICK_LABEL_SIZE = 14
COLORBAR_LABEL_SIZE = 14
COLORBAR_TICK_LABEL_SIZE = 12


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


def format_species_name(name: str) -> str:
    return str(name).replace("_", " ")


def format_species_labels(data: pd.DataFrame) -> pd.DataFrame:
    return data.rename(index=format_species_name, columns=format_species_name)


def load_membership() -> pd.DataFrame:
    membership = pd.read_csv(DATA_DIR / "host_identity_viral_cluster_membership.csv")
    required = {"Cluster_ID", "host"}
    if not required.issubset(membership.columns):
        raise ValueError(f"Membership table must contain {sorted(required)}")
    membership = membership.dropna(subset=list(required)).drop_duplicates()
    if membership.empty:
        raise ValueError("Membership table contains no valid records")
    return membership


def load_identity() -> pd.DataFrame:
    identity = pd.read_csv(DATA_DIR / "host_sequence_identity_percent.csv", index_col=0)
    identity.index = identity.index.astype(str)
    identity.columns = identity.columns.astype(str)
    if identity.index.duplicated().any() or identity.columns.duplicated().any():
        raise ValueError("Host sequence identity matrix contains duplicate labels")
    if set(identity.index) != set(identity.columns):
        raise ValueError("Identity matrix row and column labels differ")
    identity = identity.loc[identity.index, identity.index].apply(pd.to_numeric)
    if identity.isna().any().any():
        raise ValueError("Identity matrix contains missing or nonnumeric values")
    if not np.allclose(identity.to_numpy(), identity.to_numpy().T):
        raise ValueError("Identity matrix is not symmetric")
    if not np.allclose(np.diag(identity), 100):
        raise ValueError("Identity matrix diagonal must equal 100%")
    return identity


def calculate_shared_clusters(membership: pd.DataFrame, hosts: list[str]) -> pd.DataFrame:
    membership_hosts = set(membership["host"].astype(str))
    expected_hosts = set(hosts)
    if membership_hosts != expected_hosts:
        missing = sorted(expected_hosts - membership_hosts)
        extra = sorted(membership_hosts - expected_hosts)
        raise ValueError(f"Host sets differ; missing={missing}, extra={extra}")

    presence = pd.crosstab(membership["host"], membership["Cluster_ID"]).gt(0).astype(np.int16)
    shared = presence.dot(presence.T).reindex(index=hosts, columns=hosts).copy()
    for host in hosts:
        shared.loc[host, host] = 0
    if not np.allclose(shared.to_numpy(), shared.to_numpy().T):
        raise ValueError("Calculated sharing matrix is not symmetric")
    return shared



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


def spearman_statistic(x: np.ndarray, y: np.ndarray) -> float:
    rx = rank_average(np.asarray(x, dtype=float))
    ry = rank_average(np.asarray(y, dtype=float))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt(np.sum(rx ** 2) * np.sum(ry ** 2))
    return float(np.sum(rx * ry) / denom) if denom else float("nan")

def mantel_spearman(
    matrix_1: np.ndarray,
    matrix_2: np.ndarray,
    permutations: int,
    seed: int,
) -> tuple[float, float]:
    triangle = np.triu_indices_from(matrix_1, k=1)
    observed = spearman_statistic(matrix_1[triangle], matrix_2[triangle])
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(permutations):
        permutation = rng.permutation(matrix_2.shape[0])
        permuted = matrix_2[np.ix_(permutation, permutation)]
        statistic = spearman_statistic(matrix_1[triangle], permuted[triangle])
        extreme += abs(statistic) >= abs(observed)
    p_value = (extreme + 1) / (permutations + 1)
    return observed, p_value


def save_heatmap(data: pd.DataFrame, output: Path, colorbar_label: str) -> None:
    display_data = format_species_labels(data)
    values = display_data.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(17, 16))
    image = ax.imshow(values, cmap="coolwarm", aspect="equal")
    ax.set_xticks(np.arange(display_data.shape[1]))
    ax.set_yticks(np.arange(display_data.shape[0]))
    ax.set_xticklabels(display_data.columns)
    ax.set_yticklabels(display_data.index)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.72, pad=0.015)
    colorbar.ax.tick_params(labelsize=COLORBAR_TICK_LABEL_SIZE)
    colorbar.set_label(colorbar_label, fontsize=COLORBAR_LABEL_SIZE, fontfamily="Times New Roman")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=45, labelsize=SPECIES_TICK_LABEL_SIZE)
    ax.tick_params(axis="y", labelrotation=0, labelsize=SPECIES_TICK_LABEL_SIZE)
    for tick_label in ax.get_xticklabels():
        tick_label.set_horizontalalignment("right")
        tick_label.set_rotation_mode("anchor")
    for tick_label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        tick_label.set_fontfamily("Times New Roman")
        tick_label.set_fontstyle("italic")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

def main() -> None:
    configure_fonts()
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    membership = load_membership()
    identity = load_identity()
    shared = calculate_shared_clusters(membership, identity.index.tolist())
    shared_log = np.log10(shared + 1)
    for host in shared_log.index:
        shared_log.loc[host, host] = 0

    statistic, p_value = mantel_spearman(
        identity.to_numpy(),
        shared_log.to_numpy(),
        N_PERMUTATIONS,
        RANDOM_SEED,
    )

    format_species_labels(shared).to_csv(TABLES_DIR / "host_shared_cluster_matrix.csv")
    format_species_labels(shared_log).to_csv(
        TABLES_DIR / "host_shared_cluster_log10_plus1_matrix.csv",
        float_format="%.10g",
    )

    records = []
    for host_1, host_2 in combinations(identity.index, 2):
        records.append(
            {
                "Host1": format_species_name(host_1),
                "Host2": format_species_name(host_2),
                "Host_identity_percent": identity.loc[host_1, host_2],
                "Shared_clusters": int(shared.loc[host_1, host_2]),
                "Log10_shared_clusters_plus1": shared_log.loc[host_1, host_2],
            }
        )
    pd.DataFrame(records).to_csv(
        TABLES_DIR / "host_pair_identity_and_sharing.csv",
        index=False,
        float_format="%.10g",
    )
    pd.DataFrame(
        [
            {
                "Method": "Two-sided Spearman Mantel test",
                "N_hosts": len(identity),
                "N_unique_nonself_host_pairs": len(records),
                "Mantel_R": statistic,
                "P_value": p_value,
                "Permutations": N_PERMUTATIONS,
                "Random_seed": RANDOM_SEED,
                "Sharing_transformation": "log10(shared clusters + 1)",
            }
        ]
    ).to_csv(
        TABLES_DIR / "mantel_test_results.csv",
        index=False,
        float_format="%.10g",
    )

    save_heatmap(
        identity,
        FIGURES_DIR / "host_sequence_identity_heatmap.pdf",
        "Host Identity (%)",
    )
    save_heatmap(
        shared_log,
        FIGURES_DIR / "host_shared_clusters_heatmap.pdf",
        "Log10(Shared Clusters + 1)",
    )

    print(f"Hosts: {len(identity)}; unique nonself pairs: {len(records)}")
    print(f"Spearman Mantel R = {statistic:.6f}; P = {p_value:.6f}")


if __name__ == "__main__":
    main()
