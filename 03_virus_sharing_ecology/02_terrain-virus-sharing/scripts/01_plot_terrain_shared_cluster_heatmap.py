"""Create a terrain-ordered shared viral-cluster heatmap."""

import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = ROOT / "input" / "terrain_viral_cluster_membership.csv"
FIGURES_DIR = ROOT / "output" / "figures"
HEATMAP_FILE = FIGURES_DIR / "terrain_shared_cluster_heatmap.pdf"

TERRAIN_ORDER = ["Mountain", "Plain", "Hill", "Mesa"]
TERRAIN_COLORS = {
    "Mountain": "#66c2a5",
    "Plain": "#fc8d62",
    "Hill": "#8da0cb",
    "Mesa": "#e78ac3",
}


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


def load_data() -> pd.DataFrame:
    data = pd.read_csv(DATA_FILE)
    required = {"Sequence_ID", "Cluster_ID", "SampleGroup", "Terrain"}
    if not required.issubset(data.columns):
        raise ValueError(f"Input must contain {sorted(required)}")
    data = data.dropna(subset=list(required)).drop_duplicates()
    mapping_counts = data.groupby("SampleGroup")["Terrain"].nunique()
    if (mapping_counts != 1).any():
        raise ValueError("Each SampleGroup must map to exactly one Terrain")
    return data


def calculate_shared_matrix(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    sample_terrain = data.drop_duplicates("SampleGroup").set_index("SampleGroup")["Terrain"]
    unknown = sorted(set(sample_terrain) - set(TERRAIN_ORDER))
    if unknown:
        raise ValueError(f"Unexpected terrain categories: {unknown}")

    ordered_samples = sample_terrain.sort_values().index.tolist()

    presence = pd.crosstab(data["SampleGroup"], data["Cluster_ID"]).gt(0).astype(np.int16)
    shared = presence.dot(presence.T).reindex(index=ordered_samples, columns=ordered_samples).copy()
    for sample in ordered_samples:
        shared.loc[sample, sample] = 0
    return np.log10(shared + 1), sample_terrain.reindex(ordered_samples)


def draw_terrain_bar(
    ax: plt.Axes,
    sample_terrain: pd.Series,
    orientation: str,
) -> None:
    """Draw continuous terrain blocks without per-sample rendering seams."""
    terrain_values = sample_terrain.to_list()
    n_samples = len(terrain_values)
    start = 0
    for idx in range(1, n_samples + 1):
        if idx == n_samples or terrain_values[idx] != terrain_values[start]:
            terrain = terrain_values[start]
            color = TERRAIN_COLORS[terrain]
            if orientation == "vertical":
                ax.add_patch(Rectangle((0, start), 1, idx - start, color=color, linewidth=0))
            else:
                ax.add_patch(Rectangle((start, 0), idx - start, 1, color=color, linewidth=0))
            start = idx

    ax.set_xlim(0, 1 if orientation == "vertical" else n_samples)
    ax.set_ylim(n_samples if orientation == "vertical" else 0, 0 if orientation == "vertical" else 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main() -> None:
    configure_fonts()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data = load_data()
    shared_log, sample_terrain = calculate_shared_matrix(data)

    heatmap_cmap = "coolwarm"
    heatmap_values = shared_log.to_numpy()
    heatmap_vmax = float(heatmap_values.max())
    heatmap_center = float(heatmap_values.mean())

    fig = plt.figure(figsize=(34, 34))
    heatmap_ax = fig.add_axes([0.085, 0.085, 0.76, 0.76])
    left_terrain_ax = fig.add_axes([0.065, 0.085, 0.012, 0.76])
    bottom_terrain_ax = fig.add_axes([0.085, 0.065, 0.76, 0.012])
    legend_ax = fig.add_axes([0.865, 0.745, 0.12, 0.16])
    cbar_ax = fig.add_axes([0.865, 0.085, 0.028, 0.66])

    draw_terrain_bar(left_terrain_ax, sample_terrain, "vertical")

    image = heatmap_ax.imshow(shared_log.to_numpy(), cmap=heatmap_cmap, vmin=0, vmax=heatmap_vmax, aspect="equal")
    fig.colorbar(image, cax=cbar_ax)
    heatmap_ax.set_xlabel("")
    heatmap_ax.set_ylabel("")
    heatmap_ax.set_xticklabels([])
    heatmap_ax.set_yticklabels([])
    heatmap_ax.tick_params(axis="both", length=0)
    heatmap_ax.set_aspect("equal")
    cbar_ax.tick_params(labelsize=32, width=1.6, length=8)

    draw_terrain_bar(bottom_terrain_ax, sample_terrain, "horizontal")

    handles = [Patch(color=TERRAIN_COLORS[name], label=name) for name in TERRAIN_ORDER]
    legend_ax.legend(
        handles=handles,
        title="Terrain Types",
        loc="lower left",
        frameon=False,
        fontsize=34,
        title_fontsize=36,
        handlelength=1.6,
        handleheight=1.0,
        borderaxespad=0,
        labelspacing=0.5,
    )
    legend_ax.axis("off")

    try:
        fig.savefig(HEATMAP_FILE, bbox_inches="tight")
    except PermissionError as exc:
        raise PermissionError(
            f"Cannot write {HEATMAP_FILE}. Close the existing PDF if it is open and rerun."
        ) from exc
    plt.close(fig)
    print(f"Created heatmap for {len(shared_log)} sampling pools")


if __name__ == "__main__":
    main()
