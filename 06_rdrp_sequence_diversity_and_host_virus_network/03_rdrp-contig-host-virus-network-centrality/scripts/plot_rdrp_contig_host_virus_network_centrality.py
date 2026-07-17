#!/usr/bin/env python3
"""Create a portrait centrality leaderboard for manuscript panels."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib

matplotlib.use("Agg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd


REQUIRED_COLUMNS = {
    "Name",
    "Degree",
    "ClosenessCentrality",
    "NeighborhoodConnectivity",
    "Note",
}
WORKFLOW_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a ranked portrait centrality plot."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=WORKFLOW_ROOT / "input" / "rdrp_contig_host_virus_network_centrality.csv",
        help="Input CSV file (default: input/rdrp_contig_host_virus_network_centrality.csv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=WORKFLOW_ROOT / "output" / "figures" / "centrality_leaderboard_portrait_clean_legend.pdf",
        help="Output PDF file.",
    )
    parser.add_argument(
        "--top-virus",
        type=int,
        default=15,
        help="Number of virus families to show (default: 15).",
    )
    return parser.parse_args()


def configure_fonts() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_data(input_path: Path) -> pd.DataFrame:
    data = pd.read_csv(input_path)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
    return data.copy()


def marker_sizes(values: pd.Series, minimum: float, maximum: float) -> pd.Series:
    return 95 + 520 * (values - minimum) / (maximum - minimum)


def ranked_data(data: pd.DataFrame, top_virus: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    virus = (
        data[data["Note"] == "Virus"]
        .sort_values(["Degree", "ClosenessCentrality", "NeighborhoodConnectivity", "Name"], ascending=[False, False, False, True])
        .head(top_virus)
        .reset_index(drop=True)
    )
    host = (
        data[data["Note"] == "Host"]
        .sort_values(["Degree", "ClosenessCentrality", "NeighborhoodConnectivity", "Name"], ascending=[False, False, False, True])
        .reset_index(drop=True)
    )
    return virus, host


def draw_ranked_panel(
    axis: plt.Axes,
    data: pd.DataFrame,
    node_type: str,
    norm: mcolors.Normalize,
    closeness_min: float,
    closeness_max: float,
    x_max: int,
    label_size: int,
) -> plt.Collection:
    y_positions = range(len(data))
    sizes = marker_sizes(data["ClosenessCentrality"], closeness_min, closeness_max)

    axis.barh(
        y_positions,
        data["Degree"],
        height=0.58,
        color="#B8B8B8",
        edgecolor="none",
        zorder=1,
    )
    points = axis.scatter(
        data["Degree"],
        y_positions,
        s=sizes,
        c=data["NeighborhoodConnectivity"],
        cmap="Greens",
        norm=norm,
        edgecolors="black",
        linewidths=1.05,
        zorder=3,
    )

    axis.set_yticks(list(y_positions), labels=data["Name"])
    for label in axis.get_yticklabels():
        label.set_fontstyle("italic")
        label.set_fontsize(label_size)
    axis.invert_yaxis()
    axis.set_xlim(0, x_max)
    axis.set_xticks(range(0, x_max + 1, 5))
    axis.tick_params(axis="x", labelsize=13)
    axis.tick_params(axis="y", length=0, pad=4)
    axis.grid(axis="x", linestyle="--", linewidth=0.85, color="#D0D0D0")
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.2)
    axis.spines["bottom"].set_linewidth(1.2)

    return points


def draw_clean_legend_column(
    figure: plt.Figure,
    axis: plt.Axes,
    points: plt.Collection,
    closeness_min: float,
    closeness_max: float,
    norm: mcolors.Normalize,
) -> None:
    axis.axis("off")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)

    x0 = 0.10
    axis.text(
        x0,
        0.84,
        "Closeness Centrality",
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )
    values = pd.Series([closeness_max, (closeness_min + closeness_max) / 2, closeness_min])
    for y, size, label in zip(
        [0.74, 0.66, 0.59],
        marker_sizes(values, closeness_min, closeness_max),
        ["High", "Medium", "Low"],
    ):
        axis.scatter(
            [x0 + 0.13],
            [y],
            s=size,
            facecolors="white",
            edgecolors="black",
            linewidths=1.05,
        )
        axis.text(0.50, y, label, fontsize=12, va="center", ha="left")

    axis.text(
        x0,
        0.39,
        "Neighborhood\nConnectivity",
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
        linespacing=0.9,
    )
    cax = axis.inset_axes([x0 + 0.02, 0.10, 0.24, 0.24])
    colorbar = figure.colorbar(points, cax=cax, orientation="vertical")
    colorbar.set_ticks([norm.vmin, norm.vmax])
    colorbar.set_ticklabels([f"{norm.vmin:.1f}", f"{norm.vmax:.1f}"])
    colorbar.ax.tick_params(labelsize=9, length=2, pad=2)


def create_plot(data: pd.DataFrame, top_virus: int) -> plt.Figure:
    configure_fonts()
    virus, host = ranked_data(data, top_virus)
    closeness_min = data["ClosenessCentrality"].min()
    closeness_max = data["ClosenessCentrality"].max()
    norm = mcolors.Normalize(
        vmin=data["NeighborhoodConnectivity"].min(),
        vmax=data["NeighborhoodConnectivity"].max(),
    )

    figure = plt.figure(figsize=(7.2, 10.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, width_ratios=[1.0, 0.24], height_ratios=[0.56, 0.44])
    virus_axis = figure.add_subplot(grid[0, 0])
    host_axis = figure.add_subplot(grid[1, 0])
    legend_axis = figure.add_subplot(grid[:, 1])

    points = draw_ranked_panel(
        virus_axis,
        virus,
        "Virus",
        norm,
        closeness_min,
        closeness_max,
        x_max=16,
        label_size=13,
    )
    draw_ranked_panel(
        host_axis,
        host,
        "Host",
        norm,
        closeness_min,
        closeness_max,
        x_max=42,
        label_size=14,
    )

    virus_axis.set_xlabel("")
    host_axis.set_xlabel("Degree", fontsize=17)
    draw_clean_legend_column(figure, legend_axis, points, closeness_min, closeness_max, norm)
    return figure


def main() -> None:
    args = parse_args()
    data = load_data(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure = create_plot(data, args.top_virus)
    figure.savefig(args.output, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
