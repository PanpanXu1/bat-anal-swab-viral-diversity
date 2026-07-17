#!/usr/bin/env python3
"""Generate the province-host viral spectrum heatmap for samples in this study."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from matplotlib.patches import FancyBboxPatch, Patch
import numpy as np
import pandas as pd


FONT_FAMILY = "Times New Roman"
WORKFLOW_ROOT = Path(__file__).resolve().parents[1]

for font_file in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
    font_path = Path("C:/Windows/Fonts") / font_file
    if font_path.exists():
        font_manager.fontManager.addfont(str(font_path))

VIRUS_TYPE_COLORS = {
    "DNA": "#f17465",
    "Reverse Transcribing Virus": "#1b9e77",
    "RNA": "#899bd1",
}
VIRUS_TYPE_ORDER = tuple(VIRUS_TYPE_COLORS)

BUBBLE_CMAP = LinearSegmentedColormap.from_list(
    "article_bubbles", ["#5b9d32", "#d2b94f", "#f2a05f", "#d9503d"]
)
HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "article_heatmap", ["#4f9c2a", "#a8cc65", "#f2e9a7", "#e59a72", "#d64a3a"]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=WORKFLOW_ROOT / "input")
    parser.add_argument("--output-dir", type=Path, default=WORKFLOW_ROOT / "output" / "figures")
    return parser.parse_args()


def load_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0).apply(pd.to_numeric, errors="coerce")


def row_zscore(matrix: pd.DataFrame) -> pd.DataFrame:
    means = matrix.mean(axis=1)
    std = matrix.std(axis=1, ddof=0).replace(0, 1)
    return matrix.sub(means, axis=0).div(std, axis=0)


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.pdf")
    plt.close(fig)


def italicize(labels) -> None:
    for label in labels:
        label.set_fontfamily(FONT_FAMILY)
        label.set_fontstyle("italic")


def style_legend(legend) -> None:
    if legend is None:
        return
    legend.get_title().set_fontfamily(FONT_FAMILY)
    legend.get_title().set_ha("left")
    legend._legend_box.align = "left"
    for text in legend.get_texts():
        text.set_fontfamily(FONT_FAMILY)


def virus_type_cmap() -> tuple[ListedColormap, BoundaryNorm]:
    cmap = ListedColormap([VIRUS_TYPE_COLORS[name] for name in VIRUS_TYPE_ORDER])
    norm = BoundaryNorm(np.arange(len(VIRUS_TYPE_ORDER) + 1) - 0.5, cmap.N)
    return cmap, norm


def virus_type_ids(types: pd.Series) -> np.ndarray:
    unknown = sorted(set(types.dropna()) - set(VIRUS_TYPE_ORDER))
    if unknown:
        raise ValueError(f"Unknown virus type(s): {', '.join(unknown)}")
    return types.map({name: idx for idx, name in enumerate(VIRUS_TYPE_ORDER)}).to_numpy()


def virus_type_handles(types: pd.Series | None = None) -> list[Patch]:
    observed = set(VIRUS_TYPE_ORDER) if types is None else set(types.dropna())
    labels = {
        "DNA": "DNA Virus",
        "Reverse Transcribing Virus": "Reverse Transcribing Virus",
        "RNA": "RNA Virus",
    }
    return [
        Patch(facecolor=VIRUS_TYPE_COLORS[name], label=labels[name])
        for name in VIRUS_TYPE_ORDER
        if name in observed
    ]


def draw_virus_type_bar(ax: plt.Axes, types: pd.Series) -> None:
    virus_type_ids(types)
    values = types.reset_index(drop=True)
    start = 0
    for row, virus_type in enumerate(values):
        if row > 0 and virus_type != values.iat[row - 1]:
            start = row
        next_type = values.iat[row + 1] if row + 1 < len(values) else None
        if next_type == virus_type:
            continue
        height = row - start + 1
        ax.add_patch(
            plt.Rectangle(
                (0, start - 0.5),
                1,
                height,
                facecolor=VIRUS_TYPE_COLORS[virus_type],
                edgecolor="white",
                linewidth=0.35,
            )
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(len(values) - 0.5, -0.5)
    ax.axis("off")


def plot_province_host_spectrum(data_dir: Path, output_dir: Path) -> None:
    matrix = load_matrix(data_dir / "viral_spectrum_by_province_host_rpm.csv")
    row_meta = pd.read_csv(data_dir / "viral_spectrum_by_province_host_row_metadata.csv")
    col_meta = pd.read_csv(data_dir / "viral_spectrum_by_province_host_column_metadata.csv")

    row_meta = row_meta.sort_values("original_order").reset_index(drop=True)
    matrix = matrix.loc[row_meta["virus_family"], col_meta["column"]]
    values = np.log2(matrix.to_numpy() + 1)
    sizes = np.where(values > 0, values * 2.0, 0)
    bubble_norm = plt.Normalize(vmin=5, vmax=25)

    fig = plt.figure(figsize=(22.0, 11.5))
    grid = fig.add_gridspec(
        2,
        4,
        width_ratios=[0.014, 0.69, 0.09, 0.206],
        height_ratios=[0.045, 0.955],
        wspace=0.01,
        hspace=0.015,
    )
    ax_top = fig.add_subplot(grid[0, 1])
    ax_left = fig.add_subplot(grid[1, 0])
    ax = fig.add_subplot(grid[1, 1])
    ax_legend = fig.add_subplot(grid[:, 3])

    ax_top.set_xlim(-0.5, matrix.shape[1] - 0.5)
    province_counts = col_meta.groupby("province", sort=False).size()
    start = 0
    for province, count in province_counts.items():
        ax_top.add_patch(
            plt.Rectangle(
                (start - 0.5, 0.02),
                count,
                0.90,
                fill=False,
                edgecolor="#173d2e",
                linewidth=1.0,
            )
        )
        ax_top.text(
            start + (count - 1) / 2,
            0.47,
            province,
            ha="center",
            va="center",
            fontsize=12,
            fontfamily=FONT_FAMILY,
        )
        start += count
    ax_top.set_ylim(0, 1)
    ax_top.axis("off")

    draw_virus_type_bar(ax_left, row_meta["virus_type"])

    x, y = np.meshgrid(np.arange(matrix.shape[1]), np.arange(matrix.shape[0]))
    ax.scatter(
        x.ravel(),
        y.ravel(),
        s=sizes.ravel(),
        c=values.ravel(),
        cmap=BUBBLE_CMAP,
        norm=bubble_norm,
        linewidths=0,
    )
    ax.set_xlim(-0.5, matrix.shape[1] - 0.5)
    ax.set_ylim(matrix.shape[0] - 0.5, -0.5)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(
        col_meta["host_species"], rotation=58, ha="right", fontsize=8.5
    )
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(row_meta["order_family_label"], fontsize=8.5)
    ax.yaxis.tick_right()
    ax.tick_params(axis="both", length=0)
    italicize(ax.get_xticklabels())
    italicize(ax.get_yticklabels())

    province_counts = col_meta.groupby("province", sort=False).size().tolist()
    for boundary in np.cumsum(province_counts)[:-1] - 0.5:
        ax.axvline(boundary, color="#173d2e", linewidth=1.25)
    type_counts = row_meta.groupby("virus_type", sort=False).size().tolist()
    for boundary in np.cumsum(type_counts)[:-1] - 0.5:
        ax.axhline(boundary, color="#173d2e", linewidth=1.25)
    for spine in ax.spines.values():
        spine.set_color("#173d2e")
        spine.set_linewidth(1.25)

    ax_legend.axis("off")
    type_legend = ax_legend.legend(
        handles=virus_type_handles(row_meta["virus_type"]),
        title="Virus Type",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.04, 0.99),
        fontsize=11,
        title_fontsize=12,
        labelspacing=0.8,
        handlelength=1.5,
        handleheight=1.5,
        alignment="left",
    )
    style_legend(type_legend)
    ax_legend.add_artist(type_legend)

    levels = [5, 10, 15, 20, 25]
    bubble_legend = ax_legend.legend(
        handles=[
            plt.scatter(
                [],
                [],
                s=level * 2.0,
                color=BUBBLE_CMAP(bubble_norm(level)),
                label=f"{level:.2f}",
            )
            for level in levels
        ],
        title=r"Log$_2$(RPM+1)",
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(0.04, 0.015),
        fontsize=10,
        title_fontsize=11,
        labelspacing=0.55,
    )
    style_legend(bubble_legend)

    fig.subplots_adjust(left=0.018, right=0.985, top=0.975, bottom=0.29)
    save_figure(fig, output_dir, "viral_spectrum_by_province_host")


def plot_article_clustered_panel(
    matrix_path: Path,
    metadata_path: Path,
    output_dir: Path,
    stem: str,
) -> None:
    matrix = load_matrix(matrix_path)
    if matrix.isna().any().any():
        raise ValueError(f"{matrix_path.name} contains missing values.")
    meta = pd.read_csv(metadata_path).set_index("virus_family")
    matrix = matrix.loc[meta.index]
    standardized = row_zscore(matrix)
    # Columns are kept in the supplied matrix order to avoid an additional clustering dependency.

    fig = plt.figure(figsize=(10.8, 12.5))
    grid = fig.add_gridspec(
        2,
        4,
        width_ratios=[0.026, 0.57, 0.15, 0.254],
        height_ratios=[0.11, 0.89],
        wspace=0.01,
        hspace=0.012,
    )
    ax_dendro = fig.add_subplot(grid[0, 1])
    ax_type = fig.add_subplot(grid[1, 0])
    ax_heat = fig.add_subplot(grid[1, 1])
    right = grid[:, 3].subgridspec(
        5, 3,
        width_ratios=[0.20, 0.12, 0.68],
        height_ratios=[0.08, 0.27, 0.09, 0.08, 0.48],
        hspace=0.08,
        wspace=0.05,
    )
    ax_cbar_title = fig.add_subplot(right[0, :])
    ax_cbar = fig.add_subplot(right[1, 0])
    ax_gap = fig.add_subplot(right[2, :])
    ax_legend_title = fig.add_subplot(right[3, :])
    ax_legend = fig.add_subplot(right[4, :])

    ax_dendro.axis("off")

    row_types = meta.loc[matrix.index, "virus_type"]
    draw_virus_type_bar(ax_type, row_types)

    norm = plt.Normalize(vmin=-3, vmax=4)
    for row in range(standardized.shape[0]):
        for col in range(standardized.shape[1]):
            patch = FancyBboxPatch(
                (col + 0.05, row + 0.08),
                0.90,
                0.84,
                boxstyle="round,pad=0.02,rounding_size=0.18",
                linewidth=0,
                facecolor=HEATMAP_CMAP(norm(standardized.iat[row, col])),
            )
            ax_heat.add_patch(patch)

    ax_heat.set_xlim(0, standardized.shape[1])
    ax_heat.set_ylim(standardized.shape[0], 0)
    ax_heat.set_xticks(np.arange(standardized.shape[1]) + 0.5)
    ax_heat.set_xticklabels(
        standardized.columns, rotation=48, ha="right", fontsize=12.0
    )
    ax_heat.set_yticks(np.arange(standardized.shape[0]) + 0.5)
    ax_heat.set_yticklabels(standardized.index, fontsize=12.0)
    ax_heat.yaxis.tick_right()
    ax_heat.tick_params(length=0)
    italicize(ax_heat.get_xticklabels())
    italicize(ax_heat.get_yticklabels())
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=HEATMAP_CMAP), cax=ax_cbar
    )
    colorbar.set_ticks([-3, -2, -1, 0, 1, 2, 3, 4])
    colorbar.ax.tick_params(labelsize=11)
    for tick in colorbar.ax.get_yticklabels():
        tick.set_fontfamily(FONT_FAMILY)

    for axis in (ax_cbar_title, ax_gap, ax_legend_title, ax_legend):
        axis.axis("off")
    ax_cbar_title.text(
        0,
        0.05,
        r"RPM (Log$_2$)-Normalized",
        fontsize=13,
        fontfamily=FONT_FAMILY,
        ha="left",
        va="bottom",
    )
    ax_legend_title.text(
        0,
        0.05,
        "Virus Type",
        fontsize=15,
        fontfamily=FONT_FAMILY,
        ha="left",
        va="bottom",
    )
    legend = ax_legend.legend(
        handles=virus_type_handles(row_types),
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(-0.02, 1.02),
        fontsize=13,
        handlelength=1.7,
        handleheight=1.7,
        labelspacing=1.0,
    )
    style_legend(legend)

    fig.subplots_adjust(top=0.975, bottom=0.14, left=0.025, right=0.985)
    save_figure(fig, output_dir, stem)


def main() -> None:
    args = parse_args()
    matplotlib.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": 11,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    plot_province_host_spectrum(args.data_dir, args.output_dir)


if __name__ == "__main__":
    main()
