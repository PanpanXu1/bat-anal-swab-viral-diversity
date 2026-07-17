"""Draw a Spearman-correlation matrix with Shannon-environment Mantel links."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SPEARMAN_MATRIX = (
    ROOT
    / "output"
    / "tables"
    / "environmental_spearman"
    / "environmental_spearman_correlation_matrix.csv"
)
DEFAULT_SPEARMAN_TESTS = (
    ROOT
    / "output"
    / "tables"
    / "environmental_spearman"
    / "environmental_spearman_pairwise_tests.csv"
)
DEFAULT_MANTEL_TESTS = (
    ROOT / "output" / "tables" / "shannon_environment_mantel" / "shannon_environment_mantel_tests.csv"
)
DEFAULT_OUTPUT_DIR = ROOT / "output" / "figures" / "environmental_correlation_mantel_plot"

VARIABLE_ORDER = [
    "Bio19",
    "Bio18",
    "Bio17",
    "Bio16",
    "Bio15",
    "Bio14",
    "Bio13",
    "Bio12",
    "Bio11",
    "Bio10",
    "Bio9",
    "Bio8",
    "Bio7",
    "Bio6",
    "Bio5",
    "Bio4",
    "Bio3",
    "Bio2",
    "Bio1",
    "PSD",
    "GR",
    "GLH",
    "GMR",
    "HFT",
    "GDP",
    "FVC",
    "CHEQ",
    "NDVI",
    "DEM",
]

CORRELATION_BOUNDS = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
CORRELATION_COLORS = [
    "#2367b7",
    "#5c7fc2",
    "#8e97ca",
    "#bcb0d7",
    "#ddd6ea",
    "#f1d5c9",
    "#e7a7a4",
    "#de7d78",
    "#ee886f",
    "#d8583c",
]
# Matrix cells are 1 x 1 in data coordinates; square side length equals |rho|
# so the size legend and matrix use the same physical scale.
MAX_CORRELATION_SQUARE_SIZE = 1.0
ZERO_CORRELATION_LEGEND_SIZE = 0.035


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a triangular environmental Spearman-correlation plot with "
            "Shannon-environment Mantel-test links."
        )
    )
    parser.add_argument("--spearman-matrix", type=Path, default=DEFAULT_SPEARMAN_MATRIX)
    parser.add_argument("--spearman-tests", type=Path, default=DEFAULT_SPEARMAN_TESTS)
    parser.add_argument("--mantel-tests", type=Path, default=DEFAULT_MANTEL_TESTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--basename", default="environmental_correlation_mantel_plot")
    parser.add_argument(
        "--mantel-p-threshold",
        type=float,
        default=0.05,
        help="Only draw Shannon-environment Mantel links at or below this p-value.",
    )
    return parser.parse_args()


def load_inputs(
    spearman_matrix_path: Path, spearman_tests_path: Path, mantel_tests_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_csv(spearman_matrix_path, index_col=0)
    spearman_tests = pd.read_csv(spearman_tests_path)
    mantel_tests = pd.read_csv(mantel_tests_path)

    missing = set(VARIABLE_ORDER) - set(matrix.index) - set(matrix.columns)
    if missing:
        raise ValueError(f"Spearman matrix is missing variables: {sorted(missing)}")
    matrix = matrix.loc[VARIABLE_ORDER, VARIABLE_ORDER]

    required_spearman = {"Factor1", "Factor2", "p_value"}
    if not required_spearman.issubset(spearman_tests.columns):
        raise ValueError(f"Spearman tests must contain columns: {sorted(required_spearman)}")

    required_mantel = {"Environmental_variable", "Mantel_r", "p_value"}
    if not required_mantel.issubset(mantel_tests.columns):
        raise ValueError(f"Mantel tests must contain columns: {sorted(required_mantel)}")

    return matrix, spearman_tests, mantel_tests


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif", "serif"],
            "pdf.fonttype": 42,
            "font.size": 10,
            "axes.linewidth": 0.6,
            "axes.spines.left": False,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.spines.bottom": False,
        }
    )


def add_correlation_matrix(
    ax: plt.Axes,
    matrix: pd.DataFrame,
    cmap: ListedColormap,
    norm: BoundaryNorm,
) -> None:
    n_variables = len(VARIABLE_ORDER)
    ax.plot(
        [0.5, n_variables - 0.5],
        [-0.5, -0.5],
        color="#111111",
        linewidth=0.55,
        zorder=0,
    )
    ax.plot(
        [n_variables - 0.5, n_variables - 0.5],
        [-0.5, n_variables - 1.5],
        color="#111111",
        linewidth=0.55,
        zorder=0,
    )
    for position in np.arange(0.5, n_variables - 0.5, 1):
        ax.plot(
            [position, n_variables - 0.5],
            [position, position],
            color="#111111",
            linewidth=0.55,
            zorder=0,
        )
        ax.plot(
            [position, position],
            [-0.5, position],
            color="#111111",
            linewidth=0.55,
            zorder=0,
        )

    for row_index, row_variable in enumerate(VARIABLE_ORDER):
        for column_index, column_variable in enumerate(VARIABLE_ORDER):
            if column_index <= row_index:
                continue

            rho = float(matrix.loc[row_variable, column_variable])
            size = correlation_square_size(abs(rho))

            lower_left = (column_index - size / 2, row_index - size / 2)
            square = Rectangle(
                lower_left,
                size,
                size,
                facecolor=cmap(norm(rho)),
                edgecolor="#6d6d6d",
                linewidth=0.35,
            )
            ax.add_patch(square)

    diagonal_color = "#de6b3d"
    ax.scatter(
        np.arange(n_variables),
        np.arange(n_variables),
        s=30,
        color=diagonal_color,
        edgecolors="white",
        linewidths=0.45,
        zorder=5,
    )

    for index, variable in enumerate(VARIABLE_ORDER):
        ax.text(
            index,
            -0.72,
            variable,
            rotation=52,
            rotation_mode="anchor",
            ha="left",
            va="bottom",
            fontsize=9,
        )
        ax.text(n_variables + 0.12, index, variable, ha="left", va="center", fontsize=9)


def correlation_square_size(abs_rho: float) -> float:
    return MAX_CORRELATION_SQUARE_SIZE * abs_rho


def mantel_line_width(abs_r: float) -> float:
    return 2.2


def add_mantel_links(
    ax: plt.Axes, mantel_tests: pd.DataFrame, p_threshold: float
) -> None:
    n_variables = len(VARIABLE_ORDER)
    anchor = (VARIABLE_ORDER.index("Bio17"), VARIABLE_ORDER.index("FVC"))
    ax.scatter(
        [anchor[0]],
        [anchor[1]],
        s=32,
        color="#de6b3d",
        edgecolors="white",
        linewidths=0.45,
        zorder=6,
    )
    ax.text(anchor[0] - 1.45, anchor[1] + 1.15, "Shannon Index", ha="left", va="center", fontsize=13)

    mantel_by_variable = mantel_tests.set_index("Environmental_variable")
    for variable in VARIABLE_ORDER:
        if variable not in mantel_by_variable.index:
            continue
        row = mantel_by_variable.loc[variable]
        statistic = float(row["Mantel_r"])
        p_value = float(row["p_value"])
        if p_value > p_threshold:
            continue

        target_index = VARIABLE_ORDER.index(variable)
        color = "#f6b9a5" if statistic > 0 else "#9d98ca"
        ax.plot(
            [anchor[0], target_index],
            [anchor[1], target_index],
            color=color,
            linewidth=mantel_line_width(abs(statistic)),
            alpha=0.65,
            solid_capstyle="round",
            zorder=1,
        )
        midpoint_x = (anchor[0] + target_index) / 2
        midpoint_y = (anchor[1] + target_index) / 2
        ax.text(
            midpoint_x + 0.12,
            midpoint_y - 0.12,
            f"{statistic:.2f}",
            ha="center",
            va="center",
            fontsize=8.5,
            color="#333333",
            zorder=7,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.8),
        )


def add_legends(
    fig: plt.Figure,
    ax: plt.Axes,
    cmap: ListedColormap,
    norm: BoundaryNorm,
) -> None:
    cax = fig.add_axes([0.690, 0.55, 0.032, 0.24])
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    colorbar = fig.colorbar(sm, cax=cax, boundaries=CORRELATION_BOUNDS)
    colorbar.set_ticks([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
    colorbar.set_ticklabels(
        [
            "-1 - -0.8",
            "-0.8 - -0.6",
            "-0.6 - -0.4",
            "-0.4 - -0.2",
            "-0.2 - 0",
            "0 - 0.2",
            "0.2 - 0.4",
            "0.4 - 0.6",
            "0.6 - 0.8",
            "0.8 - 1",
        ]
    )
    colorbar.ax.tick_params(labelsize=9, length=0)
    colorbar.outline.set_visible(False)
    fig.text(0.690, 0.80, "Environment factors\ncorrelation", ha="left", va="bottom", fontsize=10)

    fig.text(
        0.690,
        0.500,
        "|Environment factors\ncorrelation|",
        ha="left",
        va="bottom",
        fontsize=10,
    )
    ax_position = ax.get_position()
    x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
    y_range = ax.get_ylim()[0] - ax.get_ylim()[1]
    data_unit_as_figure = min(ax_position.width / x_range, ax_position.height / y_range)
    legend_center_x = 0.705
    for center_y, rho, label in [(0.462, 1.0, "1"), (0.418, 0.5, "0.5"), (0.374, 0.0, "0")]:
        square_size = correlation_square_size(rho) * data_unit_as_figure
        if rho == 0:
            square_size = ZERO_CORRELATION_LEGEND_SIZE * data_unit_as_figure
        fig.patches.append(
            Rectangle(
                (legend_center_x - square_size / 2, center_y - square_size / 2),
                square_size,
                square_size,
                transform=fig.transFigure,
                facecolor="#2b2f35",
                edgecolor="#2b2f35",
                linewidth=0.3,
                clip_on=False,
            )
        )
        fig.text(0.735, center_y, label, ha="left", va="center", fontsize=10)

    line_handles = [
        Line2D([0], [0], color="#9d98ca", lw=5, alpha=0.65, label="Mantel r < 0"),
        Line2D([0], [0], color="#f6b9a5", lw=5, alpha=0.65, label="Mantel r > 0"),
    ]
    line_legend = fig.legend(
        handles=line_handles,
        loc="lower left",
        bbox_to_anchor=(0.690, 0.250),
        frameon=True,
        framealpha=1,
        borderpad=0.7,
        handlelength=1.6,
        fontsize=10,
    )
    line_legend.get_frame().set_edgecolor("#bdbdbd")
    line_legend.get_frame().set_linewidth(0.8)


def draw_figure(
    matrix: pd.DataFrame,
    mantel_tests: pd.DataFrame,
    output_dir: Path,
    basename: str,
    p_threshold: float,
) -> None:
    configure_matplotlib()
    cmap = ListedColormap(CORRELATION_COLORS)
    norm = BoundaryNorm(CORRELATION_BOUNDS, cmap.N, clip=True)

    n_variables = len(VARIABLE_ORDER)
    fig, ax = plt.subplots(figsize=(8.7, 8.4))
    fig.subplots_adjust(left=0.04, right=0.68, top=0.93, bottom=0.08)

    add_mantel_links(ax, mantel_tests, p_threshold)
    add_correlation_matrix(ax, matrix, cmap, norm)

    ax.set_xlim(-0.25, n_variables + 2.5)
    ax.set_ylim(n_variables + 2.2, -1.2)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.canvas.draw()
    add_legends(fig, ax, cmap, norm)

    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    matrix, spearman_tests, mantel_tests = load_inputs(
        args.spearman_matrix, args.spearman_tests, args.mantel_tests
    )
    draw_figure(
        matrix=matrix,
        mantel_tests=mantel_tests,
        output_dir=args.output_dir,
        basename=args.basename,
        p_threshold=args.mantel_p_threshold,
    )
    print(f"Figure written to: {args.output_dir}")


if __name__ == "__main__":
    main()
