from argparse import ArgumentParser
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib as mpl
mpl.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager


WORKFLOW_ROOT = Path(__file__).resolve().parents[1]

for font_file in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
    font_path = Path("C:/Windows/Fonts") / font_file
    if font_path.exists():
        font_manager.fontManager.addfont(str(font_path))

mpl.rcParams.update({
    "font.family": "Times New Roman",
    "font.serif": ["Times New Roman"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "savefig.facecolor": "white",
})


COLORS = {
    "ink": "#2C3033",
    "grid": "#E7E9EB",
    "green": "#5E8C61",
    "green_light": "#DCEBDD",
}


def parse_args():
    parser = ArgumentParser(
        description="Generate a PDF bar chart summarizing individual swabs per sequencing pool."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=WORKFLOW_ROOT / "input" / "pool_sample_metadata.xlsx",
        help="Input workbook path. Required columns: Number and pool_id.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=WORKFLOW_ROOT / "output" / "figures" / "pool_sample_count_distribution.pdf",
        help="Output PDF path.",
    )
    return parser.parse_args()

def load_pool_sizes(input_path):
    df = pd.read_excel(input_path)
    required = ["Number", "pool_id"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df.groupby("pool_id")["Number"].size().to_numpy()

def plot_pool_size_distribution(pool_sizes, output_path):
    n = len(pool_sizes)
    median = np.percentile(pool_sizes, 50)

    fig, ax = plt.subplots(figsize=(5.2, 3.35), constrained_layout=True)

    bin_edges = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 46]
    bin_labels = ["1", "2-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-35", "36-40", "41-46"]
    counts = []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        if lo == 0:
            counts.append(int((pool_sizes == 1).sum()))
        else:
            counts.append(int(((pool_sizes >= lo + 1) & (pool_sizes <= hi)).sum()))

    x = np.arange(len(bin_labels))
    bars = ax.bar(
        x,
        counts,
        width=0.72,
        color=COLORS["green_light"],
        edgecolor=COLORS["green"],
        linewidth=0.8,
    )
    for rect, count in zip(bars, counts):
        if count:
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                rect.get_height() + max(counts) * 0.025,
                str(count),
                ha="center",
                va="bottom",
                fontsize=7.5,
                color=COLORS["ink"],
            )

    ax.text(
        0.02,
        0.96,
        f"n = {n} pools\nmedian = {median:.0f} swabs",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=COLORS["ink"],
        fontsize=8,
    )

    ax.set_xlabel("Individual swabs per pool")
    ax.set_ylabel("Number of pools")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels)
    ax.set_ylim(0, max(counts) * 1.18)
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.6)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    pool_sizes = load_pool_sizes(args.input)
    plot_pool_size_distribution(pool_sizes, args.output)


if __name__ == "__main__":
    main()
