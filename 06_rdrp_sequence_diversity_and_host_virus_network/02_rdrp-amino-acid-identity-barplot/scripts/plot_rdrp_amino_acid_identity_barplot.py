import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())
import csv

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "input" / "virus_order_amino_acid_identity_counts.csv"
OUTPUT_PATH = ROOT / "output" / "figures" / "rdrp_amino_acid_identity_barplot.pdf"

COLORS = ["#FCE6C6", "#FFB3DD", "#B7E4FF"]
LEGEND_LABELS = [
    "Amino acid identity <50%",
    "50% \u2264 Amino acid identity <90%",
    "Amino acid identity \u226590%",
]

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.size": 16,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def load_data(path):
    labels = []
    values = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if len(header) < 4 or header[0] != "Virus order":
            raise ValueError("Input table must contain Virus order and three amino-acid identity count columns")
        for row in reader:
            labels.append(row[0])
            values.append(
                [
                    float(row[1]),
                    float(row[2]),
                    float(row[3]),
                ]
            )
    return labels, np.array(values)


def main():
    labels, values = load_data(DATA_PATH)
    x = np.arange(len(labels))
    width = 0.24

    fig, ax = plt.subplots(figsize=(772 / 72, 591 / 72))

    for idx, (label, color, offset) in enumerate(zip(LEGEND_LABELS, COLORS, (-width, 0, width))):
        ax.bar(
            x + offset,
            values[:, idx],
            width=width,
            label=label,
            color=color,
            edgecolor="black",
            linewidth=0.55,
        )

    ax.set_ylabel("Number of sequences", fontsize=20)
    ax.set_xlabel("Virus order", fontsize=20)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=15, fontstyle="italic")
    ax.tick_params(axis="y", labelsize=17)
    ax.set_ylim(0, max(values.max() * 1.12, 10))
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.75)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.legend(
        loc="upper left",
        frameon=False,
        fontsize=16,
        handlelength=1.3,
        borderaxespad=0.4,
    )

    fig.tight_layout(pad=1.2)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PATH, format="pdf")
    plt.close(fig)
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
