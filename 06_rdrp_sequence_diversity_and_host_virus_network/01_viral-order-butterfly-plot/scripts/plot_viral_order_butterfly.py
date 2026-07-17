"""Generate the NIPB-versus-ZOVER butterfly plot."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FuncFormatter


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "input" / "viral_order_log10_abundance_by_pool.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "figures" / "viral_order_butterfly_plot.pdf"
REQUIRED_COLUMNS = ("Order", "NIPB", "ZOVER")
NIPB_COLOR = "#FDBF6F"
ZOVER_COLOR = "#A695C5"
AXIS_LABEL_SIZE = 12
TICK_LABEL_SIZE = 11
LEGEND_FONT_SIZE = 11
PIE_TITLE_SIZE = 10
PIE_PERCENT_SIZE = 11
TIMES_NEW_ROMAN_PATHS = (
    Path("C:/Windows/Fonts/times.ttf"),
    Path("C:/Windows/Fonts/timesbd.ttf"),
    Path("/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a butterfly plot comparing NIPB and ZOVER."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input CSV file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output PDF file (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def load_data(input_path: Path) -> pd.DataFrame:
    """Load and validate the processed abundance data."""
    data = pd.read_csv(input_path)

    missing = [column for column in REQUIRED_COLUMNS if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")

    data = data.loc[:, REQUIRED_COLUMNS].copy()
    if data["Order"].isna().any():
        raise ValueError("The Order column contains missing values.")

    for column in ("NIPB", "ZOVER"):
        data[column] = pd.to_numeric(data[column], errors="raise")
        if data[column].isna().any():
            raise ValueError(f"The {column} column contains missing values.")
        if (data[column] < 0).any():
            raise ValueError(f"The {column} column contains negative values.")

    return data


def configure_fonts() -> None:
    """Configure all plot text to use Times New Roman."""
    for font_path in TIMES_NEW_ROMAN_PATHS:
        if font_path.exists():
            font_manager.fontManager.addfont(font_path)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def create_butterfly_plot(data: pd.DataFrame, output_path: Path) -> None:
    """Create and save the butterfly plot."""
    configure_fonts()

    x_positions = np.arange(len(data))
    figure, axis = plt.subplots(figsize=(14, 6.6))

    axis.bar(x_positions, -data["NIPB"], color=NIPB_COLOR, label="NIPB")
    axis.bar(x_positions, data["ZOVER"], color=ZOVER_COLOR, label="ZOVER")
    axis.axhline(0, color="black", linewidth=1)

    axis.set_xticks(x_positions)
    axis.set_xticklabels(
        data["Order"], rotation=45, ha="right", fontsize=TICK_LABEL_SIZE
    )
    for label in axis.get_xticklabels():
        label.set_fontstyle("italic")
    axis.yaxis.set_major_formatter(
        FuncFormatter(lambda value, position: f"{abs(value):g}")
    )
    axis.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)
    axis.set_ylabel("log10(Counts + 1)", fontsize=AXIS_LABEL_SIZE)
    axis.legend(frameon=False, loc="upper right", fontsize=LEGEND_FONT_SIZE)

    # Recover sequence counts from the supplied log10(counts + 1) values.
    sequence_counts = [
        np.sum(np.power(10, data["NIPB"]) - 1),
        np.sum(np.power(10, data["ZOVER"]) - 1),
    ]
    sequence_total = np.sum(sequence_counts)
    percentage_labels = [
        f"{count / sequence_total * 100:.0f}%" for count in sequence_counts
    ]
    pie_axis = axis.inset_axes([0.835, 0.62, 0.155, 0.25])
    pie_axis.pie(
        sequence_counts,
        colors=[NIPB_COLOR, ZOVER_COLOR],
        labels=percentage_labels,
        labeldistance=0.55,
        radius=0.78,
        startangle=90,
        counterclock=False,
        textprops={"fontsize": PIE_PERCENT_SIZE},
        wedgeprops={"edgecolor": "white", "linewidth": 0.8},
    )
    pie_axis.set_title("Sequence proportion", fontsize=PIE_TITLE_SIZE, pad=2)
    pie_axis.set_aspect("equal")
    pie_axis.set_frame_on(True)
    pie_axis.patch.set_facecolor("white")
    pie_axis.patch.set_edgecolor("#D0D0D0")
    pie_axis.patch.set_linewidth(0.8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(output_path, format="pdf")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    data = load_data(args.input)
    create_butterfly_plot(data, args.output)
    print(f"Saved butterfly plot to: {args.output.resolve()}")


if __name__ == "__main__":
    main()
