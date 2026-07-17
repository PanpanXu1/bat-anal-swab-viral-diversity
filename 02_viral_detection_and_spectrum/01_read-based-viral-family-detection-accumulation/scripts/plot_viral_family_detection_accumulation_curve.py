from pathlib import Path
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap


PROJECT_DIR = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_DIR / "input" / "viral_family_read_counts_by_sample.csv"
FIGURES_DIR = PROJECT_DIR / "output" / "figures"

FIGURE_FILE = FIGURES_DIR / "viral_family_detection_accumulation_curve_0_20000_reads.pdf"

MAX_ASSIGNED_READS_TO_PLOT = 20_000


def load_viral_family_read_counts(input_file: Path) -> pd.DataFrame:
    df = pd.read_csv(input_file)
    required_columns = {"sample_id", "viral_family_tax_id", "assigned_reads"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    df = df.copy()
    df["assigned_reads"] = pd.to_numeric(df["assigned_reads"], errors="raise")
    return df


def calculate_detection_accumulation(df: pd.DataFrame) -> pd.DataFrame:
    grouped_counts = (
        df.groupby(["sample_id", "assigned_reads"], as_index=False)
        .agg(viral_families_at_read_count=("viral_family_tax_id", "nunique"))
        .sort_values(["sample_id", "assigned_reads"])
    )
    grouped_counts["cumulative_viral_families"] = grouped_counts.groupby("sample_id")[
        "viral_families_at_read_count"
    ].cumsum()
    return grouped_counts


def plot_detection_accumulation(curve_df: pd.DataFrame, output_file: Path) -> None:
    for font_file in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
        font_path = Path("C:/Windows/Fonts") / font_file
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
        }
    )

    colors = LinearSegmentedColormap.from_list(
        "warm_detection_tone", ["#f6c568", "#f2b84b", "#edae49"]
    )

    sample_ids = sorted(curve_df["sample_id"].unique())
    fig, ax = plt.subplots(figsize=(10, 5.5))

    for index, sample_id in enumerate(sample_ids):
        sample_curve = curve_df[
            (curve_df["sample_id"] == sample_id)
            & (curve_df["assigned_reads"] <= MAX_ASSIGNED_READS_TO_PLOT)
        ]
        if sample_curve.empty:
            continue
        ax.plot(
            sample_curve["assigned_reads"],
            sample_curve["cumulative_viral_families"],
            color=colors(index / max(len(sample_ids) - 1, 1)),
            linewidth=0.35,
            alpha=0.72,
        )

    ax.set_xlim(0, MAX_ASSIGNED_READS_TO_PLOT)
    ax.set_ylim(0, 24)
    ax.set_xticks(range(0, MAX_ASSIGNED_READS_TO_PLOT + 1, 1000))
    ax.set_yticks(list(range(0, 11, 2)) + list(range(12, 25, 2)))
    ax.set_xlabel("Assigned reads per viral family (0-20,000)")
    ax.set_ylabel("Cumulative detected viral families")
    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", labelsize=7)

    fig.tight_layout()
    fig.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df = load_viral_family_read_counts(INPUT_FILE)
    curve_df = calculate_detection_accumulation(df)
    plot_detection_accumulation(curve_df, FIGURE_FILE)

    print(f"Saved PDF figure: {FIGURE_FILE}")


if __name__ == "__main__":
    main()
