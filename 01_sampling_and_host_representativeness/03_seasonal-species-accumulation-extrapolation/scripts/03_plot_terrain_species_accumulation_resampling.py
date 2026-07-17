"""Plot terrain-stratified species accumulation curves using resampling."""

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


SEED = 42
N_STEPS = 100
N_ITERATIONS = 500

WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = WORKFLOW_ROOT / "input" / "species_abundance_by_terrain.csv"
INEXT_METRICS_FILE = WORKFLOW_ROOT / "input" / "terrain_inext_extrapolation_metrics.csv"
OUTPUT_FILE = WORKFLOW_ROOT / "output" / "figures" / "terrain_species_accumulation_extrapolation.pdf"

TERRAIN_ORDER = ["Plain", "Mountain", "Hill", "Mesa"]
SOURCE_COLUMNS = {
    "Plain": "Plain",
    "Mountain": "Mountainous",
    "Hill": "Hill",
    "Mesa": "Mesa",
}
COLORS = {
    "Plain": "#2878B5",
    "Mountain": "#4DAF4A",
    "Hill": "#F28E2B",
    "Mesa": "#D9534F",
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


def expected_richness(counts: np.ndarray, sample_size: int) -> float:
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts > 0]
    probabilities = counts / counts.sum()
    return float(np.sum(1.0 - np.power(1.0 - probabilities, sample_size)))


def simulate_curve(counts: np.ndarray, rng: np.random.Generator):
    counts = np.asarray(counts, dtype=int)
    counts = counts[counts > 0]
    total = int(counts.sum())
    probabilities = counts / total
    sample_sizes = np.unique(np.linspace(1, total, N_STEPS).round().astype(int))

    mean = []
    lower = []
    upper = []
    for n in sample_sizes:
        richness = np.count_nonzero(
            rng.multinomial(n, probabilities, size=N_ITERATIONS), axis=1
        )
        mean.append(np.mean(richness))
        lower.append(np.quantile(richness, 0.025))
        upper.append(np.quantile(richness, 0.975))

    q_n = expected_richness(counts, total)
    q_2n = expected_richness(counts, 2 * total)
    q_plus_100 = expected_richness(counts, total + 100)
    return {
        "sample_sizes": sample_sizes,
        "mean": np.asarray(mean),
        "lower": np.asarray(lower),
        "upper": np.asarray(upper),
        "gain_to_2n": q_2n - q_n,
        "slope_per_100": q_plus_100 - q_n,
    }


def read_abundance(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, index_col=0)
    return data.apply(pd.to_numeric, errors="raise").fillna(0).astype(int)


def read_inext_metrics(path: Path) -> dict[str, tuple[float, float]]:
    if not path.exists():
        raise FileNotFoundError(
            f"Required iNEXT metrics table does not exist: {path}. "
            "Run 02_calculate_terrain_inext_extrapolation_metrics.R first."
        )
    metrics = pd.read_csv(path)
    required = {"Terrain", "Gain_to_2N", "Slope_per_100"}
    if not required.issubset(metrics.columns):
        raise ValueError(f"iNEXT metrics table must contain {sorted(required)}")
    metrics["Gain_to_2N"] = pd.to_numeric(metrics["Gain_to_2N"], errors="raise")
    metrics["Slope_per_100"] = pd.to_numeric(metrics["Slope_per_100"], errors="raise")
    return {
        row["Terrain"]: (float(row["Gain_to_2N"]), float(row["Slope_per_100"]))
        for _, row in metrics.iterrows()
    }


def plot_terrain(curves: dict[str, dict], inext_metrics: dict[str, tuple[float, float]]) -> None:
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    y_max = max(float(curve["upper"].max()) for curve in curves.values()) * 1.08
    x_max = max(int(curve["sample_sizes"].max()) for curve in curves.values())

    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for terrain in TERRAIN_ORDER:
        curve = curves[terrain]
        color = COLORS[terrain]
        ax.fill_between(
            curve["sample_sizes"],
            curve["lower"],
            curve["upper"],
            color=color,
            alpha=0.20,
            linewidth=0,
        )
    for terrain in TERRAIN_ORDER:
        curve = curves[terrain]
        color = COLORS[terrain]
        gain_to_2n, slope_per_100 = inext_metrics[terrain]
        label = f"{terrain}   Gain (2N) = {gain_to_2n:.1f} | Slope (per100) = {slope_per_100:.1f}"
        ax.plot(curve["sample_sizes"], curve["mean"], color=color, linewidth=1.4, label=label)

    ax.set_xlim(0, x_max * 1.02)
    ax.set_ylim(0, y_max)
    ax.set_xlabel("Number of Samples")
    ax.set_ylabel("Species Richness")
    ax.grid(color="#E6E6E6", linewidth=0.8, alpha=0.9)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT_FILE, format="pdf")
    plt.close(fig)


def main() -> None:
    configure_fonts()
    terrain_data = read_abundance(INPUT_FILE)
    inext_metrics = read_inext_metrics(INEXT_METRICS_FILE)
    rng = np.random.default_rng(SEED)
    curves = {}
    for terrain in TERRAIN_ORDER:
        source_column = SOURCE_COLUMNS[terrain]
        curves[terrain] = simulate_curve(terrain_data[source_column].to_numpy(), rng)
    plot_terrain(curves, inext_metrics)
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
