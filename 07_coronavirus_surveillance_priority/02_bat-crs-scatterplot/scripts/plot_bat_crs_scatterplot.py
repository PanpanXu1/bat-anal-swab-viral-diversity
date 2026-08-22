import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from mpl_toolkits.mplot3d import proj3d


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "input" / "bat_crs_scatterplot_input.csv"
OUTPUT_PATH = ROOT / "output" / "figures" / "bat_crs_scatterplot.pdf"

LABEL_FONT_SIZE = 6.4
LABEL_HEIGHT = 0.022
LABEL_CHAR_WIDTH = 0.0062
POINT_CLEARANCE = 0.018
LINE_THRESHOLD = 0.024

OFFSET_CANDIDATES = [
    (0.032, 0.020),
    (-0.032, 0.020),
    (0.032, -0.020),
    (-0.032, -0.020),
    (0.022, 0.050),
    (-0.022, 0.050),
    (0.000, 0.060),
    (0.085, 0.050),
    (-0.085, 0.050),
    (0.115, 0.040),
    (-0.115, 0.040),
    (0.048, 0.000),
    (-0.048, 0.000),
    (0.000, 0.035),
    (0.056, 0.026),
    (-0.056, 0.026),
    (0.056, -0.026),
    (-0.056, -0.026),
    (0.074, 0.000),
    (-0.074, 0.000),
]


def configure_matplotlib():
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
            "font.size": 8,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
        }
    )


def load_data(path):
    required = [
        "Host",
        "Viral Diversity_norm",
        "Zoonotic Species Count_norm",
        "ZSN_corrected_norm",
        "VS_log_norm",
        "CRS",
    ]
    df = pd.read_csv(path)
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    for column in [
        "Viral Diversity_norm",
        "Zoonotic Species Count_norm",
        "ZSN_corrected_norm",
        "VS_log_norm",
        "CRS",
    ]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df[required].copy()


def label_box(center, text):
    width = max(0.060, len(text) * LABEL_CHAR_WIDTH)
    height = LABEL_HEIGHT
    x_center, y_center = center
    return (
        x_center - width / 2,
        y_center - height / 2,
        x_center + width / 2,
        y_center + height / 2,
    )


def boxes_overlap(box_a, box_b, padding=0.003):
    return not (
        box_a[2] + padding < box_b[0]
        or box_b[2] + padding < box_a[0]
        or box_a[3] + padding < box_b[1]
        or box_b[3] + padding < box_a[1]
    )


def point_in_box(point, box, padding=POINT_CLEARANCE):
    x_point, y_point = point
    return (
        box[0] - padding <= x_point <= box[2] + padding
        and box[1] - padding <= y_point <= box[3] + padding
    )


def ccw(point_a, point_b, point_c):
    return (
        (point_c[1] - point_a[1]) * (point_b[0] - point_a[0])
        > (point_b[1] - point_a[1]) * (point_c[0] - point_a[0])
    )


def segments_intersect(seg_a, seg_b):
    a_start, a_end = seg_a
    b_start, b_end = seg_b
    shared_endpoint = (
        np.allclose(a_start, b_start)
        or np.allclose(a_start, b_end)
        or np.allclose(a_end, b_start)
        or np.allclose(a_end, b_end)
    )
    if shared_endpoint:
        return False
    return (
        ccw(a_start, b_start, b_end) != ccw(a_end, b_start, b_end)
        and ccw(a_start, a_end, b_start) != ccw(a_start, a_end, b_end)
    )


def line_start_for_box(point, box):
    x_point, y_point = point
    x_center = (box[0] + box[2]) / 2
    y_center = (box[1] + box[3]) / 2
    dx = x_center - x_point
    dy = y_center - y_point

    if abs(dx) >= abs(dy):
        x_start = box[0] if dx >= 0 else box[2]
        y_start = np.clip(y_point, box[1], box[3])
    else:
        y_start = box[1] if dy >= 0 else box[3]
        x_start = np.clip(x_point, box[0], box[2])
    return (float(x_start), float(y_start))


def candidate_label_positions(point, text):
    x_point, y_point = point
    candidates = []
    for dx, dy in OFFSET_CANDIDATES:
        center = (
            float(np.clip(x_point + dx, 0.09, 0.84)),
            float(np.clip(y_point + dy, 0.31, 0.91)),
        )
        box = label_box(center, text)
        segment_start = line_start_for_box(point, box)
        distance = float(np.hypot(center[0] - x_point, center[1] - y_point))
        candidates.append(
            {
                "center": center,
                "box": box,
                "line": (segment_start, point),
                "distance": distance,
            }
        )
    return candidates


def update_label_geometry(label, center):
    center = (
        float(np.clip(center[0], 0.09, 0.84)),
        float(np.clip(center[1], 0.31, 0.91)),
    )
    box = label_box(center, label["host"])
    line = (line_start_for_box(label["point"], box), label["point"])
    distance = float(np.hypot(center[0] - label["point"][0], center[1] - label["point"][1]))
    label["center"] = center
    label["box"] = box
    label["line"] = line
    label["distance"] = distance


def relax_label_overlaps(selected):
    for _ in range(10):
        changed = False
        for idx in range(len(selected)):
            for jdx in range(idx + 1, len(selected)):
                first = selected[idx]
                second = selected[jdx]
                if not boxes_overlap(first["box"], second["box"], padding=0.006):
                    continue

                if second["point"][0] >= first["point"][0]:
                    direction = 1
                else:
                    direction = -1
                overlap_width = min(first["box"][2], second["box"][2]) - max(
                    first["box"][0], second["box"][0]
                )
                shift_x = direction * (max(overlap_width, 0.0) + 0.035)
                new_center = (
                    second["center"][0] + shift_x,
                    second["center"][1],
                )
                if new_center[0] < 0.10 or new_center[0] > 0.83:
                    new_center = (
                        second["center"][0] - shift_x,
                        second["center"][1] + 0.035,
                    )
                update_label_geometry(second, new_center)
                changed = True
        if not changed:
            break
    return selected


def layout_labels(hosts, projected_points):
    selected = []
    point_list = list(projected_points)

    # Place more important/high-CRS hosts first, preserving input rank order.
    for host, point in zip(hosts, point_list):
        best_candidate = None
        best_score = float("inf")
        for candidate in candidate_label_positions(point, host):
            box = candidate["box"]
            score = candidate["distance"] * 100

            if box[0] < 0.05 or box[2] > 0.90 or box[1] < 0.28 or box[3] > 0.94:
                score += 1000
            if box[1] < 0.33:
                score += 250
            if box[2] > 0.86 and box[1] < 0.35:
                score += 350
            if point_in_box(point, box, padding=0.010):
                score += 500
            for other_point in point_list:
                if other_point != point and point_in_box(other_point, box):
                    score += 300
            for placed in selected:
                if boxes_overlap(box, placed["box"]):
                    score += 700
                if candidate["distance"] > LINE_THRESHOLD and placed["distance"] > LINE_THRESHOLD:
                    if segments_intersect(candidate["line"], placed["line"]):
                        score += 650

            if score < best_score:
                best_score = score
                best_candidate = candidate

        selected.append({"host": host, "point": point, **best_candidate})

    return relax_label_overlaps(selected)


def clamp_axes_position(position):
    return (
        float(np.clip(position[0], 0.08, 0.86)),
        float(np.clip(position[1], 0.30, 0.92)),
    )


def display_delta_to_axes(ax, delta_x, delta_y):
    origin = ax.transAxes.transform((0, 0))
    target = origin + np.array([delta_x, delta_y])
    return ax.transAxes.inverted().transform(target) - ax.transAxes.inverted().transform(origin)


def adjust_text_artists(fig, ax, text_records):
    for _ in range(120):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        bboxes = [
            record["artist"].get_window_extent(renderer).expanded(1.08, 1.20)
            for record in text_records
        ]
        displacements = [np.zeros(2, dtype=float) for _ in text_records]

        for idx in range(len(text_records)):
            for jdx in range(idx + 1, len(text_records)):
                box_a = bboxes[idx]
                box_b = bboxes[jdx]
                if not box_a.overlaps(box_b):
                    continue
                overlap_x = min(box_a.x1, box_b.x1) - max(box_a.x0, box_b.x0)
                overlap_y = min(box_a.y1, box_b.y1) - max(box_a.y0, box_b.y0)
                if overlap_x <= 0 or overlap_y <= 0:
                    continue

                center_a = np.array([(box_a.x0 + box_a.x1) / 2, (box_a.y0 + box_a.y1) / 2])
                center_b = np.array([(box_b.x0 + box_b.x1) / 2, (box_b.y0 + box_b.y1) / 2])
                direction = center_b - center_a
                if np.allclose(direction, 0):
                    direction = np.array([1.0, 0.0])
                direction = direction / np.linalg.norm(direction)

                if overlap_x < overlap_y:
                    push = np.array([np.sign(direction[0]) * (overlap_x / 2 + 6), 0.0])
                else:
                    push = np.array([0.0, np.sign(direction[1]) * (overlap_y / 2 + 5)])
                displacements[idx] -= push
                displacements[jdx] += push

        for idx, record in enumerate(text_records):
            box = bboxes[idx]
            point_display = ax.transAxes.transform(record["point"])
            if box.expanded(1.08, 1.35).contains(*point_display):
                box_center = np.array([(box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2])
                direction = box_center - point_display
                if np.allclose(direction, 0):
                    direction = np.array([0.0, 1.0])
                direction = direction / np.linalg.norm(direction)
                displacements[idx] += direction * 7

            axes_pos = record["artist"].get_position()
            if axes_pos[1] < 0.34:
                displacements[idx] += np.array([0.0, 5.0])
            if axes_pos[0] > 0.82 and axes_pos[1] < 0.42:
                displacements[idx] += np.array([-6.0, 4.0])

        max_move = 0.0
        for record, displacement in zip(text_records, displacements):
            if np.allclose(displacement, 0):
                continue
            displacement = np.clip(displacement, -14, 14)
            axes_delta = display_delta_to_axes(ax, displacement[0], displacement[1])
            current = np.array(record["artist"].get_position(), dtype=float)
            new_position = clamp_axes_position(current + axes_delta)
            record["artist"].set_position(new_position)
            max_move = max(max_move, float(np.linalg.norm(displacement)))

        if max_move < 0.5:
            break


def spread_low_projection_labels(text_records):
    low_records = [record for record in text_records if record["point"][1] < 0.28]
    if len(low_records) <= 1:
        return
    low_records = sorted(low_records, key=lambda record: record["point"][0])
    x_slots = np.linspace(0.22, 0.54, len(low_records))
    y_slots = [0.36, 0.40, 0.34, 0.43, 0.38]
    for index, (record, x_slot) in enumerate(zip(low_records, x_slots)):
        record["artist"].set_position(
            clamp_axes_position((float(x_slot), y_slots[index % len(y_slots)]))
        )


def refine_requested_labels(text_records):
    refined_positions = {
        "Miniopterus schreibersii": (0.82, 0.335),
        "Myotis ricketti": (0.49, 0.255),
        "Hipposideros pomona": (0.39, 0.300),
    }
    for record in text_records:
        if record["host"] in refined_positions:
            record["artist"].set_position(refined_positions[record["host"]])


def make_figure(df):
    x = df["Viral Diversity_norm"].to_numpy()
    y = df["Zoonotic Species Count_norm"].to_numpy()
    z = df["ZSN_corrected_norm"].to_numpy()
    sizes = df["VS_log_norm"].to_numpy()
    colors = df["CRS"].to_numpy()
    hosts = df["Host"].to_numpy()

    size_scale = 85
    marker_sizes = sizes * size_scale

    margin = 0.02
    x_min, x_max = x.min() - margin, x.max() + margin
    y_min, y_max = y.min() - margin, y.max() + margin
    z_min, z_max = z.min() - margin, z.max() + margin

    fig = plt.figure(figsize=(10.2, 8.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_position([0.03, 0.12, 0.77, 0.82])

    scatter = ax.scatter(
        x,
        y,
        z,
        s=marker_sizes,
        c=colors,
        cmap="viridis",
        alpha=0.65,
        edgecolors="white",
        linewidths=0.3,
    )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_locator(MultipleLocator(0.1))
        axis.set_major_formatter(FormatStrFormatter("%.1f"))
        axis.set_minor_locator(MultipleLocator(0.01))
        axis.set_tick_params(which="major", length=5, pad=1)
        axis.set_tick_params(which="minor", length=2)

    ax.grid(which="major", linestyle="-", linewidth=0.4, color="0.82")

    representative_sizes = np.linspace(sizes.min(), sizes.max(), 5)
    for value in representative_sizes:
        ax.scatter(
            [],
            [],
            [],
            s=max(value * size_scale, 8),
            c="0.45",
            alpha=0.65,
            edgecolors="white",
            linewidths=0.3,
            label=f"{value:.2f}",
        )
    legend = ax.legend(
        scatterpoints=1,
        title="VS_log_norm",
        loc="upper left",
        bbox_to_anchor=(0.02, 0.98),
        fontsize=7,
        title_fontsize=7,
    )
    ax.add_artist(legend)

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02, shrink=0.78, fraction=0.032)
    cbar.set_label("CRS", fontsize=9, labelpad=8)
    cbar.ax.tick_params(labelsize=8, length=3, pad=2)

    ax.set_xlabel("Viral Diversity_norm", labelpad=9, fontsize=9)
    ax.set_ylabel("Zoonotic Species Count_norm", labelpad=9, fontsize=9)
    ax.set_zlabel("ZSN_corrected_norm", labelpad=9, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.view_init(elev=25, azim=-60)

    projected_points = []
    for idx in range(len(hosts)):
        xi, yi, zi = x[idx], y[idx], z[idx]
        x2, y2, _ = proj3d.proj_transform(xi, yi, zi, ax.get_proj())
        point_axes = tuple(
            ax.transAxes.inverted().transform(ax.transData.transform((x2, y2)))
        )
        projected_points.append((float(point_axes[0]), float(point_axes[1])))

    label_layout = layout_labels(hosts, projected_points)
    text_records = []
    for label in label_layout:
        xt, yt = label["center"]
        artist = ax.text2D(
            xt,
            yt,
            label["host"],
            transform=ax.transAxes,
            fontsize=LABEL_FONT_SIZE,
            fontstyle="italic",
            ha="center",
            va="center",
            clip_on=False,
            zorder=30,
        )
        text_records.append({"artist": artist, "point": label["point"], "host": label["host"]})

    adjust_text_artists(fig, ax, text_records)
    spread_low_projection_labels(text_records)
    refine_requested_labels(text_records)

    for record in text_records:
        center = record["artist"].get_position()
        point_axes = record["point"]
        box = label_box(center, record["host"])
        line_start = line_start_for_box(point_axes, box)
        distance = float(np.hypot(center[0] - point_axes[0], center[1] - point_axes[1]))
        if distance > LINE_THRESHOLD:
            ax.add_artist(
                Line2D(
                    [line_start[0], point_axes[0]],
                    [line_start[1], point_axes[1]],
                    transform=ax.transAxes,
                    color="0.35",
                    linewidth=0.45,
                    clip_on=False,
                    zorder=25,
                )
            )

    return fig


def main():
    configure_matplotlib()
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df = load_data(DATA_PATH)
    fig = make_figure(df)
    fig.savefig(OUTPUT_PATH, format="pdf", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
