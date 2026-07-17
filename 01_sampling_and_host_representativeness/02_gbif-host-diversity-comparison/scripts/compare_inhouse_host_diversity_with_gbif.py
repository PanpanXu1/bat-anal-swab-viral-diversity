from __future__ import annotations

from pathlib import Path

import pandas as pd
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from openpyxl.styles import Font


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "input"
OUTPUTS = ROOT / "output"
TABLES = OUTPUTS / "tables"
FIGURES = OUTPUTS / "figures"

INHOUSE = "#D58A2A"
GBIF = "#4E79A7"
ONLY = "#C98B1F"
TEXT = "#222222"
GRID = "#E6E1DA"

THRESHOLDS = [1, 2, 3, 5, 10, 20, 30, 50, 75, 100]


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)


def register_times() -> tuple[str, str]:
    font_dir = ROOT / "fonts"
    regular_candidates = [
        font_dir / "Times_New_Roman.ttf",
        font_dir / "times.ttf",
        Path("C:/Windows/Fonts/times.ttf"),
    ]
    bold_candidates = [
        font_dir / "Times_New_Roman_Bold.ttf",
        font_dir / "timesbd.ttf",
        Path("C:/Windows/Fonts/timesbd.ttf"),
    ]
    italic_candidates = [
        font_dir / "Times_New_Roman_Italic.ttf",
        font_dir / "timesi.ttf",
        Path("C:/Windows/Fonts/timesi.ttf"),
    ]
    regular = next((p for p in regular_candidates if p.exists()), None)
    bold = next((p for p in bold_candidates if p.exists()), None)
    italic = next((p for p in italic_candidates if p.exists()), None)
    if regular:
        pdfmetrics.registerFont(TTFont("TimesNewRoman", str(regular)))
    if bold:
        pdfmetrics.registerFont(TTFont("TimesNewRoman-Bold", str(bold)))
    if italic:
        pdfmetrics.registerFont(TTFont("TimesNewRoman-Italic", str(italic)))
    return (
        "TimesNewRoman" if regular else "Times-Roman",
        "TimesNewRoman-Bold" if bold else "Times-Bold",
        "TimesNewRoman-Italic" if italic else "Times-Italic",
    )


FONT_REG, FONT_BOLD, FONT_ITALIC = register_times()


def support_class(records: int | float | None) -> str:
    if pd.isna(records) or int(records) == 0:
        return "Absent_in_GBIF"
    records = int(records)
    if records == 1:
        return "1"
    if records <= 5:
        return "2-5"
    if records <= 20:
        return "6-20"
    if records <= 50:
        return "21-50"
    return ">50"


def pdf_color(value: str):
    return HexColor(value)


def pdf_text(c, height, x, y, text, size=9, bold=False, italic=False, fill=TEXT, anchor="start"):
    face = FONT_ITALIC if italic else (FONT_BOLD if bold else FONT_REG)
    text = str(text)
    c.setFont(face, size)
    c.setFillColor(pdf_color(fill))
    width = c.stringWidth(text, face, size)
    if anchor == "middle":
        x -= width / 2
    elif anchor == "end":
        x -= width
    c.drawString(x, height - y, text)


def rotated_label(c, height, x, y, text, size=10):
    c.saveState()
    c.translate(x, height - y)
    c.rotate(90)
    c.setFont(FONT_REG, size)
    c.setFillColor(pdf_color(TEXT))
    c.drawCentredString(0, 0, text)
    c.restoreState()


def load_inputs():
    inhouse = pd.read_csv(DATA / "inhouse_host_species_by_province.csv")
    gbif = pd.read_csv(DATA / "gbif_curated_occurrence_records.csv")
    for df in [inhouse, gbif]:
        for col in ["Species", "Genus", "Family", "Province_standard"]:
            if col in df.columns:
                df[col] = df[col].fillna("").astype(str).str.strip()
    return inhouse, gbif


def compare_taxa(inhouse: pd.DataFrame, gbif: pd.DataFrame):
    inhouse_species = set(inhouse["Species"].dropna().astype(str)) - {""}
    gbif_species = set(gbif["Species"].dropna().astype(str)) - {""}
    tax = pd.concat([inhouse[["Species", "Genus", "Family"]], gbif[["Species", "Genus", "Family"]]], ignore_index=True)
    tax = tax.drop_duplicates("Species").set_index("Species")
    species_rows = []
    for species in sorted(inhouse_species | gbif_species):
        in_inhouse = species in inhouse_species
        in_gbif = species in gbif_species
        species_rows.append({
            "Rank": "Species",
            "Taxon": species,
            "Species": species,
            "Genus": tax.at[species, "Genus"] if species in tax.index else species.split()[0],
            "Family": tax.at[species, "Family"] if species in tax.index else "",
            "In_Inhouse_data": int(in_inhouse),
            "In_GBIF": int(in_gbif),
            "Comparison_status": "Shared" if in_inhouse and in_gbif else ("Inhouse_data_only" if in_inhouse else "GBIF_only"),
        })
    species_cmp = pd.DataFrame(species_rows)

    inhouse_genus = set(inhouse["Genus"].dropna().astype(str)) - {""}
    gbif_genus = set(gbif["Genus"].dropna().astype(str)) - {""}
    genus_cmp = pd.DataFrame([{
        "Rank": "Genus",
        "Taxon": genus,
        "Species": "",
        "Genus": genus,
        "Family": "",
        "In_Inhouse_data": int(genus in inhouse_genus),
        "In_GBIF": int(genus in gbif_genus),
        "Comparison_status": "Shared" if genus in inhouse_genus and genus in gbif_genus else ("Inhouse_data_only" if genus in inhouse_genus else "GBIF_only"),
    } for genus in sorted(inhouse_genus | gbif_genus)])

    summary = pd.DataFrame([
        {"Metric": "Species richness", "Inhouse data": int((species_cmp["In_Inhouse_data"] == 1).sum()), "GBIF": int((species_cmp["In_GBIF"] == 1).sum()), "Shared": int((species_cmp["Comparison_status"] == "Shared").sum()), "Inhouse_data_only": int((species_cmp["Comparison_status"] == "Inhouse_data_only").sum()), "GBIF_only": int((species_cmp["Comparison_status"] == "GBIF_only").sum())},
        {"Metric": "Genus richness", "Inhouse data": int((genus_cmp["In_Inhouse_data"] == 1).sum()), "GBIF": int((genus_cmp["In_GBIF"] == 1).sum()), "Shared": int((genus_cmp["Comparison_status"] == "Shared").sum()), "Inhouse_data_only": int((genus_cmp["Comparison_status"] == "Inhouse_data_only").sum()), "GBIF_only": int((genus_cmp["Comparison_status"] == "GBIF_only").sum())},
    ])
    return summary, species_cmp, genus_cmp


def build_species_support(gbif: pd.DataFrame, species_cmp: pd.DataFrame):
    support = (
        gbif.groupby(["Species", "Genus", "Family"], dropna=False)
        .agg(
            GBIF_records_total=("gbifID", "nunique"),
            GBIF_province_count=("Province_standard", "nunique"),
            GBIF_province_list=("Province_standard", lambda x: "; ".join(sorted(set(v for v in x.dropna().astype(str) if v)))),
        )
        .reset_index()
    )
    support["GBIF_record_support_class"] = support["GBIF_records_total"].map(support_class)
    support = support.merge(species_cmp[["Species", "In_Inhouse_data", "Comparison_status"]], on="Species", how="left")
    support["In_Inhouse_data"] = support["In_Inhouse_data"].fillna(0).astype(int)
    support["Comparison_status"] = support["Comparison_status"].fillna("GBIF_only")
    return support.sort_values(["GBIF_records_total", "Species"], ascending=[False, True])


def build_genus_support(gbif: pd.DataFrame, genus_cmp: pd.DataFrame):
    support = (
        gbif.groupby("Genus", dropna=False)
        .agg(
            GBIF_records_total=("gbifID", "nunique"),
            GBIF_species_count=("Species", "nunique"),
            GBIF_province_count=("Province_standard", "nunique"),
            GBIF_species_list=("Species", lambda x: "; ".join(sorted(set(v for v in x.dropna().astype(str) if v)))),
            GBIF_province_list=("Province_standard", lambda x: "; ".join(sorted(set(v for v in x.dropna().astype(str) if v)))),
        )
        .reset_index()
    )
    support["GBIF_record_support_class"] = support["GBIF_records_total"].map(support_class)
    support = support.merge(genus_cmp[["Taxon", "In_Inhouse_data", "Comparison_status"]].rename(columns={"Taxon": "Genus"}), on="Genus", how="left")
    support["In_Inhouse_data"] = support["In_Inhouse_data"].fillna(0).astype(int)
    support["Comparison_status"] = support["Comparison_status"].fillna("GBIF_only")
    return support.sort_values(["GBIF_records_total", "Genus"], ascending=[False, True])


def threshold_table(support: pd.DataFrame, count_name: str):
    rows = []
    for threshold in THRESHOLDS:
        sub = support[support["GBIF_records_total"] >= threshold]
        total = int(len(sub))
        covered = int(sub["In_Inhouse_data"].sum())
        rows.append({
            "GBIF_record_threshold": f">={threshold}",
            "Threshold_min_records": threshold,
            count_name: total,
            "Covered_by_inhouse_data": covered,
            "Coverage_percent": round(covered / total * 100, 1) if total else 0.0,
        })
    return pd.DataFrame(rows)


def write_workbook(summary, species_cmp, genus_cmp, species_support, species_thresholds, genus_support, genus_thresholds, gbif):
    species_aug = species_cmp.merge(
        species_support[["Species", "GBIF_records_total", "GBIF_province_count", "GBIF_record_support_class"]],
        on="Species",
        how="left",
    )
    species_aug["GBIF_records_total"] = species_aug["GBIF_records_total"].fillna(0).astype(int)
    species_aug["GBIF_province_count"] = species_aug["GBIF_province_count"].fillna(0).astype(int)
    species_aug["GBIF_species_count"] = species_aug["In_GBIF"].astype(int)
    species_aug["GBIF_record_support_class"] = species_aug["GBIF_record_support_class"].fillna("Absent_in_GBIF")

    genus_aug = genus_cmp.merge(
        genus_support[["Genus", "GBIF_records_total", "GBIF_species_count", "GBIF_province_count", "GBIF_record_support_class"]].rename(columns={"Genus": "Taxon"}),
        on="Taxon",
        how="left",
    )
    for col in ["GBIF_records_total", "GBIF_species_count", "GBIF_province_count"]:
        genus_aug[col] = genus_aug[col].fillna(0).astype(int)
    genus_aug["GBIF_record_support_class"] = genus_aug["GBIF_record_support_class"].fillna("Absent_in_GBIF")

    taxon_comparison = pd.concat([
        species_aug[["Rank", "Taxon", "Species", "Genus", "Family", "In_Inhouse_data", "In_GBIF", "Comparison_status", "GBIF_records_total", "GBIF_species_count", "GBIF_province_count", "GBIF_record_support_class"]],
        genus_aug.assign(Species="", Genus=genus_aug["Taxon"], Family="")[["Rank", "Taxon", "Species", "Genus", "Family", "In_Inhouse_data", "In_GBIF", "Comparison_status", "GBIF_records_total", "GBIF_species_count", "GBIF_province_count", "GBIF_record_support_class"]],
    ], ignore_index=True)
    only = species_aug[species_aug["Comparison_status"] == "Inhouse_data_only"].copy()
    only["Comparison_note"] = "Present in the inhouse data table and absent from the supplied curated GBIF table."
    method_counts = gbif["Province_assignment_method"].value_counts().rename_axis("Province_assignment_method").reset_index(name="Record_count")
    province_method = pd.crosstab(gbif["Province_standard"], gbif["Province_assignment_method"]).reset_index()

    with pd.ExcelWriter(TABLES / "Inhouse_data_vs_GBIF.xlsx", engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        taxon_comparison.to_excel(writer, sheet_name="Taxon_comparison", index=False)
        only.sort_values(["Genus", "Species"]).to_excel(writer, sheet_name="Inhouse_data_only_species", index=False)
        species_support.to_excel(writer, sheet_name="GBIF_species_frequency_support", index=False)
        species_thresholds.to_excel(writer, sheet_name="Coverage_by_GBIF_frequency", index=False)
        genus_support.to_excel(writer, sheet_name="GBIF_genus_frequency_support", index=False)
        genus_thresholds.to_excel(writer, sheet_name="Genus_coverage_by_frequency", index=False)
        method_counts.to_excel(writer, sheet_name="GBIF_assignment_QC", index=False)
        province_method.to_excel(writer, sheet_name="GBIF_province_QC", index=False)
        italic_headers = {"Species", "Genus", "Family", "Taxon", "GBIF_species_list"}
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            header_map = {cell.value: cell.column for cell in sheet[1]}
            for cell in sheet[1]:
                cell.style = "Headline 3"
            for header, col_idx in header_map.items():
                if header in italic_headers:
                    for row in range(2, sheet.max_row + 1):
                        sheet.cell(row=row, column=col_idx).font = Font(italic=True)
            for col in sheet.columns:
                width = max(len(str(cell.value)) if cell.value is not None else 0 for cell in col)
                sheet.column_dimensions[col[0].column_letter].width = min(max(width + 2, 12), 50)


def draw_bar_figure(summary):
    labels = ["Species", "Genus"]
    inhouse_values = summary["Inhouse data"].tolist()
    gbifs = summary["GBIF"].tolist()
    width, height = 470, 340
    ml, mr, mt, mb = 70, 28, 34, 55
    pw, ph = width - ml - mr, height - mt - mb
    yb = mt + ph
    max_v = max(inhouse_values + gbifs) * 1.18
    c = canvas.Canvas(str(FIGURES / "inhouse_data_vs_GBIF_species_genus_richness.pdf"), pagesize=(width, height))
    for i in range(6):
        val = max_v * i / 5
        y = yb - val / max_v * ph
        c.setStrokeColor(pdf_color(GRID))
        c.line(ml, height - y, ml + pw, height - y)
        pdf_text(c, height, ml - 8, y + 3, f"{val:.0f}", 8.5, anchor="end")
    c.setStrokeColor(pdf_color(TEXT))
    c.setLineWidth(1)
    c.line(ml, height - mt, ml, height - yb)
    c.line(ml, height - yb, ml + pw, height - yb)
    group_w = pw / len(labels)
    bar_w = group_w * 0.24
    for i, label in enumerate(labels):
        center = ml + group_w * i + group_w / 2
        for j, (val, color) in enumerate([(inhouse_values[i], INHOUSE), (gbifs[i], GBIF)]):
            x0 = center + (j - 1) * bar_w + (j - 0.5) * 4
            h = val / max_v * ph
            c.setFillColor(pdf_color(color))
            c.rect(x0, height - yb, bar_w, h, stroke=0, fill=1)
            pdf_text(c, height, x0 + bar_w / 2, yb - h - 5, int(val), 8.5, anchor="middle")
        pdf_text(c, height, center, yb + 20, label, 10, anchor="middle")
    rotated_label(c, height, 20, mt + ph / 2, "Taxon richness", 10)
    lx, ly = width - mr - 72, 16
    for k, (label, color) in enumerate([("Inhouse data", INHOUSE), ("GBIF", GBIF)]):
        y = ly + k * 15
        c.setFillColor(pdf_color(color))
        c.rect(lx, height - y - 8, 8, 8, stroke=0, fill=1)
        pdf_text(c, height, lx + 12, y, label, 8.5)
    c.save()


def draw_inhouse_data_only_figure(species_cmp):
    only = species_cmp[species_cmp["Comparison_status"] == "Inhouse_data_only"].copy()
    counts = only.groupby("Genus")["Species"].nunique().sort_values(ascending=False)
    labels = counts.index.tolist()
    values = counts.values.tolist()
    width, height = 480, 330
    ml, mr, mt, mb = 80, 24, 28, 68
    pw, ph = width - ml - mr, height - mt - mb
    yb = mt + ph
    max_v = max(values + [1]) * 1.30
    c = canvas.Canvas(str(FIGURES / "inhouse_data_only_species_absent_from_GBIF.pdf"), pagesize=(width, height))
    for i in range(int(max_v) + 1):
        y = yb - i / max_v * ph
        c.setStrokeColor(pdf_color(GRID))
        c.line(ml, height - y, ml + pw, height - y)
        pdf_text(c, height, ml - 8, y + 3, str(i), 8.5, anchor="end")
    c.setStrokeColor(pdf_color(TEXT))
    c.line(ml, height - mt, ml, height - yb)
    c.line(ml, height - yb, ml + pw, height - yb)
    group_w = pw / len(labels)
    bar_w = group_w * 0.52
    for label, val, i in zip(labels, values, range(len(labels))):
        center = ml + group_w * i + group_w / 2
        h = val / max_v * ph
        c.setFillColor(pdf_color(ONLY))
        c.rect(center - bar_w / 2, height - yb, bar_w, h, stroke=0, fill=1)
        pdf_text(c, height, center, yb - h - 5, int(val), 8.5, anchor="middle")
        pdf_text(c, height, center, yb + 18, label, 9.5, italic=True, anchor="middle")
    rotated_label(c, height, 22, mt + ph / 2, "Number of Inhouse-data-only species", 9.5)
    c.save()


def draw_threshold_figure(table, count_col, rank_label, filename):
    labels = table["GBIF_record_threshold"].tolist()
    coverage = table["Coverage_percent"].tolist()
    counts = table[count_col].tolist()
    width, height = 560, 350
    ml, mr, mt, mb = 74, 24, 30, 64
    pw, ph = width - ml - mr, height - mt - mb
    yb = mt + ph
    c = canvas.Canvas(str(FIGURES / filename), pagesize=(width, height))
    for val in range(0, 101, 20):
        y = yb - val / 100 * ph
        c.setStrokeColor(pdf_color(GRID))
        c.line(ml, height - y, ml + pw, height - y)
        pdf_text(c, height, ml - 8, y + 3, val, 8.5, anchor="end")
    c.setStrokeColor(pdf_color(TEXT))
    c.line(ml, height - mt, ml, height - yb)
    c.line(ml, height - yb, ml + pw, height - yb)
    x_start = ml + 20
    x_end = ml + pw - 4
    step = (x_end - x_start) / (len(labels) - 1)
    points = []
    for i, val in enumerate(coverage):
        x = x_start + i * step
        y = yb - val / 100 * ph
        points.append((x, y))
    c.setStrokeColor(pdf_color(INHOUSE))
    c.setLineWidth(2)
    for p1, p2 in zip(points[:-1], points[1:]):
        c.line(p1[0], height - p1[1], p2[0], height - p2[1])
    for x, y, val, n, label in zip([p[0] for p in points], [p[1] for p in points], coverage, counts, labels):
        c.setFillColor(pdf_color(INHOUSE))
        c.circle(x, height - y, 3.5, stroke=0, fill=1)
        pdf_text(c, height, x, y - 10, f"{val:.1f}", 8, anchor="middle")
        pdf_text(c, height, x, y + 10, f"n={n}", 7.5, fill="#555555", anchor="middle")
        pdf_text(c, height, x, yb + 20, label, 8.5, anchor="middle")
    pdf_text(c, height, ml + pw / 2, height - 22, "GBIF minimum occurrence records threshold", 9.5, anchor="middle")
    rotated_label(c, height, 22, mt + ph / 2, f"Inhouse data coverage of GBIF {rank_label} (%)", 9.5)
    c.save()


def main():
    ensure_dirs()
    inhouse, gbif = load_inputs()
    summary, species_cmp, genus_cmp = compare_taxa(inhouse, gbif)
    species_support = build_species_support(gbif, species_cmp)
    genus_support = build_genus_support(gbif, genus_cmp)
    species_thresholds = threshold_table(species_support, "GBIF_species_count")
    genus_thresholds = threshold_table(genus_support, "GBIF_genus_count")
    write_workbook(summary, species_cmp, genus_cmp, species_support, species_thresholds, genus_support, genus_thresholds, gbif)
    draw_bar_figure(summary)
    draw_inhouse_data_only_figure(species_cmp)
    draw_threshold_figure(species_thresholds, "GBIF_species_count", "species", "species_coverage_by_GBIF_record_threshold.pdf")
    draw_threshold_figure(genus_thresholds, "GBIF_genus_count", "genera", "genus_coverage_by_GBIF_record_threshold.pdf")


if __name__ == "__main__":
    main()
