from pathlib import Path

import pandas as pd


WORKFLOW_ROOT = Path(__file__).resolve().parents[1]

ENV_VARS = [
    "Annual Mean Temperature (Bio1)",
    "Mean Diurnal Range (Mean of monthly max temp - min temp) (Bio2)",
    "Isothermality (BIO2/BIO7) (x100) (Bio3)",
    "Temperature Seasonality (standard deviation x100) (Bio4)",
    "Max Temperature of Warmest Month (Bio5)",
    "Min Temperature of Coldest Month (Bio6)",
    "Temperature Annual Range (BIO5-BIO6) (Bio7)",
    "Mean Temperature of Wettest Quarter (Bio8)",
    "Mean Temperature of Driest Quarter (Bio9)",
    "Mean Temperature of Warmest Quarter (Bio10)",
    "Mean Temperature of Coldest Quarter (Bio11)",
    "Annual Precipitation (Bio12)",
    "Precipitation of Wettest Month (Bio13)",
    "Precipitation of Driest Month (Bio14)",
    "Precipitation Seasonality (Coefficient of Variation) (Bio15)",
    "Precipitation of Wettest Quarter (Bio16)",
    "Precipitation of Driest Quarter (Bio17)",
    "Precipitation of Warmest Quarter (Bio18)",
    "Precipitation of Coldest Quarter (Bio19)",
    "Global Mammal Richness (GMR)",
    "Global Railway (GR)",
    "Global Linear Hydrography (GLH)",
    "Human Footprint (HFT)",
    "Fractional Vegetation Cover (FVC)",
    "Normalized Difference Vegetation Index (NDVI)",
    "China GDP Spatial Distribution (GDP)",
    "China High-Resolution Ecological Environment Quality (CHEQ)",
    "China Population Spatial Distribution (PSD)",
    "Digital Elevation Model (DEM)",
]


def main():
    input_path = WORKFLOW_ROOT / "input" / "environmental_preprocessing_input.csv"
    output_dir = WORKFLOW_ROOT / "output" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(input_path)
    missing = [var for var in ENV_VARS if var not in data.columns]
    if missing:
        raise ValueError(f"Missing environmental variables: {missing}")

    means = data[ENV_VARS].mean()
    stds = data[ENV_VARS].std(ddof=0)
    if (stds == 0).any():
        zero_sd = stds[stds == 0].index.tolist()
        raise ValueError(f"Cannot standardize variables with zero standard deviation: {zero_sd}")

    z = (data[ENV_VARS] - means) / stds
    metadata_cols = [
        "Sample group", "Shannon index", "Season", "Year", "SampleSize (log10)",
        "latitude", "longitude", "Terrain", "Host genus",
    ]
    landuse_cols = [col for col in data.columns if col.startswith("landuseName ")]
    missing_metadata = [col for col in metadata_cols if col not in data.columns]
    if missing_metadata:
        raise ValueError(f"Missing metadata columns: {missing_metadata}")
    standardized = pd.concat(
        [data[metadata_cols].reset_index(drop=True), z.reset_index(drop=True), data[landuse_cols].reset_index(drop=True)],
        axis=1,
    )
    standardized.to_csv(output_dir / "standardized_environmental_variables_by_pool.csv", index=False)


if __name__ == "__main__":
    main()
