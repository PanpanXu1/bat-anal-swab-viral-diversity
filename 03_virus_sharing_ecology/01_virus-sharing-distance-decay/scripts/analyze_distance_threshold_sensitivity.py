"""Repeat positive-sharing correlations after geographic distance exclusions."""

from analyze_virus_sharing_distance_decay import (
    TABLES_DIR,
    format_p_value,
    format_summary,
    load_analysis,
    plot_distance_decay,
    site_label_permutation,
)


def main() -> None:
    analysis, sharing_matrix, distance_matrix = load_analysis()
    thresholds = (None, 0, 25, 50, 100)
    summary = site_label_permutation(sharing_matrix, distance_matrix, thresholds)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    format_summary(summary).to_csv(TABLES_DIR / "distance_threshold_sensitivity.csv", index=False)
    for threshold, row in zip(thresholds, summary.to_dict("records")):
        subset = analysis if threshold is None else analysis[analysis["Distance_km"] > threshold]
        plot_distance_decay(subset, row["Spearman_rho"], row["P_value"], threshold)
        print(
            f"{row['Distance_filter']}: n={row['N_site_pairs']}, "
            f"rho={row['Spearman_rho']:.2f}, P {format_p_value(row['P_value'])}"
        )


if __name__ == "__main__":
    main()
