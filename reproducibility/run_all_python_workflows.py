#!/usr/bin/env python3
"""Run all released Python workflows from the repository root."""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

PYTHON_WORKFLOWS = [
    "01_sampling_and_host_representativeness/01_pool-size-distribution/scripts/plot_pool_sample_count_distribution.py",
    "01_sampling_and_host_representativeness/02_gbif-host-diversity-comparison/scripts/compare_inhouse_host_diversity_with_gbif.py",
    "01_sampling_and_host_representativeness/03_seasonal-species-accumulation-extrapolation/scripts/03_plot_terrain_species_accumulation_resampling.py",
    "02_viral_detection_and_spectrum/01_read-based-viral-family-detection-accumulation/scripts/plot_viral_family_detection_accumulation_curve.py",
    "02_viral_detection_and_spectrum/02_province-host-viral-spectrum_in-this-study/scripts/plot_viral_spectrum.py",
    "02_viral_detection_and_spectrum/03_host-genus-viral-spectrum_in-this-study/scripts/plot_viral_spectrum.py",
    "02_viral_detection_and_spectrum/04_host-genus-viral-spectrum_previous-datasets/scripts/plot_viral_spectrum.py",
    "03_virus_sharing_ecology/01_virus-sharing-distance-decay/scripts/analyze_virus_sharing_distance_decay.py",
    "03_virus_sharing_ecology/02_terrain-virus-sharing/scripts/01_plot_terrain_shared_cluster_heatmap.py",
    "03_virus_sharing_ecology/02_terrain-virus-sharing/scripts/02_test_terrain_association_with_virus_sharing.py",
    "03_virus_sharing_ecology/03_host-identity-virus-sharing/scripts/analyze_host_identity_virus_sharing.py",
    "03_virus_sharing_ecology/04_host-taxonomy-virus-sharing/scripts/analyze_host_taxonomy_virus_sharing.py",
    "04_viral_diversity_landscape_models/02_spatial-autocorrelation/scripts/calculate_global_moran_spatial_autocorrelation.py",
    "05_environmental_association_and_gamm/01_environmental-correlation-analysis/scripts/01_calculate_environmental_spearman_correlations.py",
    "05_environmental_association_and_gamm/01_environmental-correlation-analysis/scripts/02_run_shannon_environment_mantel_tests.py",
    "05_environmental_association_and_gamm/01_environmental-correlation-analysis/scripts/03_plot_environmental_correlation_mantel.py",
    "05_environmental_association_and_gamm/02_environmental-variable-standardization/scripts/standardize_environmental_variables.py",
    "05_environmental_association_and_gamm/03_environmental-correlation-vif-screening/scripts/screen_environmental_correlation_vif.py",
    "05_environmental_association_and_gamm/05_environmental-final-variable-selection/scripts/define_parsimonious_environmental_predictors.py",
    "06_rdrp_sequence_diversity_and_host_virus_network/01_viral-order-butterfly-plot/scripts/plot_viral_order_butterfly.py",
    "06_rdrp_sequence_diversity_and_host_virus_network/02_rdrp-amino-acid-identity-barplot/scripts/plot_rdrp_amino_acid_identity_barplot.py",
    "06_rdrp_sequence_diversity_and_host_virus_network/03_rdrp-contig-host-virus-network-centrality/scripts/plot_rdrp_contig_host_virus_network_centrality.py",
    "07_coronavirus_surveillance_priority/01_bat-host-surveillance-priority-assessment/scripts/01_calculate_entropy_weighted_crs.py",
    "07_coronavirus_surveillance_priority/01_bat-host-surveillance-priority-assessment/scripts/02_select_prioritized_bat_hosts.py",
    "07_coronavirus_surveillance_priority/02_bat-crs-scatterplot/scripts/plot_bat_crs_scatterplot.py",
]


def run_script(relative_path: str, index: int, total: int) -> int:
    script_path = REPOSITORY_ROOT / relative_path
    if not script_path.exists():
        print(f"[{index}/{total}] MISSING {relative_path}", flush=True)
        return 1

    start = time.time()
    print(f"[{index}/{total}] RUN {relative_path}", flush=True)
    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(script_path.parent),
        text=True,
    )
    elapsed = time.time() - start
    status = "OK" if result.returncode == 0 else "FAIL"
    print(f"[{index}/{total}] {status} {elapsed:.1f}s {relative_path}", flush=True)
    return result.returncode


def main() -> int:
    failures = []
    for index, relative_path in enumerate(PYTHON_WORKFLOWS, start=1):
        return_code = run_script(relative_path, index, len(PYTHON_WORKFLOWS))
        if return_code != 0:
            failures.append((relative_path, return_code))

    if failures:
        print("\nPython workflow failures:")
        for relative_path, return_code in failures:
            print(f"- {relative_path}: exit {return_code}")
        return 1

    print("\nAll released Python workflows completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
