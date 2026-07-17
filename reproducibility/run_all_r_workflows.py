#!/usr/bin/env python3
"""Run all released R workflows from the repository root."""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

R_WORKFLOWS = [
    "01_sampling_and_host_representativeness/03_seasonal-species-accumulation-extrapolation/scripts/01_plot_seasonal_species_accumulation_extrapolation.R",
    "01_sampling_and_host_representativeness/03_seasonal-species-accumulation-extrapolation/scripts/02_calculate_terrain_inext_extrapolation_metrics.R",
    "04_viral_diversity_landscape_models/01_adjusted-shannon-model/scripts/fit_adjusted_shannon_terrain_model.R",
    "05_environmental_association_and_gamm/04_environmental-linear-nonlinear-form-assessment/scripts/assess_environmental_model_forms.R",
    "05_environmental_association_and_gamm/06_environmental-driver-gamm-final-model/scripts/fit_parsimonious_gamm.R",
]


def find_rscript() -> str | None:
    return shutil.which("Rscript")


def run_script(rscript: str, relative_path: str, index: int, total: int) -> int:
    script_path = REPOSITORY_ROOT / relative_path
    if not script_path.exists():
        print(f"[{index}/{total}] MISSING {relative_path}", flush=True)
        return 1

    start = time.time()
    print(f"[{index}/{total}] RUN {relative_path}", flush=True)
    result = subprocess.run(
        [rscript, str(script_path)],
        cwd=str(script_path.parent),
        text=True,
    )
    elapsed = time.time() - start
    status = "OK" if result.returncode == 0 else "FAIL"
    print(f"[{index}/{total}] {status} {elapsed:.1f}s {relative_path}", flush=True)
    return result.returncode


def main() -> int:
    rscript = find_rscript()
    if rscript is None:
        print("Rscript was not found on PATH. Install R 4.3 or later and retry.")
        return 1

    failures = []
    for index, relative_path in enumerate(R_WORKFLOWS, start=1):
        return_code = run_script(rscript, relative_path, index, len(R_WORKFLOWS))
        if return_code != 0:
            failures.append((relative_path, return_code))

    if failures:
        print("\nR workflow failures:")
        for relative_path, return_code in failures:
            print(f"- {relative_path}: exit {return_code}")
        print("\nCheck environment.yml and workflow-level r_requirements.txt files for required R packages.")
        return 1

    print("\nAll released R workflows completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
