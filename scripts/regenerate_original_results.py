"""
Regenerate ALL original results (figs 1-24 with 10 pg default SPION loading).
Fills in missing figures from the original results/ folder.
"""

import sys
import time
from pathlib import Path

proj_root = Path(__file__).parent.parent
sys.path.insert(0, str(proj_root))

# Import all figure modules
from stent_capture.figures import (
    fig01_single_strut, fig02_ring_heatmaps, fig03_gradient_vs_distance,
    fig04_magnetisation_sweep, fig05_strut_dimensions, fig06_n_struts,
    fig07_gradient_contours, fig08_force_parameter, fig09_axial_profile,
    fig10_convergence, fig11_rz_heatmap,
    fig12_external_field_comparison,
    fig13_force_parameter, fig14_force_vs_distance, fig15_drag_vs_velocity,
    fig16_capture_map,
    fig17_spion_loading_sweep, fig18_single_trajectory, fig19_trajectory_bundle,
    fig20_capture_efficiency, fig21_static_vs_trajectory,
    fig22_concentration_field, fig23_concentration_vs_distance, fig24_time_to_threshold,
)

# List of (module, description, stage)
FIGURES = [
    (fig01_single_strut, "Single strut B-field", "Stage 1"),
    (fig02_ring_heatmaps, "Ring heatmaps", "Stage 1"),
    (fig03_gradient_vs_distance, "Gradient vs distance", "Stage 1"),
    (fig04_magnetisation_sweep, "Magnetisation sweep", "Stage 1"),
    (fig05_strut_dimensions, "Strut dimensions", "Stage 1"),
    (fig06_n_struts, "Number of struts", "Stage 1"),
    (fig07_gradient_contours, "Gradient contours", "Stage 1"),
    (fig08_force_parameter, "Force parameter", "Stage 1"),
    (fig09_axial_profile, "Axial profile", "Stage 1"),
    (fig10_convergence, "Convergence", "Stage 1"),
    (fig11_rz_heatmap, "R-Z heatmap", "Stage 1"),
    (fig12_external_field_comparison, "External field comparison", "Stage 1"),
    (fig13_force_parameter, "Force parameter (expanded)", "Stage 2"),
    (fig14_force_vs_distance, "Force vs distance", "Stage 2"),
    (fig15_drag_vs_velocity, "Drag vs velocity", "Stage 2"),
    (fig16_capture_map, "Capture map", "Stage 2"),
    (fig17_spion_loading_sweep, "SPION loading sweep", "Stage 3"),
    (fig18_single_trajectory, "Single trajectory", "Stage 3"),
    (fig19_trajectory_bundle, "Trajectory bundle", "Stage 3"),
    (fig20_capture_efficiency, "Capture efficiency", "Stage 3"),
    (fig21_static_vs_trajectory, "Static vs trajectory", "Stage 3"),
    (fig22_concentration_field, "Concentration field (VEGF)", "Stage 3c"),
    (fig23_concentration_vs_distance, "Concentration vs distance", "Stage 3c"),
    (fig24_time_to_threshold, "Time to threshold", "Stage 3c"),
]

def main():
    print(f"\n{'='*80}")
    print("REGENERATING ALL ORIGINAL RESULTS (DEFAULT 10 PG SPION LOADING)")
    print(f"Output directory: {proj_root / 'results'}")
    print(f"{'='*80}\n")

    t_start = time.time()
    success_count = 0
    fail_count = 0

    for i, (module, description, stage) in enumerate(FIGURES, 1):
        try:
            fig_name = module.__name__.split('.')[-1]
            print(f"[{i:2d}/24] {fig_name:35s} | {description:35s} | {stage:10s}", end=" ... ", flush=True)

            t0 = time.time()
            module.main()
            elapsed = time.time() - t0

            print(f"OK ({elapsed:.1f}s)")
            success_count += 1

        except Exception as e:
            print(f"FAILED: {e}")
            fail_count += 1

    total_time = time.time() - t_start

    print(f"\n{'='*80}")
    print(f"RESULTS: {success_count} OK  |  {fail_count} FAILED  |  Total time: {total_time/60:.1f} min")
    print(f"Output saved to: {proj_root / 'results'}")
    print(f"{'='*80}\n")

if __name__ == "__main__":
    main()
