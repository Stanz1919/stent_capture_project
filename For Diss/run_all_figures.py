"""
run_all_figures.py
==================
Standalone driver that regenerates every dissertation appendix figure
(fig7, fig7b, fig7extra, fig8, fig9, fig10, fig11, fig12) into
``For Diss/results/``.

Usage (from inside the ``For Diss`` folder)::

    python run_all_figures.py

Long-running steps (fig11 and fig12) integrate many trajectories; allow
several minutes on a laptop-class CPU. All figures are saved as PNG+PDF.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from stent_capture.figures import (
    fig07_comsol_gradient_validation,
    fig07extra_comsol_multigeometry,
    fig08_force_parameter,
    fig09_spion_loading_sweep,
    fig10_single_trajectory,
    fig11_static_vs_trajectory,
    fig12_capture_efficiency,
)


FIGURES = [
    ("Fig 7 / 7b — COMSOL gradient validation",    fig07_comsol_gradient_validation.main),
    ("Fig 7 (extra) — COMSOL multi-geometry",      fig07extra_comsol_multigeometry.main),
    ("Fig 8 — Force parameter",                    fig08_force_parameter.main),
    ("Fig 9 — SPION loading sweep",                fig09_spion_loading_sweep.main),
    ("Fig 10 — Single-cell trajectory",            fig10_single_trajectory.main),
    ("Fig 11 — Static vs trajectory capture",      fig11_static_vs_trajectory.main),
    ("Fig 12 — Capture efficiency curves",         fig12_capture_efficiency.main),
]


def main() -> None:
    t_start = time.time()
    for name, fn in FIGURES:
        print()
        print("=" * 72)
        print(name)
        print("=" * 72)
        t0 = time.time()
        fn()
        print(f"  ({time.time() - t0:.1f} s)")
    print()
    print("=" * 72)
    print(f"All figures done in {time.time() - t_start:.0f} s")
    print("=" * 72)


if __name__ == "__main__":
    main()
