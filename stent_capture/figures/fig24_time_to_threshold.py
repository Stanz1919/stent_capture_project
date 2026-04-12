"""
fig24 — Time to reach therapeutic VEGF threshold
=================================================
Bar chart: time for VEGF concentration to reach 5 ng/mL at 50, 100, and
200 µm from the nearest captured cell.  Uses enhanced (100×) secretion
scenario (VEGF-transfected cells), since basal secretion does not reach
the therapeutic threshold.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import make_ring, save_fig, DEFAULTS
from stent_capture.paracrine.transport import ParacrineField, L_DIFFUSION
from stent_capture.paracrine.secretion import VEGFSource, Q_CELL_G_PER_S
from stent_capture.paracrine.therapeutic import (
    time_to_threshold,
    C_THERAPEUTIC_LOW,
)


def _dense_cell_positions(ring, z_centre, Lx, n_per_strut=40):
    base = VEGFSource.cell_positions_on_ring(ring, np.full(ring.n_struts, z_centre))
    rng = np.random.default_rng(42)
    sigma_spread = 50e-6
    all_pos = []
    for bx, bz in base:
        offsets = rng.normal(scale=sigma_spread, size=(n_per_strut, 2))
        pts = np.column_stack([bx + offsets[:, 0], bz + offsets[:, 1]])
        pts[:, 0] = pts[:, 0] % Lx
        all_pos.append(pts)
    return np.vstack(all_pos)


def generate(show: bool = False) -> None:
    ring = make_ring()
    R = DEFAULTS["R"]
    Lx = 2.0 * np.pi * R
    Lz = 4.0 * DEFAULTS["L"]
    z_centre = Lz / 2.0

    field = ParacrineField(Lx=Lx, Lz=Lz, Nx=150, Nz=150)
    src = VEGFSource(q_cell=100.0 * Q_CELL_G_PER_S)
    cell_pos = _dense_cell_positions(ring, z_centre, Lx, n_per_strut=40)
    field.source = src.source_field(field.X, field.Z, cell_pos)

    t_final = 5.0 / field.k_deg
    print(f"  Transient solve: t_final = {t_final:.0f} s ({t_final/3600:.1f} h) …",
          flush=True)

    C_snaps, t_snaps = field.solve_transient(t_final=t_final, n_snapshots=200)

    ref_cell = cell_pos[0]
    probe_dists = np.array([50e-6, 100e-6, 200e-6])
    probes = np.column_stack([
        np.full(len(probe_dists), ref_cell[0]),
        ref_cell[1] + probe_dists,
    ])

    t_thresh = time_to_threshold(
        C_snaps, t_snaps, field.X, field.Z, probes,
        threshold=C_THERAPEUTIC_LOW,
    )

    # ── Plot ──────────────────────────────────────────────��─────────
    labels = [f"{d*1e6:.0f} µm" for d in probe_dists]
    t_min = t_thresh / 60.0

    fig, ax = plt.subplots(figsize=(6, 5))
    colors = ["#2ecc71", "#f39c12", "#e74c3c"]
    bars = ax.bar(labels, t_min, color=colors, edgecolor="k", linewidth=0.6)

    for bar, val in zip(bars, t_min):
        if np.isfinite(val):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    f"{val:.1f} min", ha="center", va="bottom", fontsize=9)
        else:
            y_mid = ax.get_ylim()[1] * 0.5 if ax.get_ylim()[1] > 0 else 10
            ax.text(bar.get_x() + bar.get_width() / 2, y_mid,
                    "not reached", ha="center", va="center", fontsize=9,
                    color="#e74c3c", fontweight="bold")

    ax.set_xlabel("Distance from nearest captured cell")
    ax.set_ylabel("Time to reach threshold  (min)")
    ax.set_title(
        f"Fig 24 — Time to {C_THERAPEUTIC_LOW} ng/mL threshold\n"
        f"(VEGF-enhanced 100× secretion, 320 cells)"
    )

    ax.text(0.97, 0.97,
            f"$L_D$ = {L_DIFFUSION*1e6:.0f} µm\n"
            f"$t_{{1/2}}$ = 60 min\n"
            f"q = 100× basal",
            transform=ax.transAxes, fontsize=8,
            ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    fig.tight_layout()
    save_fig(fig, "fig24_time_to_threshold")

    for d, t in zip(probe_dists, t_thresh):
        status = f"{t/60:.1f} min" if np.isfinite(t) else "NOT REACHED"
        print(f"    {d*1e6:.0f} µm  ->  {status}")
    print("  [OK] fig24_time_to_threshold saved")

    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    generate()
