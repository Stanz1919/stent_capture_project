"""
fig23 — VEGF concentration vs distance from captured cells
===========================================================
Radial C(r) profiles for basal and enhanced (100×) secretion,
with the Ozawa et al. (2004) therapeutic window shaded.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import make_ring, save_fig, DEFAULTS
from stent_capture.paracrine.transport import ParacrineField, L_DIFFUSION
from stent_capture.paracrine.secretion import VEGFSource, Q_CELL_G_PER_S
from stent_capture.paracrine.therapeutic import (
    concentration_vs_distance,
    therapeutic_zone_radius,
    C_THERAPEUTIC_LOW,
    C_THERAPEUTIC_HIGH,
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

    cell_pos = _dense_cell_positions(ring, z_centre, Lx, n_per_strut=40)

    fig, ax = plt.subplots(figsize=(7, 5))
    r_max_plot = 3.0 * L_DIFFUSION

    for mult, label, color, ls in [
        (1, "Basal (1×)", "#3498db", "-"),
        (100, "Enhanced (100×)", "#e74c3c", "-"),
    ]:
        field = ParacrineField(Lx=Lx, Lz=Lz, Nx=300, Nz=300)
        src = VEGFSource(q_cell=mult * Q_CELL_G_PER_S)
        field.source = src.source_field(field.X, field.Z, cell_pos)
        C = field.solve_steady_state()

        r_centres, C_mean = concentration_vs_distance(
            C, field.X, field.Z, cell_pos, n_bins=100, r_max=r_max_plot,
        )
        r_um = r_centres * 1e6
        ax.plot(r_um, C_mean, color=color, ls=ls, lw=1.8, label=label)

        r_tz = therapeutic_zone_radius(C, field.X, field.Z, cell_pos,
                                       threshold=C_THERAPEUTIC_LOW)
        if r_tz > 0:
            ax.axvline(r_tz * 1e6, color=color, ls=":", lw=1.0, alpha=0.7)
            ax.annotate(f"{r_tz*1e6:.0f} µm", xy=(r_tz*1e6, C_THERAPEUTIC_LOW),
                        xytext=(r_tz*1e6 + 20, C_THERAPEUTIC_LOW * 1.5),
                        fontsize=8, color=color,
                        arrowprops=dict(arrowstyle="->", color=color, lw=0.8))

    # Therapeutic window
    ax.axhspan(C_THERAPEUTIC_LOW, C_THERAPEUTIC_HIGH,
               color="#2ecc71", alpha=0.12,
               label=f"Therapeutic ({C_THERAPEUTIC_LOW}–{C_THERAPEUTIC_HIGH} ng/mL)")
    ax.axhline(C_THERAPEUTIC_LOW, color="#2ecc71", ls="--", lw=0.8)
    ax.axhline(C_THERAPEUTIC_HIGH, color="#2ecc71", ls="--", lw=0.8)

    ax.axvline(L_DIFFUSION * 1e6, color="gray", ls="-.", lw=0.8, alpha=0.5,
               label=f"$L_D$ = {L_DIFFUSION*1e6:.0f} µm")

    ax.set_xlabel("Distance from nearest captured cell  (µm)")
    ax.set_ylabel("VEGF concentration  (ng mL$^{-1}$)")
    ax.set_title("Fig 23 — Concentration vs distance (320 cells)")
    ax.set_xlim(0, r_max_plot * 1e6)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    save_fig(fig, "fig23_concentration_vs_distance")
    print("  [OK] fig23_concentration_vs_distance saved")
    if show:
        plt.show()
    plt.close(fig)


def main():
    """Entry point for regenerate_original_results."""
    generate()


if __name__ == "__main__":
    generate()
