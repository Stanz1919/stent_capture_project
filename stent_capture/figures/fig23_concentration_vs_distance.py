"""
fig23 — VEGF concentration vs distance from captured cells
===========================================================
Radial C(r) profiles for normal and threshold secretion rates
(Ozawa et al. 2004 categories), showing concentration decay with distance.

Normal: 37.5 ng/10^6 cells/day (midpoint of normal angiogenesis range)
Threshold: 85 ng/10^6 cells/day (midpoint of threshold angiogenesis range)
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import make_ring, save_fig, DEFAULTS
from stent_capture.paracrine.transport import ParacrineField, L_DIFFUSION
from stent_capture.paracrine.secretion import VEGFSource
from stent_capture.paracrine.therapeutic import concentration_vs_distance


def q_cell_from_ozawa_rate(population_rate_ng_per_million_per_day: float) -> float:
    """Convert Ozawa population secretion rate (ng/10^6 cells/day) to per-cell rate (g/s)."""
    rate_g_per_cell_per_s = (population_rate_ng_per_million_per_day * 1e-9) / 1e6 / 86400
    return rate_g_per_cell_per_s


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

    cell_pos = _dense_cell_positions(ring, z_centre, Lx, n_per_strut=80)

    fig, ax = plt.subplots(figsize=(9, 6))
    r_max_plot = 600e-6  # 600 um (most values within first 600 um)

    for ozawa_rate, label, color, ls in [
        (37.5, "Normal (37.5 ng/M)", "#3498db", "-"),
        (85.0, "Threshold (85 ng/M)", "#e74c3c", "-"),
    ]:
        q_cell = q_cell_from_ozawa_rate(ozawa_rate)
        field = ParacrineField(Lx=Lx, Lz=Lz, Nx=600, Nz=600)
        src = VEGFSource(q_cell=q_cell)
        field.source = src.source_field(field.X, field.Z, cell_pos)
        C = field.solve_steady_state()

        r_centres, C_mean = concentration_vs_distance(
            C, field.X, field.Z, cell_pos, n_bins=100, r_max=r_max_plot,
        )
        r_um = r_centres * 1e6
        ax.plot(r_um, C_mean, color=color, ls=ls, lw=2.0, label=label)

    ax.axvline(L_DIFFUSION * 1e6, color="gray", ls="-.", lw=1.0, alpha=0.6,
               label=f"Diffusion length L_D = {L_DIFFUSION*1e6:.0f} um")

    ax.set_xlabel("Distance from nearest captured cell (um)")
    ax.set_ylabel("VEGF concentration (ng/mL)")
    ax.set_title("Figure 23: VEGF Concentration Profile vs Distance (640 cells)")
    ax.set_xlim(0, r_max_plot * 1e6)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9, loc="upper right", framealpha=0.95)
    ax.grid(True, alpha=0.3, linestyle=':')

    fig.tight_layout(pad=1.0)
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
