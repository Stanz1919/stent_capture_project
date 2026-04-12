"""
fig22 — Steady-state VEGF concentration field
==============================================
Side-by-side comparison: basal secretion (left) vs VEGF-enhanced cells (right,
100× secretion — modelling transfected ECs as in Chorny 2007 / Polyak 2008).

Demonstrates that basal EC secretion (Stefanini 2008) is insufficient for
therapeutic VEGF levels with ~320 captured cells, but VEGF-enhanced cells
reach the 5 ng/mL therapeutic threshold (Ozawa 2004).
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import make_ring, save_fig, DEFAULTS
from stent_capture.paracrine.transport import ParacrineField, L_DIFFUSION
from stent_capture.paracrine.secretion import VEGFSource, Q_CELL_G_PER_S
from stent_capture.paracrine.therapeutic import C_THERAPEUTIC_LOW


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
    n_total = len(cell_pos)

    # ── Basal scenario ──────────────────────────────────────────────
    field_b = ParacrineField(Lx=Lx, Lz=Lz, Nx=250, Nz=250)
    src_b = VEGFSource(q_cell=Q_CELL_G_PER_S)
    field_b.source = src_b.source_field(field_b.X, field_b.Z, cell_pos)
    C_b = field_b.solve_steady_state()

    # ── Enhanced scenario (100×) ────────────────────────────────────
    field_e = ParacrineField(Lx=Lx, Lz=Lz, Nx=250, Nz=250)
    src_e = VEGFSource(q_cell=100.0 * Q_CELL_G_PER_S)
    field_e.source = src_e.source_field(field_e.X, field_e.Z, cell_pos)
    C_e = field_e.solve_steady_state()

    # ── Plot ────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    x_um = field_b.x * 1e6
    z_um = field_b.z * 1e6

    for ax, C, title, mult in [
        (ax1, C_b, "Basal secretion (1×)", 1),
        (ax2, C_e, "VEGF-enhanced (100×)", 100),
    ]:
        vmax = max(np.percentile(C, 99.5), C_THERAPEUTIC_LOW * 2)
        im = ax.pcolormesh(z_um, x_um, C, shading="auto", cmap="inferno",
                           vmin=0, vmax=vmax)
        fig.colorbar(im, ax=ax, pad=0.02, label="ng mL$^{-1}$")

        ax.scatter(cell_pos[:, 1] * 1e6, cell_pos[:, 0] * 1e6,
                   marker=".", s=3, color="cyan", alpha=0.4)

        if np.max(C) >= C_THERAPEUTIC_LOW:
            ax.contour(z_um, x_um, C, levels=[C_THERAPEUTIC_LOW],
                       colors=["#2ecc71"], linewidths=1.2, linestyles="--")

        ax.set_xlabel("Axial position  (µm)")
        ax.set_title(f"{title}\n$C_{{max}}$ = {np.max(C):.2f} ng/mL")
        ax.text(0.03, 0.03, f"{n_total} cells",
                transform=ax.transAxes, fontsize=7, color="white",
                bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.5))

    ax1.set_ylabel("Circumferential position  (µm)")

    fig.suptitle("Fig 22 — Steady-state VEGF field on vessel-wall surface",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    save_fig(fig, "fig22_concentration_field")
    print(f"  [OK] fig22_concentration_field saved  "
          f"({n_total} cells, basal C_max = {np.max(C_b):.3f}, "
          f"enhanced C_max = {np.max(C_e):.2f} ng/mL)")
    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    generate()
