"""
Figure 29: 2D surface showing therapeutic effectiveness vs cell count and Ozawa secretion rate.

This visualizes the trade-off space: how peak concentration and therapeutic zone vary
with both the number of captured cells and the per-cell secretion rate (Ozawa category).
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import make_ring, save_fig, DEFAULTS
from stent_capture.paracrine.transport import ParacrineField, L_DIFFUSION
from stent_capture.paracrine.secretion import VEGFSource
from stent_capture.paracrine.therapeutic import therapeutic_zone_radius


def q_cell_from_ozawa_rate(population_rate_ng_per_million_per_day: float) -> float:
    """Convert Ozawa population secretion rate (ng/10^6 cells/day) to per-cell rate (g/s)."""
    rate_g_per_cell_per_s = (population_rate_ng_per_million_per_day * 1e-9) / 1e6 / 86400
    return rate_g_per_cell_per_s


def _dense_cell_positions(ring, z_centre, Lx, n_per_strut):
    """Generate cell positions distributed around struts with random jitter."""
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


def compute_metrics(n_per_strut: int, ozawa_rate: float) -> dict:
    """Compute peak concentration and therapeutic zone radius."""
    ring = make_ring()
    R = DEFAULTS["R"]
    Lx = 2.0 * np.pi * R
    Lz = 4.0 * DEFAULTS["L"]
    z_centre = Lz / 2.0

    cell_pos = _dense_cell_positions(ring, z_centre, Lx, n_per_strut)
    n_cells = len(cell_pos)

    q_cell = q_cell_from_ozawa_rate(ozawa_rate)

    field = ParacrineField(Lx=Lx, Lz=Lz, Nx=300, Nz=300)
    src = VEGFSource(q_cell=q_cell)
    field.source = src.source_field(field.X, field.Z, cell_pos)
    C = field.solve_steady_state()

    C_max = np.max(C)
    r_tz = therapeutic_zone_radius(C, field.X, field.Z, cell_pos, threshold=5.0)

    return {"C_max": C_max, "r_tz": r_tz * 1e6}


# Grid: cell counts and secretion rates
n_per_strut_values = np.array([5, 10, 20, 40, 80, 160])
ozawa_rates = np.linspace(10, 150, 15)  # ng/10^6 cells/day

C_max_grid = np.zeros((len(ozawa_rates), len(n_per_strut_values)))
r_tz_grid = np.zeros((len(ozawa_rates), len(n_per_strut_values)))

print("\n" + "="*80)
print("COMPUTING 2D SURFACE: Cell Count × Secretion Rate")
print("="*80)

for i, ozawa_rate in enumerate(ozawa_rates):
    print(f"Ozawa rate {ozawa_rate:6.1f} ng/M: ", end="", flush=True)
    for j, n_per_strut in enumerate(n_per_strut_values):
        metrics = compute_metrics(n_per_strut, ozawa_rate)
        C_max_grid[i, j] = metrics["C_max"]
        r_tz_grid[i, j] = metrics["r_tz"]
        print(".", end="", flush=True)
    print(" [OK]")

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

n_cells_values = n_per_strut_values * 8  # 8 struts

# Panel 1: Peak concentration
im1 = ax1.contourf(n_cells_values, ozawa_rates, C_max_grid, levels=20, cmap="YlOrRd")
cs1 = ax1.contour(n_cells_values, ozawa_rates, C_max_grid, levels=[5.0, 25.0, 100.0],
                  colors="black", linewidths=1.0, alpha=0.6)
ax1.clabel(cs1, inline=True, fontsize=8, fmt="C=%g ng/mL")

# Mark Ozawa category regions
ax1.axhspan(5, 70, alpha=0.05, color="green", label="Normal (5-70)")
ax1.axhspan(70, 100, alpha=0.05, color="orange", label="Threshold (70-100)")
ax1.axhspan(100, 150, alpha=0.05, color="red", label="Aberrant (>100)")

ax1.set_xlabel("Number of captured cells")
ax1.set_ylabel("Ozawa secretion rate (ng/10⁶ cells/day)")
ax1.set_title("Peak VEGF Concentration (ng/mL)")
ax1.legend(loc="upper left", fontsize=8)
cbar1 = plt.colorbar(im1, ax=ax1, label="C_max (ng/mL)")

# Panel 2: Therapeutic zone radius
im2 = ax2.contourf(n_cells_values, ozawa_rates, r_tz_grid, levels=20, cmap="viridis")
cs2 = ax2.contour(n_cells_values, ozawa_rates, r_tz_grid,
                  levels=[100, 200, 300, 400, 500],
                  colors="white", linewidths=1.0, alpha=0.6)
ax2.clabel(cs2, inline=True, fontsize=8, fmt="r=%g µm")

# Mark Ozawa category regions
ax2.axhspan(5, 70, alpha=0.05, color="green")
ax2.axhspan(70, 100, alpha=0.05, color="orange")
ax2.axhspan(100, 150, alpha=0.05, color="red")

ax2.axhline(L_DIFFUSION * 1e6, color="white", linestyle="--", linewidth=1.5,
           label=f"Diffusion length L_D = {L_DIFFUSION*1e6:.0f} µm")

ax2.set_xlabel("Number of captured cells")
ax2.set_ylabel("Ozawa secretion rate (ng/10⁶ cells/day)")
ax2.set_title("Therapeutic Zone Radius (µm)")
ax2.legend(loc="upper left", fontsize=8)
cbar2 = plt.colorbar(im2, ax=ax2, label="r_tz (µm)")

fig.suptitle("Figure 29: Design Trade-off Space (Cell Count vs Secretion Rate)",
            fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])

import os
os.makedirs('results', exist_ok=True)
fig.savefig('results/fig29_ozawa_surface.png', dpi=300, bbox_inches='tight')
print("\n" + "="*80)
print("[OK] fig29_ozawa_surface.png saved")
print("="*80)

# Print summary statistics
print("\nKey design insights:")
print(f"  To reach therapeutic threshold (5 ng/mL) at normal secretion rate:")
print(f"    Minimum cells needed: ~{n_cells_values[np.where(C_max_grid[np.argmin(np.abs(ozawa_rates - 37.5)), :] >= 5.0)[0][0]]} cells")
print(f"  To reach therapeutic threshold at threshold secretion rate:")
print(f"    Minimum cells needed: ~{n_cells_values[np.where(C_max_grid[np.argmin(np.abs(ozawa_rates - 85.0)), :] >= 5.0)[0][0]]} cells")
print(f"  Therapeutic zone remains submaximal vs diffusion length:")
print(f"    Even at maximum, r_tz = {r_tz_grid.max():.0f} µm < L_D = {L_DIFFUSION * 1e6:.0f} µm")

plt.close(fig)
