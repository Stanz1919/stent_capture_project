"""
Parametric sweep: cell count vs therapeutic effectiveness, constrained by Ozawa secretion rates.

Ozawa et al. (2004) defines cell population secretion rates:
- Normal angiogenesis:  5–70 ng/10^6 cells/day
- Threshold region:     70–100 ng/10^6 cells/day
- Aberrant angiogenesis: >100 ng/10^6 cells/day

This script computes required per-cell q_cell for each category, then sweeps
across cell counts to show tissue-level VEGF concentrations.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import make_ring, DEFAULTS
from stent_capture.paracrine.transport import ParacrineField, L_DIFFUSION
from stent_capture.paracrine.secretion import VEGFSource
from stent_capture.paracrine.therapeutic import therapeutic_zone_radius


MW_VEGF = 45000  # g/mol
NA = 6.022e23    # molecules/mol


def q_cell_from_ozawa_rate(population_rate_ng_per_million_per_day: float) -> float:
    """
    Convert Ozawa population secretion rate (ng/10^6 cells/day) to per-cell rate (g/s).

    Parameters
    ----------
    population_rate_ng_per_million_per_day : float
        Rate in ng per 10^6 cells per day (Ozawa units)

    Returns
    -------
    q_cell : float
        Per-cell secretion rate in g/s
    """
    # Convert: ng/10^6 cells/day -> g/cell/s
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
    """
    Compute peak concentration and therapeutic zone radius for given cell count and Ozawa rate.
    """
    ring = make_ring()
    R = DEFAULTS["R"]
    Lx = 2.0 * np.pi * R
    Lz = 4.0 * DEFAULTS["L"]
    z_centre = Lz / 2.0

    cell_pos = _dense_cell_positions(ring, z_centre, Lx, n_per_strut)
    n_cells = len(cell_pos)

    # Compute per-cell rate from Ozawa category
    q_cell = q_cell_from_ozawa_rate(ozawa_rate)

    # Solve PDE
    field = ParacrineField(Lx=Lx, Lz=Lz, Nx=300, Nz=300)
    src = VEGFSource(q_cell=q_cell)
    field.source = src.source_field(field.X, field.Z, cell_pos)
    C = field.solve_steady_state()

    # Metrics
    C_max = np.max(C)
    r_tz = therapeutic_zone_radius(C, field.X, field.Z, cell_pos, threshold=5.0)

    return {
        "n_cells": n_cells,
        "n_per_strut": n_per_strut,
        "ozawa_rate": ozawa_rate,
        "q_cell": q_cell,
        "C_max": C_max,
        "r_tz": r_tz,
    }


# Define Ozawa categories
ozawa_categories = {
    "Normal (37.5 ng/M)": 37.5,           # Midpoint of 5–70
    "Threshold (85 ng/M)": 85.0,          # Midpoint of 70–100
    "Aberrant (150 ng/M)": 150.0,         # Typical aberrant level
}

# Sweep cell counts
n_per_strut_values = [10, 20, 40, 80, 160]

results = {cat: [] for cat in ozawa_categories}

print("\n" + "="*80)
print("OZAWA-CONSTRAINED PARAMETRIC SWEEP")
print("="*80)

for cat_name, ozawa_rate in ozawa_categories.items():
    print(f"\n{cat_name}:")
    print(f"  {'Cells':>5} | {'C_max (ng/mL)':>15} | {'Therapeutic radius (um)':>24}")
    print(f"  {'-'*5}-+-{'-'*15}-+-{'-'*24}")

    for n_per_strut in n_per_strut_values:
        metrics = compute_metrics(n_per_strut, ozawa_rate)
        results[cat_name].append(metrics)
        print(f"  {metrics['n_cells']:>5} | {metrics['C_max']:>15.3f} | {metrics['r_tz']*1e6:>24.1f}")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

colors = {"Normal (37.5 ng/M)": "#3498db", "Threshold (85 ng/M)": "#f39c12", "Aberrant (150 ng/M)": "#e74c3c"}
markers = {"Normal (37.5 ng/M)": "o", "Threshold (85 ng/M)": "s", "Aberrant (150 ng/M)": "^"}

for cat_name in ozawa_categories:
    data = results[cat_name]
    n_cells = [m["n_cells"] for m in data]
    C_max = [m["C_max"] for m in data]
    r_tz = [m["r_tz"] * 1e6 for m in data]

    ax1.plot(n_cells, C_max, color=colors[cat_name], marker=markers[cat_name],
             linestyle="-", lw=2.0, markersize=8, label=cat_name)
    ax2.plot(n_cells, r_tz, color=colors[cat_name], marker=markers[cat_name],
             linestyle="-", lw=2.0, markersize=8, label=cat_name)

ax1.axhline(5.0, color="gray", linestyle="--", lw=1.0, alpha=0.5, label="Therapeutic threshold (5 ng/mL)")
ax1.set_xlabel("Number of captured cells")
ax1.set_ylabel("Peak VEGF concentration (ng/mL)")
ax1.set_title("Peak Concentration vs Cell Count")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=9)
ax1.set_xlim(0, max([m["n_cells"] for m in sum(results.values(), [])]) * 1.1)

ax2.axhline(L_DIFFUSION * 1e6, color="gray", linestyle="-.", lw=1.0, alpha=0.5, label=f"Diffusion length ({L_DIFFUSION*1e6:.0f} µm)")
ax2.set_xlabel("Number of captured cells")
ax2.set_ylabel("Therapeutic zone radius (µm)")
ax2.set_title("Therapeutic Zone vs Cell Count")
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=9)
ax2.set_xlim(0, max([m["n_cells"] for m in sum(results.values(), [])]) * 1.1)

fig.suptitle("Figure 28: Ozawa-Constrained Parametric Sweep", fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])

import os
os.makedirs('results', exist_ok=True)
fig.savefig('results/fig28_ozawa_parametric_sweep.png', dpi=300, bbox_inches='tight')
print("\n" + "="*80)
print("[OK] fig28_ozawa_parametric_sweep.png saved")
print("="*80)
plt.close(fig)
