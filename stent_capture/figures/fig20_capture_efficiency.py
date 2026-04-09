"""
Fig 20 — Quantitative capture efficiency curves.

Two-panel figure sweeping the two key design parameters while holding all
others fixed.  Injection line spans the magnetically active near-wall region:
r = 1.20–1.45 mm, y = 0, z = −2 mm  (20 evenly-spaced injection points).

(a) Efficiency vs mean blood velocity  (v̄ = 0.01–0.50 m/s, log x-axis)
    Fixed SPION loading: 200 pg.
    Shows the critical velocity above which capture falls sharply.

(b) Efficiency vs SPION loading per cell  (5–300 pg, log x-axis)
    Fixed flow: v̄ = 0.05 m/s (MCA-representative, slow end).
    Shows the loading threshold below which capture efficiency collapses.

Parameters
----------
- Cell radius: 10 µm, χ = 2.0
- Stent: 8 struts, R = 1.5 mm, t = 80 µm, w = 100 µm, L = 500 µm, M = 1 MA/m
- B0: 0.5 T axial (+z)
- Vessel radius: 1.54 mm, η = 4 mPa·s
- Injection: 20 cells, r ∈ [1.20, 1.45] mm, y = z_start = 0/−2 mm
- z_end = 2 mm, max_time = 1.0 s

Run standalone::

    python -m stent_capture.figures.fig20_capture_efficiency
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.simulation.capture_efficiency import (
    capture_efficiency_vs_velocity,
    capture_efficiency_vs_loading,
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

_B0_Z    = 0.5
_R_VES   = 1.54e-3
_N_CELLS = 20

# Near-wall injection line: magnetically active region
_LINE_START = np.array([1.20e-3, 0.0, -2e-3])
_LINE_END   = np.array([1.45e-3, 0.0, -2e-3])

_KW_TRAJ = dict(z_end=2e-3, max_time=1.0)

# Sweep ranges (8 points each)
_VELOCITIES = np.geomspace(0.01, 0.50, 8)        # m/s
_LOADINGS   = np.array([5, 10, 20, 40, 80, 150, 200, 300], dtype=float) * 1e-15  # kg

# Reference markers
_V_REF    = 0.05    # MCA slow-end reference (m/s)
_V_MCA    = 0.20    # MCA mean (Aaslid 1982) (m/s)
_M_REF    = 10e-15  # Polyak 2008 default (kg)
_M_HIGH   = 200e-15 # Stage 3a capture case (kg)

_COLOR_V = "#2980b9"   # blue
_COLOR_L = "#27ae60"   # green


# ---------------------------------------------------------------------------
# Compute sweeps
# ---------------------------------------------------------------------------

def _run_sweeps():
    ring = make_ring()
    ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    # (a) Velocity sweep — fixed 200 pg loading
    cell_200 = SPIONLabelledCell(spion_mass_per_cell=200e-15)
    print("  (a) Velocity sweep …")
    res_v = capture_efficiency_vs_velocity(
        cell_200, tf, ring,
        velocities=_VELOCITIES,
        line_start=_LINE_START,
        line_end=_LINE_END,
        n_cells_per_velocity=_N_CELLS,
        vessel_radius=_R_VES,
        **_KW_TRAJ,
    )

    # (b) Loading sweep — fixed v = 0.05 m/s
    flow_ref = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_REF)

    def cell_factory(m_kg: float) -> SPIONLabelledCell:
        return SPIONLabelledCell(spion_mass_per_cell=m_kg)

    print("  (b) Loading sweep …")
    res_l = capture_efficiency_vs_loading(
        cell_factory, tf, flow_ref, ring,
        loadings_kg=_LOADINGS,
        line_start=_LINE_START,
        line_end=_LINE_END,
        n_cells_per_loading=_N_CELLS,
        **_KW_TRAJ,
    )

    return res_v, res_l


# ---------------------------------------------------------------------------
# Make figure
# ---------------------------------------------------------------------------

def make_figure():
    res_v, res_l = _run_sweeps()

    fig, (ax_v, ax_l) = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(wspace=0.32, left=0.08, right=0.97, top=0.88, bottom=0.14)

    # -----------------------------------------------------------------------
    # (a) Efficiency vs velocity
    # -----------------------------------------------------------------------
    v_ms   = res_v['velocities']
    eff_v  = res_v['efficiencies']
    nc_v   = res_v['n_captured']
    n_tot  = res_v['n_total']

    ax_v.semilogx(v_ms, eff_v, '-o', color=_COLOR_V, lw=2.0, ms=7,
                  markerfacecolor='white', markeredgewidth=1.8,
                  zorder=5)

    # Annotate each point with n_captured / n_total
    for v, eff, nc in zip(v_ms, eff_v, nc_v):
        ax_v.annotate(f'{nc}/{n_tot}', xy=(v, eff),
                      xytext=(0, 6), textcoords='offset points',
                      ha='center', fontsize=7, color=_COLOR_V)

    # Reference lines
    ax_v.axvline(_V_REF, color='gray', ls='--', lw=1.0, alpha=0.8,
                 label=f'v̄ = {_V_REF*100:.0f} cm/s (Stage 3a ref.)')
    ax_v.axvline(_V_MCA, color='dimgray', ls=':', lw=1.0, alpha=0.8,
                 label=f'v̄ = {_V_MCA*100:.0f} cm/s (MCA mean, Aaslid 1982)')

    ax_v.set_xlabel('Mean blood velocity $\\bar{v}$ (m/s)', fontsize=10)
    ax_v.set_ylabel('Capture efficiency (fraction)', fontsize=10)
    ax_v.set_ylim(-0.05, 1.05)
    ax_v.set_xlim(v_ms[0] * 0.8, v_ms[-1] * 1.3)
    ax_v.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2g'))
    ax_v.legend(fontsize=7, loc='upper right')
    ax_v.grid(True, which='both', alpha=0.25)
    ax_v.set_title('(a) Efficiency vs mean blood velocity\n'
                   '(200 pg SPION loading, B0 = 0.5 T)', fontsize=10)

    # -----------------------------------------------------------------------
    # (b) Efficiency vs loading
    # -----------------------------------------------------------------------
    m_pg   = res_l['loadings_pg']
    eff_l  = res_l['efficiencies']
    nc_l   = res_l['n_captured']
    n_tot_l = res_l['n_total']

    ax_l.semilogx(m_pg, eff_l, '-s', color=_COLOR_L, lw=2.0, ms=7,
                  markerfacecolor='white', markeredgewidth=1.8,
                  zorder=5)

    # Annotate
    for m, eff, nc in zip(m_pg, eff_l, nc_l):
        ax_l.annotate(f'{nc}/{n_tot_l}', xy=(m, eff),
                      xytext=(0, 6), textcoords='offset points',
                      ha='center', fontsize=7, color=_COLOR_L)

    # Reference lines
    ax_l.axvline(_M_REF * 1e15, color='gray', ls='--', lw=1.0, alpha=0.8,
                 label=f'{_M_REF*1e15:.0f} pg (Polyak 2008 default)')
    ax_l.axvline(_M_HIGH * 1e15, color='dimgray', ls=':', lw=1.0, alpha=0.8,
                 label=f'{_M_HIGH*1e15:.0f} pg (Stage 3a capture case)')

    ax_l.set_xlabel('SPION loading per cell (pg)', fontsize=10)
    ax_l.set_ylabel('Capture efficiency (fraction)', fontsize=10)
    ax_l.set_ylim(-0.05, 1.05)
    ax_l.set_xlim(m_pg[0] * 0.7, m_pg[-1] * 1.4)
    ax_l.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.4g'))
    ax_l.legend(fontsize=7, loc='upper left')
    ax_l.grid(True, which='both', alpha=0.25)
    ax_l.set_title('(b) Efficiency vs SPION loading\n'
                   '($\\bar{v}$ = 0.05 m/s, B0 = 0.5 T)', fontsize=10)

    # -----------------------------------------------------------------------
    # Shared suptitle
    # -----------------------------------------------------------------------
    r_in = (DEFAULTS["R"] - DEFAULTS["t"] / 2) * 1e3
    fig.suptitle(
        "Fig 20 — Capture efficiency: 20 cells injected along r = 1.20–1.45 mm "
        f"(magnetically active near-wall region; lumen inner radius = {r_in:.2f} mm).  "
        "Each fraction label shows n$_{\\mathrm{captured}}$ / n$_{\\mathrm{total}}$.",
        fontsize=9,
    )
    return fig, res_v, res_l


def main():
    print("  Fig 20: Capture efficiency curves…")
    fig, res_v, res_l = make_figure()
    fig.savefig(OUT / "fig20_capture_efficiency.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig20_capture_efficiency.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig20_capture_efficiency saved")

    # Print summary table
    print("\n  Velocity sweep (200 pg):")
    print(f"  {'v_mean (m/s)':>14}  {'n_capt':>6}  {'efficiency':>10}")
    for v, nc, eff in zip(res_v['velocities'], res_v['n_captured'],
                          res_v['efficiencies']):
        print(f"  {v:>14.4f}  {nc:>6}  {eff:>10.3f}")

    print("\n  Loading sweep (v = 0.05 m/s):")
    print(f"  {'loading (pg)':>14}  {'n_capt':>6}  {'efficiency':>10}")
    for m, nc, eff in zip(res_l['loadings_pg'], res_l['n_captured'],
                          res_l['efficiencies']):
        print(f"  {m:>14.1f}  {nc:>6}  {eff:>10.3f}")


if __name__ == "__main__":
    main()
