"""
Fig 15 — Stokes drag on a cell in Poiseuille blood flow.

2x2 panel figure:

(a) Stokes drag vs distance from vessel wall, for mean velocities
    0.05, 0.1, 0.2, 0.5 m/s. Y-axis in piconewtons.
(b) Velocity profile across vessel cross-section (linear r-axis).
(c) Drag at the stent inner surface (r ~ R - t/2) vs mean velocity,
    with horizontal reference lines showing F_mag at the same location.
(d) Wall shear stress vs mean velocity.

Vessel: R = 1.54 mm, viscosity = 4 mPa·s.
Cell: 10 um radius.
Stent: same as Stage 1 defaults.

Run standalone::

    python -m stent_capture.figures.fig15_drag_vs_velocity
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT, COLORS_CODE_CALIBRATED, COLORS_THRESHOLD,
    COLORS_THRESHOLD_HIGH, COLORS_MARKER_REFERENCE
)
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force
from stent_capture.physics.hydrodynamics import BloodFlow, stokes_drag

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_V_CASES  = [0.05, 0.1, 0.2, 0.5]   # mean velocities (m/s)
# Unified color palette: blue (low), green (medium), orange (moderate-high), red (high)
_COLORS_V = [
    COLORS_CODE_DEFAULT,          # 0.05 m/s - low velocity (blue)
    COLORS_THRESHOLD,             # 0.1 m/s - medium velocity (green)
    COLORS_CODE_CALIBRATED,       # 0.2 m/s - moderate-high velocity (dark orange)
    COLORS_THRESHOLD_HIGH,        # 0.5 m/s - high velocity (red)
]
_R_VES    = 1.54e-3   # vessel radius
_CELL     = SPIONLabelledCell()


def _make_flow(v_mean: float) -> BloodFlow:
    return BloodFlow(vessel_radius=_R_VES, mean_velocity=v_mean)


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    R_ves = _R_VES
    R_stent = DEFAULTS["R"]
    t = DEFAULTS["t"]

    # -----------------------------------------------------------------------
    # Panel (a): drag vs distance from vessel wall
    # Points are r = R_ves - d, i.e. d = R_ves - r from wall inward
    # -----------------------------------------------------------------------
    d_from_wall = np.linspace(0.5e-6, R_ves * 0.99, 300)   # 0.5 um to ~lumen centre
    r_vals = R_ves - d_from_wall
    pts_a = np.column_stack([r_vals, np.zeros_like(r_vals), np.zeros_like(r_vals)])

    drag_profiles = {}
    for v_mean in _V_CASES:
        flow = _make_flow(v_mean)
        F = stokes_drag(_CELL, flow, pts_a)
        drag_profiles[v_mean] = np.linalg.norm(F, axis=1) * 1e12   # pN

    # -----------------------------------------------------------------------
    # Panel (b): velocity profile across vessel diameter
    # -----------------------------------------------------------------------
    r_profile = np.linspace(0, R_ves, 300)
    v_profiles = {}
    for v_mean in _V_CASES:
        flow = _make_flow(v_mean)
        pts_b = np.column_stack([r_profile, np.zeros_like(r_profile),
                                  np.zeros_like(r_profile)])
        v_profiles[v_mean] = flow.velocity_at(pts_b)[:, 2] * 100   # cm/s

    # -----------------------------------------------------------------------
    # Panel (c): drag at stent inner surface vs mean velocity
    # also: F_mag at same location for B0 = 0, 0.5 T
    # -----------------------------------------------------------------------
    r_inner = R_stent - t / 2   # stent inner surface
    v_range = np.linspace(0.01, 0.6, 60)
    pt_c = np.array([[r_inner, 0.0, 0.0]])

    drag_vs_v = []
    for v_mean in v_range:
        flow = _make_flow(v_mean)
        F = stokes_drag(_CELL, flow, pt_c)
        drag_vs_v.append(float(np.linalg.norm(F)) * 1e12)
    drag_vs_v = np.array(drag_vs_v)

    # F_mag reference lines at stent inner surface
    ring = make_ring()
    ring.assume_saturation = True
    fmag_ref = {}
    for B0 in [0.0, 0.5, 1.0]:
        vec = np.array([0.0, 0.0, B0])
        ext = UniformExternalField(vec) if B0 > 0 else None
        tf = TotalField(ring, ext)
        F = magnetic_force(_CELL, tf, pt_c)
        fmag_ref[B0] = float(np.linalg.norm(F)) * 1e12

    # -----------------------------------------------------------------------
    # Panel (d): wall shear stress vs mean velocity
    # -----------------------------------------------------------------------
    wss_vs_v = np.array([
        _make_flow(v).wall_shear_stress for v in v_range
    ])
    # Wall shear stress reference ranges for cerebral arteries
    # Atheroprotective: 1.5-7 Pa; Atherogenic (low): < 0.4 Pa

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # --- Panel (a) ---
    for v_mean, col in zip(_V_CASES, _COLORS_V):
        ax_a.semilogy(d_from_wall * 1e6, drag_profiles[v_mean],
                      color=col, lw=2, label=f"v_mean = {v_mean:.2f} m/s")
    ax_a.set_xlabel("Distance from vessel wall (um)")
    ax_a.set_ylabel("Stokes drag |F_drag| (pN)")
    ax_a.set_title("(a) Drag on 10 um cell vs distance from wall")
    ax_a.set_xlim(0, R_ves * 1e6)
    ax_a.legend(fontsize=8)
    # Mark stent inner surface distance from wall
    d_inner_from_wall = R_ves - r_inner
    ax_a.axvline(d_inner_from_wall * 1e6, color=COLORS_MARKER_REFERENCE, ls=":",
                 lw=1.2, label=f"Stent inner surface ({d_inner_from_wall*1e6:.0f} um)")
    ax_a.legend(fontsize=7)

    # --- Panel (b) ---
    for v_mean, col in zip(_V_CASES, _COLORS_V):
        ax_b.plot(r_profile * 1e3, v_profiles[v_mean],
                  color=col, lw=2, label=f"v_mean = {v_mean:.2f} m/s")
    ax_b.axvline(r_inner * 1e3, color=COLORS_MARKER_REFERENCE, ls=":", lw=1.2)
    ax_b.text(r_inner * 1e3 + 0.02, ax_b.get_ylim()[1] * 0.5 if ax_b.get_ylim()[1] > 0 else 10,
              "Stent\ninner\nsurface", fontsize=7, color=COLORS_MARKER_REFERENCE)
    ax_b.set_xlabel("Radial position r (mm)")
    ax_b.set_ylabel("v_z (cm/s)")
    ax_b.set_title("(b) Poiseuille velocity profile")
    ax_b.set_xlim(0, R_ves * 1e3)
    ax_b.legend(fontsize=8)

    # --- Panel (c) ---
    ax_c.semilogy(v_range * 100, drag_vs_v, color=COLORS_MARKER_REFERENCE, lw=2.5,
                 label="Stokes drag at stent inner surface")
    ref_colors = [COLORS_MARKER_REFERENCE, COLORS_CODE_CALIBRATED, COLORS_THRESHOLD_HIGH]
    for (B0, col) in zip([0.0, 0.5, 1.0], ref_colors):
        ax_c.axhline(fmag_ref[B0], color=col, ls="--", lw=1.5,
                     label=f"|F_mag| at stent inner surface (B0={B0:.1f}T)")
    ax_c.set_xlabel("Mean blood velocity (cm/s)")
    ax_c.set_ylabel("Force (pN)")
    ax_c.set_title("(c) Drag vs velocity — crossing F_mag lines = capture limit")
    ax_c.legend(fontsize=7)
    ax_c.set_xlim(v_range[0] * 100, v_range[-1] * 100)
    ax_c.grid(True, which="both", alpha=0.3)

    # --- Panel (d) ---
    ax_d.plot(v_range * 100, wss_vs_v, color=COLORS_CODE_DEFAULT, lw=2)
    ax_d.axhspan(1.5, 7.0, alpha=0.08, color=COLORS_THRESHOLD, label="Atheroprotective (1.5-7 Pa)")
    ax_d.axhspan(0.0, 0.4, alpha=0.1, color=COLORS_THRESHOLD_HIGH, label="Atherogenic (< 0.4 Pa)")
    ax_d.set_xlabel("Mean blood velocity (cm/s)")
    ax_d.set_ylabel("Wall shear stress (Pa)")
    ax_d.set_title("(d) Wall shear stress vs mean velocity")
    ax_d.legend(fontsize=8)
    ax_d.grid(True, alpha=0.3)
    ax_d.set_xlim(v_range[0] * 100, v_range[-1] * 100)
    ax_d.set_ylim(0)

    fig.suptitle(
        "Hydrodynamic drag on a 10 um endothelial cell in Poiseuille blood flow\n"
        "(vessel R = 1.54 mm, viscosity = 4 mPa·s, density = 1060 kg/m³ — "
        "cerebral arterial / MCA-representative; Aaslid et al. 1982)\n"
        "Stent inner surface at r = R - t/2 = 1.46 mm",
        fontsize=11, y=0.98,
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 15: Drag vs velocity...")
    fig = make_figure()
    fig.savefig(OUT / "fig15_drag_vs_velocity.png")
    fig.savefig(OUT / "fig15_drag_vs_velocity.pdf")
    plt.close(fig)
    print("  [OK] fig15_drag_vs_velocity saved")


if __name__ == "__main__":
    main()
