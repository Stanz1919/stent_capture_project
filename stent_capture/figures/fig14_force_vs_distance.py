"""
Fig 14 — Magnetic force on SPION-labelled cell vs distance from stent.

2x2 panel figure:

(a) |F_mag| vs distance, axial B0 (0, 0.1, 0.5, 1.0 T). Log y-axis, pN.
(b) |F_mag| vs distance, transverse B0 (+x, parallel to M). Log y-axis, pN.
(c) |F_mag| vs circumferential angle at fixed r = r_outer + 200 um,
    B0 = 0.5 T axial. Polar plot showing 8-strut periodicity.
    B0 = 0 reference overlaid as dashed.
(d) |F_mag| at 200 um vs SPION load (1-100 pg), three B0 values.

Default cell: 10 um radius, 200 pg iron oxide, chi = 2.0.
Default stent: R=1.5mm, 8 struts, 100x80 um, M=1 MA/m, assume_saturation=True.

Run standalone::

    python -m stent_capture.figures.fig14_force_vs_distance
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT, COLORS_CODE_CALIBRATED, COLORS_THRESHOLD_HIGH,
    COLORS_MARKER_REFERENCE
)
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_B0_CASES = [0.0, 0.1, 0.5, 1.0, 1.5]   # T — 1.5 T = MRI strength (COMSOL reference)
# Unified color palette: gray (baseline), blue (small), orange (moderate), red (strong), gray (MRI)
_COLORS = [
    COLORS_MARKER_REFERENCE,      # 0.0 T - baseline reference (gray)
    COLORS_CODE_DEFAULT,          # 0.1 T - small field (blue)
    COLORS_CODE_CALIBRATED,       # 0.5 T - moderate field (dark orange)
    COLORS_THRESHOLD_HIGH,        # 1.0 T - strong field (red)
    COLORS_MARKER_REFERENCE,      # 1.5 T - MRI reference (gray)
]
_CELL     = SPIONLabelledCell()       # default: 10 um, 200 pg, chi=2.0


def _make_tf(B0_z: float = 0.0, B0_x: float = 0.0) -> TotalField:
    # Adaptive M: at B0 = 1.5 T use COMSOL-calibrated M = 2.20 MA/m
    B0_magnitude = B0_z if B0_z > 0 else (B0_x if B0_x > 0 else 0.0)
    ring = make_ring(B0_magnitude=B0_magnitude)
    ring.assume_saturation = True
    vec = np.array([B0_x, 0.0, B0_z])
    ext = UniformExternalField(vec) if np.any(vec != 0) else None
    return TotalField(ring, ext)


def _force_profile(tf: TotalField, d: np.ndarray) -> np.ndarray:
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2
    pts = np.column_stack([d + r_outer, np.zeros_like(d), np.zeros_like(d)])
    F = magnetic_force(_CELL, tf, pts)
    return np.linalg.norm(F, axis=1) * 1e12   # pN


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    d_profile = np.linspace(5e-6, 1.5e-3, 200)
    d_um = d_profile * 1e6

    # Panels (a) and (b): force profiles
    fp_axial  = {B0: _force_profile(_make_tf(B0_z=B0),  d_profile) for B0 in _B0_CASES}
    fp_trans  = {B0: _force_profile(_make_tf(B0_x=B0),  d_profile) for B0 in _B0_CASES}

    # Panel (c): polar plot — |F_mag| vs angle at r_outer + 200 um
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_obs = R + t / 2 + 200e-6
    n_phi = 360
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    obs_polar = np.column_stack([
        r_obs * np.cos(phi),
        r_obs * np.sin(phi),
        np.zeros(n_phi),
    ])
    tf_B15 = _make_tf(B0_z=1.5)   # 1.5 T = MRI strength (COMSOL reference)
    tf_noB = _make_tf()
    F_polar_B15 = np.linalg.norm(magnetic_force(_CELL, tf_B15, obs_polar), axis=1) * 1e12
    F_polar_no  = np.linalg.norm(magnetic_force(_CELL, tf_noB, obs_polar), axis=1) * 1e12

    # Panel (d): force at 200 um vs SPION mass
    d_200 = 200e-6
    pt_200 = np.array([[R + t / 2 + d_200, 0.0, 0.0]])
    mass_pg = np.logspace(0, 2, 40)   # 1-100 pg
    fp_mass: dict[float, np.ndarray] = {}
    for B0 in [0.0, 0.5, 1.0, 1.5]:
        tf = _make_tf(B0_z=B0)
        forces = []
        for m_pg in mass_pg:
            cell = SPIONLabelledCell(spion_mass_per_cell=m_pg * 1e-15)
            F = magnetic_force(cell, tf, pt_200)
            forces.append(float(np.linalg.norm(F)) * 1e12)
        fp_mass[B0] = np.array(forces)

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    (ax_a, ax_b), (ax_c_cart, ax_d) = axes

    def _profile_panel(ax, data_dict, title, panel_label):
        for B0, col in zip(_B0_CASES, _COLORS):
            lbl = f"B0 = {B0:.1f} T"
            ax.semilogy(d_um, data_dict[B0], color=col, lw=2, label=lbl)
        ax.set_xlabel("Distance from stent surface (um)")
        ax.set_ylabel("|F_mag| (pN)")
        ax.set_title(f"({panel_label}) {title}")
        ax.set_xlim(0, 1500)
        ax.set_ylim(1e-4, 1e5)
        ax.legend(fontsize=8)

    _profile_panel(ax_a, fp_axial, "Axial B0 (+z)", "a")
    _profile_panel(ax_b, fp_trans, "Transverse B0 (+x, parallel to M)", "b")

    # Panel (c): polar plot in a Cartesian axes (polar projection)
    ax_c_cart.remove()
    ax_c = fig.add_subplot(2, 2, 3, projection="polar")
    ax_c.plot(phi, F_polar_B15, color=_COLORS[4], lw=2, label="B0 = 1.5 T axial (MRI)")
    ax_c.plot(phi, F_polar_no,  color=_COLORS[0], lw=1.5, ls="--",
              label="B0 = 0 (reference)")
    # Mark strut positions dynamically from DEFAULTS["n_struts"]
    strut_angles = np.linspace(0, 2 * np.pi, DEFAULTS["n_struts"], endpoint=False)
    r_max_polar = max(F_polar_B15.max(), F_polar_no.max()) * 1.15
    for sa in strut_angles:
        ax_c.axvline(sa, color=COLORS_MARKER_REFERENCE, lw=0.8, alpha=0.5)
    ax_c.set_title("(c) |F_mag| vs angle at r_outer + 200 um\n(B0 = 1.5 T axial, MRI)",
                   pad=15, fontsize=10)
    ax_c.set_theta_zero_location("E")
    ax_c.set_theta_direction(1)
    ax_c.legend(fontsize=7, loc="lower right")
    # Annotate strut locations
    ax_c.set_xticks(strut_angles)
    ax_c.set_xticklabels([f"{int(np.degrees(a))}°" for a in strut_angles], fontsize=7)

    # Panel (d): force vs SPION mass
    mass_colors = [_COLORS[0], _COLORS[2], _COLORS[3], _COLORS[4]]
    for (B0, col) in zip([0.0, 0.5, 1.0, 1.5], mass_colors):
        lbl = f"B0 = {B0:.1f} T" + (" (MRI)" if B0 == 1.5 else "")
        ax_d.loglog(mass_pg, fp_mass[B0], color=col, lw=2, label=lbl)
    ax_d.axvline(10, color=COLORS_MARKER_REFERENCE, ls=":", lw=1.2, label="Default (200 pg)")
    ax_d.set_xlabel("SPION load (pg iron oxide per cell)")
    ax_d.set_ylabel("|F_mag| at 200 um (pN)")
    ax_d.set_title("(d) Force vs SPION load (at 200 um, axial B0)")
    ax_d.legend(fontsize=8)
    ax_d.grid(True, which="both", alpha=0.3)

    fig.suptitle(
        "Magnetic force on SPION-labelled cell near magnetised stent\n"
        "(cell: 10 um, 200 pg iron oxide, chi = 2.0; "
        "stent: 12 struts / V2-2C, R = 1.5 mm, M = 1.0 MA/m, assume_saturation = True)",
        fontsize=11, y=0.98,
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 14: Force vs distance...")
    fig = make_figure()
    fig.savefig(OUT / "fig14_force_vs_distance.png")
    fig.savefig(OUT / "fig14_force_vs_distance.pdf")
    plt.close(fig)
    print("  [OK] fig14_force_vs_distance saved")


if __name__ == "__main__":
    main()
