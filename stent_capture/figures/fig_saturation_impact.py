"""
SPION Saturation Impact Comparison.

Three-panel figure comparing constant-susceptibility (linear) vs Langevin saturation
model for SPIONs. Shows the dramatic reduction in effective chi at high fields typical
of MRI-strength stent targeting.

Panels
------
(a) chi_eff vs |B|: Analytical Langevin curve vs constant baseline
(b) Force vs distance: Radial profile along strut-aligned axis (constant chi vs Langevin)
(c) Static capture distance vs SPION loading: Constant chi vs Langevin model

Parameters
-----------
- SPION: chi_0 = 2.0 (low-field limit), M_sat = 446 kA/m (magnetite, Furlani & Ng 2006)
- Stent: 8 struts, R = 1.5 mm, 80 µm thick, M = 2.20 MA/m (COMSOL-calibrated at B0=1.5T)
- B₀ = 1.5 T axial (MRI strength)
- Blood flow: v = 0.2 m/s (MCA mean, Aaslid 1982)

Run standalone::

    python -m stent_capture.figures.fig_saturation_impact
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import (
    SPIONLabelledCell, magnetic_force, MU_0, _chi_effective,
)
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.physics.capture_criterion import capture_distance

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

_B0_Z = 1.5  # T — MRI strength
_R_VES = 1.54e-3  # m
_V_MEAN = 0.2  # m/s
_M_SAT_MAGNETITE = 446e3  # A/m

_COLOR_LANGEVIN = "#2980b9"  # blue
_COLOR_CONSTANT = "#e67e22"  # orange

# ---------------------------------------------------------------------------
# Panel (a): Analytical chi_eff vs |B|
# ---------------------------------------------------------------------------


def _make_panel_a():
    """Analytical curve: chi_eff(B) using Langevin model."""
    fig, ax = plt.subplots(figsize=(5, 4.5))

    B_range = np.linspace(0, 3.5, 200)
    chi_0 = 2.0

    # Langevin saturation
    chi_langevin = _chi_effective(chi_0, _M_SAT_MAGNETITE, B_range)

    # Constant baseline
    chi_constant = np.full_like(B_range, chi_0)

    ax.plot(B_range, chi_langevin, "-", color=_COLOR_LANGEVIN, lw=2.5,
            label="Langevin (M_sat=446 kA/m)", zorder=5)
    ax.plot(B_range, chi_constant, "--", color=_COLOR_CONSTANT, lw=2.0,
            label="Constant (linear model)", zorder=4)

    # Mark key operating points
    ax.axvline(1.5, color="gray", ls=":", lw=1.0, alpha=0.6)
    ax.text(1.5, 1.9, "B₀=1.5T\n(external)", ha="center", fontsize=8, color="gray")

    ax.axvline(3.0, color="dimgray", ls=":", lw=1.0, alpha=0.6)
    ax.text(3.0, 1.9, "B≈3T\n(stent surface)", ha="center", fontsize=8, color="dimgray")

    # Annotations at key points
    chi_at_1p5 = float(_chi_effective(chi_0, _M_SAT_MAGNETITE, np.array([1.5]))[0])
    chi_at_3 = float(_chi_effective(chi_0, _M_SAT_MAGNETITE, np.array([3.0]))[0])

    ax.scatter([1.5], [chi_at_1p5], color=_COLOR_LANGEVIN, s=80, zorder=6, edgecolor="black", linewidth=0.5)
    ax.text(1.5, chi_at_1p5 - 0.25, f"{chi_at_1p5:.2f}", ha="center", fontsize=8, fontweight="bold")

    ax.scatter([3.0], [chi_at_3], color=_COLOR_LANGEVIN, s=80, zorder=6, edgecolor="black", linewidth=0.5)
    ax.text(3.0, chi_at_3 - 0.25, f"{chi_at_3:.2f}", ha="center", fontsize=8, fontweight="bold")

    ax.set_xlabel("Magnetic field magnitude |B| (T)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Effective susceptibility χ_eff", fontsize=10, fontweight="bold")
    ax.set_xlim(0, 3.5)
    ax.set_ylim(0, 2.3)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_title("(a) Susceptibility saturation curve", fontsize=10, fontweight="bold")

    return fig, ax


# ---------------------------------------------------------------------------
# Panel (b): Radial force profile
# ---------------------------------------------------------------------------


def _make_panel_b():
    """Force vs distance: Langevin vs constant chi."""
    fig, ax = plt.subplots(figsize=(5, 4.5))

    # Setup stent and field
    ring = make_ring(B0_magnitude=_B0_Z)
    ring.assume_saturation = True
    tf_langevin = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    # For constant chi comparison
    ring_const = make_ring(B0_magnitude=None, M=1.0e6)  # Use baseline M at B0=0
    ring_const.assume_saturation = True
    tf_constant = TotalField(ring_const, UniformExternalField([0.0, 0.0, _B0_Z]))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    # Radial profile: distance from outer stent surface
    d_range = np.linspace(5e-6, 400e-6, 80)
    x_pos = r_outer + d_range

    # Cell: 50 pg SPION loading
    cell_langevin = SPIONLabelledCell(spion_mass_per_cell=50e-15,
                                     spion_sat_magnetization=_M_SAT_MAGNETITE)
    cell_constant = SPIONLabelledCell(spion_mass_per_cell=50e-15,
                                     spion_sat_magnetization=None)

    pts = np.column_stack([x_pos, np.zeros_like(x_pos), np.zeros_like(x_pos)])

    F_langevin = magnetic_force(cell_langevin, tf_langevin, pts)
    F_constant = magnetic_force(cell_constant, tf_constant, pts)

    F_mag_langevin = np.linalg.norm(F_langevin, axis=1) * 1e12  # pN
    F_mag_constant = np.linalg.norm(F_constant, axis=1) * 1e12  # pN

    ax.semilogy(d_range * 1e6, F_mag_langevin, "-o", color=_COLOR_LANGEVIN, lw=2.0, ms=4,
               label="Langevin (sat=446 kA/m)", zorder=5)
    ax.semilogy(d_range * 1e6, F_mag_constant, "--s", color=_COLOR_CONSTANT, lw=2.0, ms=4,
               label="Constant (linear model)", zorder=4)

    # Ratio annotation at 100 µm
    idx_100 = np.argmin(np.abs(d_range - 100e-6))
    ratio = F_mag_langevin[idx_100] / F_mag_constant[idx_100]
    ax.text(100, F_mag_langevin[idx_100] * 0.6, f"{ratio:.2f}x", fontsize=8, fontweight="bold")

    ax.set_xlabel("Distance from stent outer surface (µm)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Magnetic force |F| (pN)", fontsize=10, fontweight="bold")
    ax.set_xlim(0, 400)
    ax.legend(fontsize=9, loc="lower left")
    ax.grid(True, alpha=0.3, linestyle="--", which="both")
    ax.set_title("(b) Force profile (50 pg, B₀=1.5T)", fontsize=10, fontweight="bold")

    return fig, ax


# ---------------------------------------------------------------------------
# Panel (c): Static capture distance vs SPION loading
# ---------------------------------------------------------------------------


def _make_panel_c():
    """Static capture distance: Langevin vs constant chi."""
    fig, ax = plt.subplots(figsize=(5, 4.5))

    # Setup
    ring = make_ring(B0_magnitude=_B0_Z)
    ring.assume_saturation = True
    tf_langevin = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    ring_const = make_ring(B0_magnitude=None, M=1.0e6)
    ring_const.assume_saturation = True
    tf_constant = TotalField(ring_const, UniformExternalField([0.0, 0.0, _B0_Z]))

    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MEAN)

    # SPION loadings matching fig21 panel (a)
    loadings_pg = np.array([10., 30., 50., 100., 200.])
    loadings_kg = loadings_pg * 1e-15

    d_langevin = []
    d_constant = []

    for m_kg in loadings_kg:
        cell_langevin = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=_M_SAT_MAGNETITE)
        cell_constant = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=None)

        d_l = capture_distance(cell_langevin, tf_langevin, flow, direction="inward") * 1e6
        d_c = capture_distance(cell_constant, tf_constant, flow, direction="inward") * 1e6

        d_langevin.append(d_l)
        d_constant.append(d_c)

    d_langevin = np.array(d_langevin)
    d_constant = np.array(d_constant)

    ax.semilogx(loadings_pg, d_langevin, "-o", color=_COLOR_LANGEVIN, lw=2.0, ms=7,
               markerfacecolor="white", markeredgewidth=1.5,
               label="Langevin (sat=446 kA/m)", zorder=5)
    ax.semilogx(loadings_pg, d_constant, "--s", color=_COLOR_CONSTANT, lw=2.0, ms=7,
               markerfacecolor="white", markeredgewidth=1.5,
               label="Constant (linear model)", zorder=4)

    # Annotate values
    for m, dl, dc in zip(loadings_pg, d_langevin, d_constant):
        ax.text(m, dl + 5, f"{dl:.0f}", ha="center", fontsize=7, color=_COLOR_LANGEVIN)
        ax.text(m, dc - 10, f"{dc:.0f}", ha="center", fontsize=7, color=_COLOR_CONSTANT)

    ax.set_xlabel("SPION loading per cell (pg)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Static capture distance (µm)", fontsize=10, fontweight="bold")
    ax.set_xlim(7, 250)
    ax.set_ylim(0, max(d_constant) * 1.2)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.4g"))
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_title(f"(c) Static capture distance (v̄={_V_MEAN} m/s, B₀=1.5T)", fontsize=10, fontweight="bold")

    return fig, ax


# ---------------------------------------------------------------------------
# Main figure assembly
# ---------------------------------------------------------------------------


def make_figure():
    """Create three-panel comparison figure."""
    fig_a, ax_a = _make_panel_a()
    fig_b, ax_b = _make_panel_b()
    fig_c, ax_c = _make_panel_c()

    # Create combined figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Copy data from subplots to combined figure (this is a workaround)
    # Better approach: recreate directly on combined axes
    plt.close(fig_a)
    plt.close(fig_b)
    plt.close(fig_c)

    # Now recreate on combined axes
    # Panel (a)
    B_range = np.linspace(0, 3.5, 200)
    chi_0 = 2.0
    chi_langevin = _chi_effective(chi_0, _M_SAT_MAGNETITE, B_range)
    chi_constant = np.full_like(B_range, chi_0)

    axes[0].plot(B_range, chi_langevin, "-", color=_COLOR_LANGEVIN, lw=2.5,
                label="Langevin (M_sat=446 kA/m)", zorder=5)
    axes[0].plot(B_range, chi_constant, "--", color=_COLOR_CONSTANT, lw=2.0,
                label="Constant (linear model)", zorder=4)
    axes[0].axvline(1.5, color="gray", ls=":", lw=1.0, alpha=0.6)
    axes[0].axvline(3.0, color="dimgray", ls=":", lw=1.0, alpha=0.6)
    chi_at_1p5 = float(_chi_effective(chi_0, _M_SAT_MAGNETITE, np.array([1.5]))[0])
    chi_at_3 = float(_chi_effective(chi_0, _M_SAT_MAGNETITE, np.array([3.0]))[0])
    axes[0].scatter([1.5], [chi_at_1p5], color=_COLOR_LANGEVIN, s=80, zorder=6, edgecolor="black", linewidth=0.5)
    axes[0].scatter([3.0], [chi_at_3], color=_COLOR_LANGEVIN, s=80, zorder=6, edgecolor="black", linewidth=0.5)
    axes[0].set_xlabel("Magnetic field magnitude |B| (T)", fontsize=10, fontweight="bold")
    axes[0].set_ylabel("Effective susceptibility χ_eff", fontsize=10, fontweight="bold")
    axes[0].set_xlim(0, 3.5)
    axes[0].set_ylim(0, 2.3)
    axes[0].legend(fontsize=9, loc="upper right")
    axes[0].grid(True, alpha=0.3, linestyle="--")
    axes[0].set_title("(a) Susceptibility saturation curve", fontsize=10, fontweight="bold")

    # Panel (b) — force profile
    ring = make_ring(B0_magnitude=_B0_Z)
    ring.assume_saturation = True
    tf_langevin = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    ring_const = make_ring(B0_magnitude=None, M=1.0e6)
    ring_const.assume_saturation = True
    tf_constant = TotalField(ring_const, UniformExternalField([0.0, 0.0, _B0_Z]))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    d_range = np.linspace(5e-6, 400e-6, 80)
    x_pos = r_outer + d_range

    cell_langevin = SPIONLabelledCell(spion_mass_per_cell=50e-15,
                                     spion_sat_magnetization=_M_SAT_MAGNETITE)
    cell_constant = SPIONLabelledCell(spion_mass_per_cell=50e-15,
                                     spion_sat_magnetization=None)

    pts = np.column_stack([x_pos, np.zeros_like(x_pos), np.zeros_like(x_pos)])

    F_langevin = magnetic_force(cell_langevin, tf_langevin, pts)
    F_constant = magnetic_force(cell_constant, tf_constant, pts)

    F_mag_langevin = np.linalg.norm(F_langevin, axis=1) * 1e12
    F_mag_constant = np.linalg.norm(F_constant, axis=1) * 1e12

    axes[1].semilogy(d_range * 1e6, F_mag_langevin, "-o", color=_COLOR_LANGEVIN, lw=2.0, ms=4,
                    label="Langevin (sat=446 kA/m)", zorder=5)
    axes[1].semilogy(d_range * 1e6, F_mag_constant, "--s", color=_COLOR_CONSTANT, lw=2.0, ms=4,
                    label="Constant (linear model)", zorder=4)
    axes[1].set_xlabel("Distance from stent outer surface (µm)", fontsize=10, fontweight="bold")
    axes[1].set_ylabel("Magnetic force |F| (pN)", fontsize=10, fontweight="bold")
    axes[1].set_xlim(0, 400)
    axes[1].legend(fontsize=9, loc="lower left")
    axes[1].grid(True, alpha=0.3, linestyle="--", which="both")
    axes[1].set_title("(b) Force profile (50 pg, B₀=1.5T)", fontsize=10, fontweight="bold")

    # Panel (c) — static capture distance
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MEAN)
    loadings_pg = np.array([10., 30., 50., 100., 200.])
    loadings_kg = loadings_pg * 1e-15

    d_langevin = []
    d_constant = []

    for m_kg in loadings_kg:
        cell_langevin = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=_M_SAT_MAGNETITE)
        cell_constant = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=None)

        d_l = capture_distance(cell_langevin, tf_langevin, flow, direction="inward") * 1e6
        d_c = capture_distance(cell_constant, tf_constant, flow, direction="inward") * 1e6

        d_langevin.append(d_l)
        d_constant.append(d_c)

    d_langevin = np.array(d_langevin)
    d_constant = np.array(d_constant)

    axes[2].semilogx(loadings_pg, d_langevin, "-o", color=_COLOR_LANGEVIN, lw=2.0, ms=7,
                    markerfacecolor="white", markeredgewidth=1.5,
                    label="Langevin (sat=446 kA/m)", zorder=5)
    axes[2].semilogx(loadings_pg, d_constant, "--s", color=_COLOR_CONSTANT, lw=2.0, ms=7,
                    markerfacecolor="white", markeredgewidth=1.5,
                    label="Constant (linear model)", zorder=4)
    axes[2].set_xlabel("SPION loading per cell (pg)", fontsize=10, fontweight="bold")
    axes[2].set_ylabel("Static capture distance (µm)", fontsize=10, fontweight="bold")
    axes[2].set_xlim(7, 250)
    axes[2].set_ylim(0, max(d_constant) * 1.2)
    axes[2].xaxis.set_major_formatter(mticker.FormatStrFormatter("%.4g"))
    axes[2].legend(fontsize=9, loc="upper left")
    axes[2].grid(True, alpha=0.3, linestyle="--")
    axes[2].set_title(f"(c) Static capture distance (v̄={_V_MEAN} m/s, B₀=1.5T)", fontsize=10, fontweight="bold")

    fig.suptitle(
        "SPION Saturation Impact: Langevin vs Linear Susceptibility Model\n"
        "Blue: Langevin saturation (realistic for high fields) | Orange: Constant χ (linear approximation)\n"
        "At MRI strength (1.5T), Langevin model gives 5.7× lower susceptibility; force scales accordingly.",
        fontsize=11, fontweight="bold", y=1.00,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def main():
    print("  Fig: SPION saturation impact comparison...")
    fig = make_figure()
    fig.savefig(OUT / "fig_saturation_impact.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "fig_saturation_impact.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig_saturation_impact saved")


if __name__ == "__main__":
    main()
