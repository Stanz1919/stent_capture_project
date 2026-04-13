"""
Fig 13 — Force parameter: |B_total| * |∇|B_total|| vs distance.

The cell-capture force on a superparamagnetic particle is:

    F = (V_p * chi_eff / mu_0) * |B_total| * |nabla||B_total||
      = (V_p * chi_eff / (2*mu_0)) * nabla(|B_total|^2)

The force parameter FP = |B_total| * |∇|B_total|| (units T²/m) is the
field-geometry contribution to the force, independent of particle properties.
A larger FP means stronger capture force for a given SPION-labelled cell.

This figure shows why applying B0 enables cell capture even though it reduces
|∇|B_total|| (see fig12):

- |∇|B_total|| decreases when B0 is axial (B_total rotates towards z,
  suppressing the projection of dB_stent/dx onto B_total).
- But |B_total| increases from ~30 mT to ~500 mT (×17 at 200 µm).
- The product FP = |B_total| * |∇|B_total|| therefore increases substantially,
  explaining the experimental observation that external fields enhance capture.

4-panel layout:

(a) FP vs distance, log y-axis — axial B0 (0, 0.1, 0.5, 1.0 T)
(b) FP vs distance, log y-axis — transverse B0 parallel to M (+x, 0, 0.1, 0.5, 1.0 T)
(c) FP at fixed distance (200 µm) vs B0 magnitude, axial and transverse
(d) Enhancement ratio FP(B0) / FP(B0=0) vs distance for both directions

Default stent: R=1.5 mm, w=100 µm, t=80 µm, L=500 µm, M=1.0 MA/m, 8 struts.
assume_saturation=True throughout.

Run standalone::

    python -m stent_capture.figures.fig13_force_parameter
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_B0_CASES = [0.0, 0.1, 0.5, 1.0, 1.5]   # T — 1.5 T = MRI strength (COMSOL reference)


def _force_parameter_profile(
    B0_vec: np.ndarray,
    d: np.ndarray,
) -> np.ndarray:
    """
    Compute FP = |B_total| * |∇|B_total|| along the through-strut radial line.

    Parameters
    ----------
    B0_vec : (3,) array — external field vector in Tesla (zero-vector for B0=0)
    d      : 1-D array — distances from stent outer surface (m)

    Returns
    -------
    FP : 1-D array — force parameter in T²/m
    """
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    ring = make_ring()
    ring.assume_saturation = True

    ext = UniformExternalField(B0_vec) if np.any(B0_vec != 0) else None
    tf  = TotalField(ring, ext)

    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    obs_z = np.zeros_like(d)
    pts   = np.column_stack([obs_x, obs_y, obs_z])

    B_mag = np.linalg.norm(tf.field_at(pts), axis=1)
    G     = tf.grad_B(obs_x, obs_y, obs_z)
    return B_mag * G


# ---------------------------------------------------------------------------
# Colour palette: one colour per B0 magnitude, shared across panels
# ---------------------------------------------------------------------------

_COLORS = ["#333333", "#2980b9", "#e67e22", "#c0392b", "#8e44ad"]  # 0, 0.1, 0.5, 1.0, 1.5 T
_LINESTYLES_AX    = ["-",  "-",  "-",  "-",  "-"]
_LINESTYLES_TRANS = ["--", "--", "--", "--", "--"]


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    d_profile = np.linspace(5e-6, 1.5e-3, 200)
    d_um = d_profile * 1e6

    # Compute force parameter profiles for all B0 values and both orientations
    fp_axial  = {}   # keyed by B0 value (float)
    fp_trans  = {}

    for B0 in _B0_CASES:
        fp_axial[B0] = _force_parameter_profile(
            np.array([0.0, 0.0, B0]), d_profile
        )
        fp_trans[B0] = _force_parameter_profile(
            np.array([B0, 0.0, 0.0]), d_profile
        )

    # -----------------------------------------------------------------------
    # Panel (c): FP at 200 µm vs B0 magnitude (finer sweep)
    # -----------------------------------------------------------------------
    B0_range = np.linspace(0.0, 1.5, 46)   # extended to MRI strength
    idx_200 = np.argmin(np.abs(d_profile - 200e-6))

    fp_vs_B0_axial  = np.array([
        _force_parameter_profile(np.array([0.0, 0.0, b]), d_profile)[idx_200]
        for b in B0_range
    ])
    fp_vs_B0_trans  = np.array([
        _force_parameter_profile(np.array([b, 0.0, 0.0]), d_profile)[idx_200]
        for b in B0_range
    ])

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # --- Panel (a): axial B0 ---
    for B0, col in zip(_B0_CASES, _COLORS):
        label = f"B0 = {B0:.1f} T (axial)"
        ax_a.semilogy(d_um, fp_axial[B0], color=col, lw=2.0, label=label)
    ax_a.set_xlabel("Distance from stent surface (µm)")
    ax_a.set_ylabel("|B_total| · |∇|B_total|| (T²/m)")
    ax_a.set_title("(a) Force parameter — axial B0 (+z)")
    ax_a.set_xlim(0, 1500)
    ax_a.set_ylim(1e-2, 1e4)
    ax_a.legend(fontsize=8)
    ax_a.text(
        0.97, 0.05,
        "Higher = stronger capture force\nfor same SPION parameters",
        transform=ax_a.transAxes,
        ha="right", va="bottom", fontsize=7, color="#555555",
        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8),
    )

    # --- Panel (b): transverse B0 ---
    for B0, col in zip(_B0_CASES, _COLORS):
        label = f"B0 = {B0:.1f} T (transverse, +x)"
        ax_b.semilogy(d_um, fp_trans[B0], color=col, lw=2.0, ls="--", label=label)
    ax_b.set_xlabel("Distance from stent surface (µm)")
    ax_b.set_ylabel("|B_total| · |∇|B_total|| (T²/m)")
    ax_b.set_title("(b) Force parameter — transverse B0 (+x, parallel to M)")
    ax_b.set_xlim(0, 1500)
    ax_b.set_ylim(1e-2, 1e4)
    ax_b.legend(fontsize=8)

    # --- Panel (c): FP at 200 µm vs B0 magnitude ---
    ax_c.plot(B0_range, fp_vs_B0_axial,  "b-",  lw=2, label="Axial (+z)")
    ax_c.plot(B0_range, fp_vs_B0_trans,  "r--", lw=2, label="Transverse (+x)")
    ax_c.axhline(fp_axial[0.0][idx_200], color="k", ls=":", lw=1.2,
                 alpha=0.6, label="B0 = 0 baseline")
    ax_c.axvline(1.5, color="#8e44ad", ls="--", lw=1.2, alpha=0.7, label="1.5 T (MRI / COMSOL)")
    ax_c.set_xlabel("B0 magnitude (T)")
    ax_c.set_ylabel("|B_total| · |∇|B_total|| (T²/m)")
    ax_c.set_title("(c) Force parameter at 200 µm vs B0 magnitude")
    ax_c.set_xlim(0, 1.5)
    ax_c.legend(fontsize=9)
    ax_c.grid(True, alpha=0.3)

    # --- Panel (d): enhancement ratio FP(B0) / FP(B0=0) vs distance ---
    fp0 = fp_axial[0.0]
    fp0_safe = np.where(fp0 > 0, fp0, np.nan)   # avoid divide-by-zero at origin
    for B0, col in zip(_B0_CASES[1:], _COLORS[1:]):   # skip B0=0 (ratio=1)
        ratio_ax   = fp_axial[B0] / fp0_safe
        ratio_tr   = fp_trans[B0] / fp0_safe
        ax_d.semilogy(d_um, ratio_ax, color=col, lw=2.0, ls="-",
                      label=f"B0={B0:.1f}T axial")
        ax_d.semilogy(d_um, ratio_tr, color=col, lw=2.0, ls="--",
                      label=f"B0={B0:.1f}T transverse")
    ax_d.axhline(1.0, color="k", ls=":", lw=1.2, alpha=0.5, label="No enhancement")
    ax_d.set_xlabel("Distance from stent surface (µm)")
    ax_d.set_ylabel("FP(B0) / FP(B0 = 0)")
    ax_d.set_title("(d) Force parameter enhancement ratio")
    ax_d.set_xlim(0, 1500)
    ax_d.legend(fontsize=7, ncol=2)
    ax_d.text(
        0.97, 0.97,
        "Ratio > 1: B0 helps capture\n"
        "Ratio < 1: B0 hinders capture",
        transform=ax_d.transAxes,
        ha="right", va="top", fontsize=7, color="#555555",
        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8),
    )

    fig.suptitle(
        "Force Parameter FP = |B_total| · |∇|B_total||  (proportional to SPION capture force)\n"
        "(12 struts / V2-2C, R = 1.5 mm, M = 1.0 MA/m, assume_saturation = True,\n"
        "through-strut radial profile at z = 0)\n"
        "Stage 2 will convert FP to physical force F = (V_p·χ_eff/μ₀) · FP",
        fontsize=11, y=1.02,
    )
    plt.tight_layout()
    return fig


def main():
    print("  Fig 13: Force parameter...")
    fig = make_figure()
    fig.savefig(OUT / "fig13_force_parameter.png")
    fig.savefig(OUT / "fig13_force_parameter.pdf")
    plt.close(fig)
    print("  [OK] fig13_force_parameter saved")


if __name__ == "__main__":
    main()
