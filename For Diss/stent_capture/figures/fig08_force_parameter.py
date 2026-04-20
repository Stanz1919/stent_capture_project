"""
Fig 8 — Force parameter: |B_total| * |∇|B_total|| vs distance.

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
- But |B_total| increases from ~30 mT to ~500 mT (x17 at 200 um).
- The product FP = |B_total| * |∇|B_total|| therefore increases substantially,
  explaining the experimental observation that external fields enhance capture.

4-panel layout:

(a) FP vs distance, log y-axis — axial B0 (0, 0.1, 0.5, 1.0 T)
(b) FP vs distance, log y-axis — transverse B0 parallel to M (+x, 0, 0.1, 0.5, 1.0 T)
(c) FP at fixed distance (200 um) vs B0 magnitude, axial and transverse
(d) Enhancement ratio FP(B0) / FP(B0=0) vs distance for both directions

Default stent: R=1.5 mm, w=100 um, t=80 um, L=500 um, M=1.0 MA/m, 8 struts.
assume_saturation=True throughout.

Run standalone::

    python -m stent_capture.figures.fig08_force_parameter
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT, COLORS_CODE_CALIBRATED, COLORS_THRESHOLD_HIGH,
    COLORS_MARKER_REFERENCE
)
from stent_capture.physics.external_field import TotalField, UniformExternalField


_B0_CASES = [0.0, 0.1, 0.5, 1.0, 1.5]


def _force_parameter_profile(
    B0_vec: np.ndarray,
    d: np.ndarray,
) -> np.ndarray:
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    B0_magnitude = np.linalg.norm(B0_vec)
    ring = make_ring(B0_magnitude=B0_magnitude)
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



_COLORS = [
    COLORS_MARKER_REFERENCE,
    COLORS_CODE_DEFAULT,
    COLORS_CODE_CALIBRATED,
    COLORS_THRESHOLD_HIGH,
    COLORS_MARKER_REFERENCE,
]
_LINESTYLES_AX    = ["-",  "-",  "-",  "-",  "-"]
_LINESTYLES_TRANS = ["--", "--", "--", "--", "--"]



def make_figure():
    d_profile = np.linspace(5e-6, 1.5e-3, 200)
    d_um = d_profile * 1e6

    fp_axial  = {}
    fp_trans  = {}

    for B0 in _B0_CASES:
        fp_axial[B0] = _force_parameter_profile(
            np.array([0.0, 0.0, B0]), d_profile
        )
        fp_trans[B0] = _force_parameter_profile(
            np.array([B0, 0.0, 0.0]), d_profile
        )

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(14, 11))
    gs = GridSpec(2, 4, figure=fig)
    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:4])
    ax_c = fig.add_subplot(gs[1, 1:3])

    for B0, col in zip(_B0_CASES, _COLORS):
        label = f"B0 = {B0:.1f} T (axial)"
        ax_a.semilogy(d_um, fp_axial[B0], color=col, lw=2.0, label=label)
    ax_a.set_xlabel("Distance from stent surface (um)")
    ax_a.set_ylabel("|B_total| · |∇|B_total|| (T²/m)")
    ax_a.set_title("(a)")
    ax_a.set_xlim(0, 1500)
    ax_a.set_ylim(1e-2, 1e4)
    ax_a.legend(fontsize=10)

    for B0, col in zip(_B0_CASES, _COLORS):
        label = f"B0 = {B0:.1f} T (transverse, +x)"
        ax_b.semilogy(d_um, fp_trans[B0], color=col, lw=2.0, ls="--", label=label)
    ax_b.set_xlabel("Distance from stent surface (um)")
    ax_b.set_ylabel("|B_total| · |∇|B_total|| (T²/m)")
    ax_b.set_title("(b)")
    ax_b.set_xlim(0, 1500)
    ax_b.set_ylim(1e-2, 1e4)
    ax_b.legend(fontsize=10)

    fp0 = fp_axial[0.0]
    fp0_safe = np.where(fp0 > 0, fp0, np.nan)
    for B0, col in zip(_B0_CASES[1:], _COLORS[1:]):
        ratio_ax = fp_axial[B0] / fp0_safe
        ratio_tr = fp_trans[B0] / fp0_safe
        ax_c.semilogy(d_um, ratio_ax, color=col, lw=2.0, ls="-",
                      label=f"B0={B0:.1f}T axial")
        ax_c.semilogy(d_um, ratio_tr, color=col, lw=2.0, ls="--",
                      label=f"B0={B0:.1f}T transverse")
    ax_c.axhline(1.0, color="k", ls=":", lw=1.2, alpha=0.5, label="No enhancement")
    ax_c.set_xlabel("Distance from stent surface (um)")
    ax_c.set_ylabel("FP(B0) / FP(B0 = 0)")
    ax_c.set_title("(c)")
    ax_c.set_xlim(0, 1500)
    ax_c.legend(fontsize=9, ncol=2)

    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 8: Force parameter...")
    fig = make_figure()
    fig.savefig(OUT / "fig8_force_parameter.png")
    fig.savefig(OUT / "fig8_force_parameter.pdf")
    plt.close(fig)
    print("  [OK] fig8_force_parameter saved")


if __name__ == "__main__":
    main()
