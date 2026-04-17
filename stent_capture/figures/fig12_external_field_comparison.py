"""
Fig 12 — External field comparison: effect of B0 on gradient magnitude.

2x2 panel figure:

(a) Gradient through-strut profile — B0 = 0  (stent field alone)
(b) Gradient through-strut profile — B0 = 0.5 T axial (+z)
(c) Gradient through-strut profile — B0 = 0.5 T transverse (+x, parallel to M)
(d) Capture distance vs B0 magnitude (0 to 1 T, axial direction) based on
    the gradient threshold |∇|B_total|| ≥ threshold

Physics note — WHY |∇|B_total|| is REDUCED by axial B0
-------------------------------------------------------
When B0 is axial (+z) and the stent magnetisation is radial (+x), the total
field at an observation point near the strut is dominated by B0:

    B_total  ~  [B_stent_x, 0, B0_z]   (B_stent_x << B0_z at large B0)

The gradient of the magnitude is:

    ∂|B_total|/∂x = (B_total · ∂B_stent/∂x) / |B_total|
                   ~  (B_stent_x / B0_z) * ∂B_stent_x/∂x  + small terms

Because B_total is rotated towards z, the projection of the radial gradient
∂B_stent/∂x onto B_total is suppressed by the factor B_stent_x / B0_z ~ 0.06.
Hence |∇|B_total|| is correctly reduced — this is NOT a numerical error.

However, this does NOT mean B0 reduces cell capture.  The relevant metric for
the force on a SPION-labelled cell is the force parameter FP = |B_total|*|∇|B_total||
(see fig13).  When B0 = 0.5 T, |B_total| increases ~17x while |∇|B_total||
decreases ~5x, giving a net ~3x gain in FP at 200 um and ~55x at 500 um.

This figure shows |∇|B_total|| to document the physics correctly.
Panel (d) computes capture distance from gradient-threshold crossings; it is
provided for completeness but should be interpreted alongside fig13 which
shows the physically relevant force parameter.

Default stent: R=1.5 mm, w=100 um, t=80 um, L=500 um, M=1.0 MA/m, 8 struts.
assume_saturation=True throughout.

Run standalone::

    python -m stent_capture.figures.fig12_external_field_comparison
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT,
    COLORS_THRESHOLD, COLORS_THRESHOLD_MEDIUM, COLORS_THRESHOLD_HIGH,
    COLORS_MARKER_REFERENCE
)
from stent_capture.physics.external_field import TotalField, UniformExternalField


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _radial_profile(
    B0_vec: np.ndarray,
    d: np.ndarray,
) -> np.ndarray:
    """
    Compute |∇|B_total|| along the through-strut (+x) radial line at z=0.

    Parameters
    ----------
    B0_vec : (3,) array — external field vector in Tesla
    d      : 1-D array — distances from stent outer surface (m)

    Returns
    -------
    G : 1-D array — |∇|B_total|| in T/m
    """
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    # Adaptive M: at B0 = 1.5 T use COMSOL-calibrated M
    B0_magnitude = np.linalg.norm(B0_vec)
    ring = make_ring(B0_magnitude=B0_magnitude)
    ring.assume_saturation = True

    ext = UniformExternalField(B0_vec) if np.any(B0_vec != 0) else None
    tf  = TotalField(ring, ext)

    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    obs_z = np.zeros_like(d)

    return tf.grad_B(obs_x, obs_y, obs_z)


def _capture_distance_vs_B0(
    B0_direction: np.ndarray,
    B0_magnitudes: np.ndarray,
    d: np.ndarray,
    threshold: float,
) -> np.ndarray:
    """
    Compute capture distance at `threshold` as a function of B0 magnitude.

    Parameters
    ----------
    B0_direction : unit vector (3,) for the direction of B0
    B0_magnitudes : 1-D array of B0 magnitudes in Tesla
    d             : 1-D distance array from stent surface (m)
    threshold     : gradient threshold in T/m

    Returns
    -------
    cap_dist_um : 1-D array — capture distances in um, same length as B0_magnitudes
    """
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    obs_z = np.zeros_like(d)
    d_um  = d * 1e6

    cap_dists = []
    for B0_mag in B0_magnitudes:
        B0_vec = B0_mag * B0_direction
        G = _radial_profile(B0_vec, d)
        above = np.where(G >= threshold)[0]
        cap_dists.append(float(d_um[above[-1]]) if len(above) > 0 else 0.0)

    return np.array(cap_dists)


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    # Radial profiles
    d_profile = np.linspace(5e-6, 1.5e-3, 200)
    d_um = d_profile * 1e6

    G_noB0  = _radial_profile(np.zeros(3),          d_profile)   # panel (a)
    G_axial = _radial_profile(np.array([0, 0, 1.5]), d_profile)  # panel (b) — 1.5 T MRI
    G_trans = _radial_profile(np.array([1.5, 0, 0]), d_profile)  # panel (c) — 1.5 T MRI

    # Capture distance vs B0 magnitude — panel (d), extended to MRI strength
    B0_range = np.linspace(0.0, 1.5, 31)
    cap_dist_data: dict[str, np.ndarray] = {}

    for lbl, thr in THRESHOLDS.items():
        cap_dist_data[lbl] = _capture_distance_vs_B0(
            np.array([0.0, 0.0, 1.0]),
            B0_range,
            np.linspace(5e-6, 1.5e-3, 200),
            thr,
        )

    # Unified threshold colors (green → orange → red for low → medium → high)
    threshold_colors = {
        "40 T/m": COLORS_THRESHOLD,
        "100 T/m": COLORS_THRESHOLD_MEDIUM,
        "300 T/m": COLORS_THRESHOLD_HIGH,
    }

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # Annotation explaining the physics of gradient reduction
    _physics_note = (
        "Note: axial B0 rotates B_total towards z,\n"
        "suppressing the projection of dB_stent/dx\n"
        "onto B_total. Gradient reduction is correct\n"
        "physics, NOT a numerical error.\n"
        "See fig13 for the force parameter |B|*|∇B|."
    )

    def _profile_ax(ax, G, title, panel_label, add_note=False):
        ax.semilogy(d_um, G_noB0, color=COLORS_MARKER_REFERENCE, ls="--", lw=1.5, alpha=0.7,
                   label="B0 = 0 (reference)")
        ax.semilogy(d_um, G, color=COLORS_CODE_DEFAULT, lw=2.0, label=title)
        for lbl, val in THRESHOLDS.items():
            c = threshold_colors[lbl]
            ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.7, label=lbl)
        ax.set_xlabel("Distance from stent surface (um)")
        ax.set_ylabel("|∇|B_total|| (T/m)")
        ax.set_title(f"({panel_label}) {title}")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 1500)
        ax.set_ylim(0.1, 1e4)
        if add_note:
            ax.text(
                0.97, 0.97, _physics_note,
                transform=ax.transAxes,
                ha="right", va="top",
                fontsize=7, color="#555555",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8),
            )

    _profile_ax(ax_a, G_noB0,  "B0 = 0 (stent field only)", "a")
    _profile_ax(ax_b, G_axial, "B0 = 1.5 T axial (+z, MRI)", "b", add_note=True)
    _profile_ax(ax_c, G_trans, "B0 = 1.5 T transverse (+x, parallel to M)", "c")

    # Panel (d): capture distance vs B0 — note in subtitle that grad metric
    # underestimates capture benefit; refer reader to fig13
    ax_d.set_title(
        "(d) Capture distance vs B0 magnitude (axial field)\n"
        r"[threshold on $|\nabla|B_{total}||$ — see fig13 for force parameter]",
        fontsize=10,
    )
    for lbl, dists in cap_dist_data.items():
        c = threshold_colors[lbl]
        ax_d.plot(B0_range, dists, "o-", color=c, lw=2, ms=5, label=lbl)
    ax_d.axvline(1.5, color=COLORS_MARKER_REFERENCE, ls="--", lw=1.5, alpha=0.8,
                label="1.5 T (MRI / COMSOL)")
    ax_d.set_xlabel("B0 magnitude (T)")
    ax_d.set_ylabel("Capture distance (um)")
    ax_d.legend(fontsize=9)
    ax_d.set_xlim(0, 1.5)
    ax_d.grid(True, alpha=0.3)

    fig.suptitle(
        "Effect of Uniform External Field on Gradient Magnitude |∇|B_total||\n"
        "(12 struts / V2-2C, R = 1.5 mm, M = 1.0 MA/m, assume_saturation = True,\n"
        "through-strut radial profile at z = 0 — for force parameter see fig13)",
        fontsize=12, y=0.98,
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 12: External field comparison...")
    fig = make_figure()
    fig.savefig(OUT / "fig12_external_field_comparison.png")
    fig.savefig(OUT / "fig12_external_field_comparison.pdf")
    plt.close(fig)
    print("  [OK] fig12_external_field_comparison saved")


if __name__ == "__main__":
    main()
