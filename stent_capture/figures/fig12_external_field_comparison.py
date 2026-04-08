"""
Fig 12 — External field comparison: effect of B0 on gradient and capture distance.

2×2 panel figure:

(a) Gradient through-strut profile — B0 = 0  (stent field alone)
(b) Gradient through-strut profile — B0 = 0.5 T axial (+z)
(c) Gradient through-strut profile — B0 = 0.5 T transverse (+x, parallel to M)
(d) Capture distance vs B0 magnitude (0 to 1 T) for three threshold levels

Physics note
------------
Axial B0 (+z) is perpendicular to the radial strut magnetisation (+x), so it
maximally changes the direction of B_total and tends to INCREASE the gradient
of |B_total| near the strut.  Transverse B0 (+x) is parallel to M; it raises
|B_total| nearly uniformly in the radial direction, which can REDUCE the
relative spatial variation and hence the gradient of the magnitude.

The capture distance for each threshold is the outermost radial position (from
the stent outer surface) at which |∇|B_total|| ≥ threshold.

Default stent: R=1.5 mm, w=100 µm, t=80 µm, L=500 µm, M=1.0 MA/m, 8 struts.
The stent is assumed saturated (M = M_sat) at all B0 values shown here.
StentRing is constructed once per B0 magnitude with assume_saturation=True
to signal this modelling intent clearly.

Run standalone::

    python -m stent_capture.figures.fig12_external_field_comparison
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring
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

    ring = make_ring()
    ring.assume_saturation = True   # document saturation assumption in output

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
    cap_dist_um : 1-D array — capture distances in µm, same length as B0_magnitudes
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

    G_noB0  = _radial_profile(np.zeros(3),         d_profile)   # panel (a)
    G_axial = _radial_profile(np.array([0,0,0.5]),  d_profile)   # panel (b)
    G_trans = _radial_profile(np.array([0.5,0,0]),  d_profile)   # panel (c)

    # Capture distance vs B0 magnitude — panel (d)
    B0_range = np.linspace(0.0, 1.0, 21)
    cap_dist_data: dict[str, np.ndarray] = {}

    for lbl, thr in THRESHOLDS.items():
        # Use axial direction (most relevant for clinical MRI)
        cap_dist_data[lbl] = _capture_distance_vs_B0(
            np.array([0.0, 0.0, 1.0]),   # axial +z direction
            B0_range,
            np.linspace(5e-6, 1.5e-3, 200),
            thr,
        )

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # Common profile plot settings
    def _profile_ax(ax, G, title, panel_label):
        ax.semilogy(d_um, G_noB0, "k--", lw=1.5, alpha=0.5, label="B0 = 0 (reference)")
        ax.semilogy(d_um, G,      "b-",  lw=2.0, label=title)
        for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
            ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.7, label=lbl)
        ax.set_xlabel("Distance from stent surface (µm)")
        ax.set_ylabel("|∇|B_total|| (T/m)")
        ax.set_title(f"({panel_label}) {title}")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 1500)
        ax.set_ylim(0.1, 1e4)

    _profile_ax(ax_a, G_noB0,  "B0 = 0 (stent field only)", "a")
    _profile_ax(ax_b, G_axial, "B0 = 0.5 T axial (+z)",     "b")
    _profile_ax(ax_c, G_trans, "B0 = 0.5 T transverse (+x, parallel to M)", "c")

    # Panel (d): capture distance vs B0
    ax_d.set_title("(d) Capture distance vs B0 magnitude (axial field)")
    for (lbl, dists), c in zip(cap_dist_data.items(), TH_COLORS):
        ax_d.plot(B0_range, dists, "o-", color=c, lw=2, ms=5, label=lbl)
    ax_d.set_xlabel("B0 magnitude (T)")
    ax_d.set_ylabel("Capture distance (µm)")
    ax_d.legend(fontsize=9)
    ax_d.set_xlim(0, 1.0)
    ax_d.grid(True, alpha=0.3)

    fig.suptitle(
        "Effect of Uniform External Field on Gradient and Capture Distance\n"
        "(8 struts, R = 1.5 mm, M = 1.0 MA/m, assume_saturation = True,\n"
        "through-strut radial profile at z = 0)",
        fontsize=13, y=1.01,
    )
    plt.tight_layout()
    return fig


def main():
    print("  Fig 12: External field comparison…")
    fig = make_figure()
    fig.savefig(OUT / "fig12_external_field_comparison.png")
    fig.savefig(OUT / "fig12_external_field_comparison.pdf")
    plt.close(fig)
    print("  [OK] fig12_external_field_comparison saved")


if __name__ == "__main__":
    main()
