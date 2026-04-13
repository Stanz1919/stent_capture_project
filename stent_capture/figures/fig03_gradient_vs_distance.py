"""
Fig 3 — Gradient vs distance: code comparison with COMSOL reference.

2-panel figure:

(a) Log-scale through-strut gradient profile.
    Blue  solid : code default (M = 1.0 MA/m, B0 = 0).
    Orange solid: code calibrated (M = 0.619 MA/m, B0 = 0) — calibrated to
                  match COMSOL's 100 T/m crossing at 240 µm.
    Coloured stars: COMSOL V2-2C reference threshold crossings
                    (300 T/m @ 120 µm, 100 T/m @ 240 µm, 40 T/m @ 380 µm).
    Horizontal coloured lines: threshold levels.

(b) Log-scale between-struts profile, same M comparison, no COMSOL overlay
    (COMSOL reports through-strut values only).

Physics note:
  COMSOL models a soft ferromagnet (mu_r = 2) that concentrates B0 axially
  inside the stent wall, producing dB_concentrated/dr with NO suppression.
  The code models a permanently-magnetised ring (radial M) superposed with a
  separate uniform B0.  At B0 = 0 both approaches produce a radial near-
  surface gradient; calibrating M gives good agreement.  Adding B0 = 1.5 T
  suppresses |nabla|B_total|| in the code (B_total rotates towards z), so the
  two models diverge — this is documented in fig12 and fig13.

Run standalone::

    python -m stent_capture.figures.fig03_gradient_vs_distance
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import (
    DEFAULTS, THRESHOLDS, TH_COLORS, OUT,
    make_ring, make_ring_comsol,
    M_COMSOL_EFF, COMSOL_CROSSINGS,
)


def make_figure():
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    d = np.linspace(5e-6, 2e-3, 300)
    z = np.zeros_like(d)
    d_um = d * 1e6

    # Default M = 1.0 MA/m
    ring_def = make_ring()
    G_def_through = ring_def.grad_B(d + r_outer, np.zeros_like(d), z)

    # Calibrated M = 0.619 MA/m (COMSOL-matched, B0 = 0)
    ring_cal = make_ring_comsol()
    G_cal_through = ring_cal.grad_B(d + r_outer, np.zeros_like(d), z)

    # Between-struts profiles
    angle_bet = np.pi / DEFAULTS["n_struts"]
    r_bet = d + r_outer
    G_def_between = ring_def.grad_B(
        r_bet * np.cos(angle_bet), r_bet * np.sin(angle_bet), z
    )
    G_cal_between = ring_cal.grad_B(
        r_bet * np.cos(angle_bet), r_bet * np.sin(angle_bet), z
    )

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(15, 6))

    # -----------------------------------------------------------------------
    # Panel (a): through-strut, log scale — COMSOL comparison
    # -----------------------------------------------------------------------
    ax_a.semilogy(d_um, G_def_through, "b-", lw=2.0,
                  label=f"Code: M = 1.0 MA/m, B0 = 0")
    ax_a.semilogy(d_um, G_cal_through, color="darkorange", lw=2.0,
                  label=f"Code calibrated: M = {M_COMSOL_EFF/1e6:.3f} MA/m, B0 = 0")

    # Horizontal threshold lines
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax_a.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7)

    # COMSOL V2-2C reference crossing points — coloured stars
    for (lbl, d_cross), c in zip(COMSOL_CROSSINGS.items(), TH_COLORS):
        th_val = THRESHOLDS[lbl]
        d_cross_um = d_cross * 1e6
        # Vertical guide line at COMSOL crossing distance
        ax_a.axvline(d_cross_um, color=c, ls="--", lw=0.9, alpha=0.45)
        # Star marker at the crossing point
        ax_a.scatter([d_cross_um], [th_val],
                     marker="*", s=160, color=c, zorder=6, edgecolors="k", lw=0.4,
                     label=f"COMSOL V2-2C: {lbl} @ {d_cross_um:.0f} µm")

    ax_a.set_xlabel("Distance from stent surface (µm)")
    ax_a.set_ylabel("|∇|B|| (T/m)")
    ax_a.set_title("(a) Through-strut profile — COMSOL comparison (log scale)")
    ax_a.set_xlim(0, 1000)
    ax_a.set_ylim(0.1, 1e4)
    ax_a.legend(fontsize=8, loc="upper right")

    # Explanatory annotation
    ax_a.text(
        0.03, 0.04,
        "Stars = COMSOL V2-2C (FEM, mu_r=2, B0=1.5T)\n"
        "Orange = code calibrated to COMSOL 100 T/m crossing\n"
        "Blue   = code default (M = 1.0 MA/m)\n"
        "Both code curves use B0 = 0 (no suppression)",
        transform=ax_a.transAxes,
        ha="left", va="bottom", fontsize=7, color="#444444",
        bbox=dict(boxstyle="round,pad=0.35", fc="lightyellow", alpha=0.85),
    )

    # -----------------------------------------------------------------------
    # Panel (b): between-struts, log scale
    # -----------------------------------------------------------------------
    ax_b.semilogy(d_um, G_def_between, "b-", lw=2.0,
                  label="Code: M = 1.0 MA/m, B0 = 0")
    ax_b.semilogy(d_um, G_cal_between, color="darkorange", lw=2.0,
                  label=f"Code calibrated: M = {M_COMSOL_EFF/1e6:.3f} MA/m, B0 = 0")
    ax_b.semilogy(d_um, G_def_through, "b-", lw=1.0, alpha=0.25,
                  label="Through-strut (M=1.0, ref.)")

    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax_b.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)

    ax_b.set_xlabel("Distance from stent surface (µm)")
    ax_b.set_ylabel("|∇|B|| (T/m)")
    ax_b.set_title("(b) Between-struts profile (log scale)")
    ax_b.set_xlim(0, 1000)
    ax_b.set_ylim(0.1, 1e4)
    ax_b.legend(fontsize=8, loc="upper right")

    ax_b.text(
        0.97, 0.97,
        "Between-strut gradient is weaker;\nno COMSOL reference available.",
        transform=ax_b.transAxes,
        ha="right", va="top", fontsize=7, color="#555555",
        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8),
    )

    fig.suptitle(
        f"Gradient vs Distance from Stent Surface — Analytical 3-D Model vs COMSOL\n"
        f"({DEFAULTS['n_struts']} struts / V2-2C, B0 = 0 for code curves; "
        f"COMSOL: FEM, mu_r = 2, B0 = 1.5 T)\n"
        f"Calibrated M = {M_COMSOL_EFF/1e6:.3f} MA/m matches COMSOL 100 T/m crossing exactly",
        fontsize=11, y=1.02,
    )
    plt.tight_layout()
    return fig


def main():
    print("  Fig 3: Gradient vs distance (COMSOL comparison)...")
    fig = make_figure()
    fig.savefig(OUT / "fig3_gradient_vs_distance.png")
    fig.savefig(OUT / "fig3_gradient_vs_distance.pdf")
    plt.close(fig)
    print("  [OK] fig3_gradient_vs_distance saved")


if __name__ == "__main__":
    main()
