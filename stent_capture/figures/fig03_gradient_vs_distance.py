"""
Fig 3 — Gradient vs distance: analytical profile with COMSOL reference.

2-panel figure:

(a) Log-scale through-strut gradient profile (M = 1.0 MA/m, B0 = 0).
    Coloured circles: COMSOL V2-2C threshold crossing distances
                      (300 T/m @ 120 um, 100 T/m @ 240 um, 40 T/m @ 380 um)
                      — from FEM at B0 = 1.5 T (MRI), shown for reference.
    Horizontal coloured lines: threshold levels.

(b) Log-scale between-struts profile (M = 1.0 MA/m, B0 = 0), no COMSOL
    overlay (COMSOL reports through-strut values only).

Run standalone::

    python -m stent_capture.figures.fig03_gradient_vs_distance
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import (
    DEFAULTS, THRESHOLDS, OUT,
    make_ring, COMSOL_CROSSINGS,
)
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT,
    COLORS_THRESHOLD, COLORS_THRESHOLD_MEDIUM, COLORS_THRESHOLD_HIGH,
    COLORS_MARKER_REFERENCE
)


def make_figure():
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    d = np.linspace(5e-6, 2e-3, 300)
    z = np.zeros_like(d)
    d_um = d * 1e6

    # Default M = 1.0 MA/m, B0 = 0
    ring_def = make_ring()
    G_def_through = ring_def.grad_B(d + r_outer, np.zeros_like(d), z)

    # Between-struts profile
    angle_bet = np.pi / DEFAULTS["n_struts"]
    r_bet = d + r_outer
    G_def_between = ring_def.grad_B(
        r_bet * np.cos(angle_bet), r_bet * np.sin(angle_bet), z
    )

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(15, 6))

    # Color mapping for thresholds
    threshold_colors = {
        "40 T/m": COLORS_THRESHOLD,
        "100 T/m": COLORS_THRESHOLD_MEDIUM,
        "300 T/m": COLORS_THRESHOLD_HIGH,
    }

    # -----------------------------------------------------------------------
    # Panel (a): through-strut, log scale — COMSOL comparison
    # -----------------------------------------------------------------------
    ax_a.semilogy(d_um, G_def_through, color=COLORS_CODE_DEFAULT, lw=2.0,
                  label="Analytical model: M = 1.0 MA/m, B0 = 0")

    # Horizontal threshold lines
    for lbl, val in THRESHOLDS.items():
        c = threshold_colors[lbl]
        ax_a.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)

    # COMSOL V2-2C reference crossing points — circles (instead of stars)
    for lbl, d_cross in COMSOL_CROSSINGS.items():
        th_val = THRESHOLDS[lbl]
        d_cross_um = d_cross * 1e6
        c = threshold_colors[lbl]
        # Vertical guide line at COMSOL crossing distance
        ax_a.axvline(d_cross_um, color=c, ls="--", lw=0.9, alpha=0.45)
        # Circle marker at the crossing point
        ax_a.scatter([d_cross_um], [th_val],
                     marker="o", s=120, color=c, zorder=6, edgecolors="black", lw=0.8,
                     label=f"COMSOL: {lbl} @ {d_cross_um:.0f} um")

    ax_a.set_xlabel("Distance from stent surface (um)")
    ax_a.set_ylabel("|grad_B (T/m)")
    ax_a.set_title("(a) Through-strut profile — COMSOL comparison (log scale)")
    ax_a.set_xlim(0, 1000)
    ax_a.set_ylim(0.1, 1e4)
    ax_a.legend(fontsize=8, loc="upper right")

    # Explanatory annotation
    ax_a.text(
        0.03, 0.04,
        "Circles = COMSOL V2-2C reference crossing distances\n"
        "(FEM, mu_r=2, B0=1.5 T — shown for scale)\n"
        "Analytical model uses B0 = 0, M = 1.0 MA/m",
        transform=ax_a.transAxes,
        ha="left", va="bottom", fontsize=7, color="#444444",
        bbox=dict(boxstyle="round,pad=0.35", fc="lightyellow", alpha=0.85),
    )

    # -----------------------------------------------------------------------
    # Panel (b): between-struts, log scale
    # -----------------------------------------------------------------------
    ax_b.semilogy(d_um, G_def_between, color=COLORS_CODE_DEFAULT, lw=2.0,
                  label="Between-struts: M = 1.0 MA/m, B0 = 0")
    ax_b.semilogy(d_um, G_def_through, color=COLORS_CODE_DEFAULT, lw=1.0, alpha=0.25,
                  label="Through-strut (reference)")

    for lbl, val in THRESHOLDS.items():
        c = threshold_colors[lbl]
        ax_b.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)

    ax_b.set_xlabel("Distance from stent surface (um)")
    ax_b.set_ylabel("|grad_B (T/m)")
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
        f"Gradient vs Distance from Stent Surface — Analytical 3-D Model\n"
        f"({DEFAULTS['n_struts']} struts / V2-2C, M = 1.0 MA/m, B0 = 0)\n"
        f"COMSOL reference circles from FEM at B0 = 1.5 T (shown for scale)",
        fontsize=11, y=0.98,
    )
    plt.tight_layout(pad=1.0)
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
