"""
Fig 6 — Effect of number of struts on gradient and angular uniformity.

Run standalone::

    python -m stent_capture.figures.fig06_n_struts
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import (
    DEFAULTS, THRESHOLDS, OUT, make_ring,
    COMSOL_CROSSINGS,
)
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT, COLORS_THRESHOLD, COLORS_THRESHOLD_MEDIUM, COLORS_THRESHOLD_HIGH,
    COLORS_MARKER_REFERENCE
)


def make_figure():
    # Unified threshold colors (green → orange → red for low → medium → high)
    threshold_colors = {
        "40 T/m": COLORS_THRESHOLD,
        "100 T/m": COLORS_THRESHOLD_MEDIUM,
        "300 T/m": COLORS_THRESHOLD_HIGH,
    }

    # Match COMSOL variants: V1 (6 cells), V2-2C (12 cells, headline), V3-2C (18 cells)
    n_vals  = [6, 12, 18]
    labels  = ["V1 (6 cells)", "V2-2C (12 cells) — headline", "V3-2C (18 cells)"]
    # Unified colors: gray (reference), blue (standard), gray (variant)
    colors  = [COLORS_MARKER_REFERENCE, COLORS_CODE_DEFAULT, COLORS_MARKER_REFERENCE]

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    d = np.linspace(5e-6, 1e-3, 200)
    z = np.zeros_like(d)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    for ns, col, lbl in zip(n_vals, colors, labels):
        ring = make_ring(n_struts=ns)
        G = ring.grad_B(d + r_outer, np.zeros_like(d), z)
        ax.semilogy(d * 1e6, G, color=col, lw=1.8, label=lbl)
    for lbl, val in THRESHOLDS.items():
        c = threshold_colors[lbl]
        ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.5)

    # COMSOL V2-2C reference threshold crossing points (circles instead of stars)
    for lbl, d_cross in COMSOL_CROSSINGS.items():
        th_val = THRESHOLDS[lbl]
        c = threshold_colors[lbl]
        ax.scatter([d_cross * 1e6], [th_val],
                   marker="o", s=100, color=c, zorder=6, edgecolors="black", lw=0.8)
        ax.axvline(d_cross * 1e6, color=c, ls="--", lw=0.8, alpha=0.35)

    ax.set_xlabel("Distance from stent surface (um)")
    ax.set_ylabel("|grad_B (T/m)")
    ax.set_title("(a) Gradient through strut — COMSOL variants\n(circles = COMSOL V2-2C reference crossings)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 1000)
    ax.set_ylim(0.1, 1e4)

    # Angular variation at fixed distance
    fixed_d = 200e-6
    r_eval = r_outer + fixed_d
    n_angles = 200
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)

    ax = axes[1]
    for ns, col, lbl in zip(n_vals, colors, labels):
        ring = make_ring(n_struts=ns)
        obs_x = r_eval * np.cos(angles)
        obs_y = r_eval * np.sin(angles)
        obs_z = np.zeros_like(angles)
        G = ring.grad_B(obs_x, obs_y, obs_z)
        ax.plot(np.degrees(angles), G, color=col, lw=1.5, label=lbl)
    ax.set_xlabel("Angle (degrees)")
    ax.set_ylabel("|grad_B (T/m)")
    ax.set_title(f"(b) Angular variation at {fixed_d*1e6:.0f} um from surface")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 360)

    fig.suptitle(
        "Effect of Cell Count on Gradient Distribution — COMSOL variants\n"
        "V1 (6), V2-2C (12, headline), V3-2C (18)  |  Analytical 3-D model",
        fontsize=12, y=0.98,
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 6: Number of struts study…")
    fig = make_figure()
    fig.savefig(OUT / "fig6_n_struts.png")
    fig.savefig(OUT / "fig6_n_struts.pdf")
    plt.close(fig)
    print("  [OK] fig6_n_struts saved")


if __name__ == "__main__":
    main()
