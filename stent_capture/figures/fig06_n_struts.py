"""
Fig 6 — Effect of number of struts on gradient and angular uniformity.

Run standalone::

    python -m stent_capture.figures.fig06_n_struts
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import (
    DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring,
    COMSOL_CROSSINGS,
)


def make_figure():
    # Match COMSOL variants: V1 (6 cells), V2-2C (12 cells), V3-2C (18 cells)
    n_vals  = [6, 12, 18]
    labels  = ["V1 (6 cells)", "V2-2C (12 cells) — headline", "V3-2C (18 cells)"]
    colors  = ["#e74c3c", "#2980b9", "#27ae60"]

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
    for (lbl_t, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.5)

    # COMSOL V2-2C reference threshold crossing points (stars)
    for (lbl, d_cross), c in zip(COMSOL_CROSSINGS.items(), TH_COLORS):
        th_val = THRESHOLDS[lbl]
        ax.scatter([d_cross * 1e6], [th_val],
                   marker="*", s=160, color=c, zorder=6, edgecolors="k", lw=0.4)
        ax.axvline(d_cross * 1e6, color=c, ls="--", lw=0.8, alpha=0.35)

    ax.set_xlabel("Distance from stent surface (µm)")
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title("(a) Gradient through strut — COMSOL variants\n(stars = COMSOL V2-2C reference crossings)")
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
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title(f"(b) Angular variation at {fixed_d*1e6:.0f} µm from surface")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 360)

    fig.suptitle(
        "Effect of Cell Count on Gradient Distribution — COMSOL variants\n"
        "V1 (6), V2-2C (12, headline), V3-2C (18)  |  Analytical 3-D model",
        fontsize=12, y=1.02,
    )
    plt.tight_layout()
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
