"""
Fig 6 — Effect of number of struts on gradient and angular uniformity.

Run standalone::

    python -m stent_capture.figures.fig06_n_struts
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring


def make_figure():
    n_vals = [4, 6, 8, 10, 12]
    colors = plt.cm.tab10(np.linspace(0, 0.5, len(n_vals)))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    d = np.linspace(5e-6, 1e-3, 200)
    z = np.zeros_like(d)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    for ns, col in zip(n_vals, colors):
        ring = make_ring(n_struts=ns)
        G = ring.grad_B(d + r_outer, np.zeros_like(d), z)
        ax.semilogy(d * 1e6, G, color=col, lw=1.8, label=f"{ns} struts")
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.5)
    ax.set_xlabel("Distance from stent surface (µm)")
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title("(a) Gradient through strut")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 1000)
    ax.set_ylim(0.1, 1e4)

    # Angular variation at fixed distance
    fixed_d = 200e-6
    r_eval = r_outer + fixed_d
    n_angles = 200
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)

    ax = axes[1]
    for ns, col in zip(n_vals, colors):
        ring = make_ring(n_struts=ns)
        obs_x = r_eval * np.cos(angles)
        obs_y = r_eval * np.sin(angles)
        obs_z = np.zeros_like(angles)
        G = ring.grad_B(obs_x, obs_y, obs_z)
        ax.plot(np.degrees(angles), G, color=col, lw=1.5, label=f"{ns} struts")
    ax.set_xlabel("Angle (degrees)")
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title(f"(b) Angular variation at {fixed_d*1e6:.0f} µm from surface")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 360)

    fig.suptitle("Effect of Number of Struts on Gradient Distribution — 3-D model",
                 fontsize=13, y=1.02)
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
