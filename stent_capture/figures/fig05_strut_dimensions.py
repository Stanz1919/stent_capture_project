"""
Fig 5 — Effect of strut thickness and width on gradient.

Run standalone::

    python -m stent_capture.figures.fig05_strut_dimensions
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring


def make_figure():
    R = DEFAULTS["R"]
    thicknesses = [40e-6, 60e-6, 80e-6, 100e-6, 120e-6]
    widths       = [50e-6, 80e-6, 100e-6, 150e-6, 200e-6]
    colors_t = plt.cm.plasma(np.linspace(0.15, 0.85, len(thicknesses)))
    colors_w = plt.cm.cividis(np.linspace(0.15, 0.85, len(widths)))

    d = np.linspace(5e-6, 1e-3, 200)
    z = np.zeros_like(d)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    for thick, col in zip(thicknesses, colors_t):
        r_outer = R + thick / 2
        ring = make_ring(t=thick)
        G = ring.grad_B(d + r_outer, np.zeros_like(d), z)
        ax.semilogy(d * 1e6, G, color=col, lw=1.8,
                    label=f"t = {thick*1e6:.0f} µm")
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.5)
    ax.set_xlabel("Distance from strut surface (µm)")
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title(f"(a) Varying thickness (w = {DEFAULTS['w']*1e6:.0f} µm)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 1000)
    ax.set_ylim(0.1, 1e4)

    ax = axes[1]
    r_outer = R + DEFAULTS["t"] / 2
    for w, col in zip(widths, colors_w):
        ring = make_ring(w=w)
        G = ring.grad_B(d + r_outer, np.zeros_like(d), z)
        ax.semilogy(d * 1e6, G, color=col, lw=1.8,
                    label=f"w = {w*1e6:.0f} µm")
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.5)
    ax.set_xlabel("Distance from strut surface (µm)")
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title(f"(b) Varying width (t = {DEFAULTS['t']*1e6:.0f} µm)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 1000)
    ax.set_ylim(0.1, 1e4)

    fig.suptitle("Effect of Strut Dimensions on Field Gradient — 3-D model",
                 fontsize=13, y=1.02)
    plt.tight_layout()
    return fig


def main():
    print("  Fig 5: Strut dimensions sweep…")
    fig = make_figure()
    fig.savefig(OUT / "fig5_strut_dimensions.png")
    fig.savefig(OUT / "fig5_strut_dimensions.pdf")
    plt.close(fig)
    print("  [OK] fig5_strut_dimensions saved")


if __name__ == "__main__":
    main()
