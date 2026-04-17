"""
Fig 10 — 3-D model convergence: gradient vs radial distance for increasing L.

As strut length L → ∞ the 3-D Akoun & Yonnet model converges to the
infinite-cylinder (2-D) limit.  This figure demonstrates that convergence
by sweeping L from 200 um to 100 mm and showing the midplane (z = 0) radial
gradient profile.  The L = 100 mm curve serves as the effective infinite-length
reference.

In the original script this figure compared the 2-D surface-charge model
against the 3-D model; that comparison is no longer meaningful now that the
package uses the 3-D model exclusively.

Run standalone::

    python -m stent_capture.figures.fig10_convergence
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_THRESHOLD, COLORS_THRESHOLD_MEDIUM, COLORS_THRESHOLD_HIGH,
    COLORS_MARKER_REFERENCE
)


def make_figure():
    # Unified threshold colors (green → orange → red for low → medium → high)
    threshold_colors = {
        "40 T/m": COLORS_THRESHOLD,
        "100 T/m": COLORS_THRESHOLD_MEDIUM,
        "300 T/m": COLORS_THRESHOLD_HIGH,
    }

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    d = np.linspace(5e-6, 1.5e-3, 250)
    r_outer = R + t / 2
    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    obs_z = np.zeros_like(d)
    d_um = d * 1e6

    # L = 100 mm acts as the practical infinite-length limit
    ring_ref = make_ring(L=100e-3)
    G_ref = ring_ref.grad_B(obs_x, obs_y, obs_z)

    L_vals = [200e-6, 500e-6, 1e-3, 5e-3, 20e-3]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(L_vals)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for ax, yscale in zip(axes, ["log", "linear"]):
        fn = ax.semilogy if yscale == "log" else ax.plot
        fn(d_um, G_ref, color=COLORS_MARKER_REFERENCE, lw=2.5,
           label="L = 100 mm  (infinite-length limit)", zorder=5)
        for L, col in zip(L_vals, colors):
            ring = make_ring(L=L)
            G = ring.grad_B(obs_x, obs_y, obs_z)
            fn(d_um, G, "--", color=col, lw=1.8,
               label=f"L = {L*1e3:.1f} mm")
        for _, val in THRESHOLDS.items():
            lbl = list(THRESHOLDS.keys())[list(THRESHOLDS.values()).index(val)]
            c = threshold_colors[lbl]
            ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.5)
        ax.set_xlabel("Distance from stent surface (um)")
        ax.set_ylabel("|grad_B (T/m)")
        ax.set_xlim(0, 1500)
        ax.legend(fontsize=8)

    axes[0].set_title("(a) Log scale")
    axes[0].set_ylim(0.1, 1e4)
    axes[1].set_title("(b) Linear scale")
    axes[1].set_ylim(0, 2500)

    fig.suptitle(
        "3-D Model Convergence: Effect of Strut Length on Midplane Gradient\n"
        "(z = 0, through-strut radial profile)",
        fontsize=13, y=0.98,
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 10: 3-D convergence check…")
    fig = make_figure()
    fig.savefig(OUT / "fig10_convergence.png")
    fig.savefig(OUT / "fig10_convergence.pdf")
    plt.close(fig)
    print("  [OK] fig10_convergence saved")


if __name__ == "__main__":
    main()
