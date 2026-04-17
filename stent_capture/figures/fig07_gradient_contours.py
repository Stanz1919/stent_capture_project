"""
Fig 7 — Cross-section with gradient contours at capture thresholds.

Run standalone::

    python -m stent_capture.figures.fig07_gradient_contours
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CONCENTRATION,
    COLORS_THRESHOLD, COLORS_THRESHOLD_MEDIUM, COLORS_THRESHOLD_HIGH
)


def make_figure():
    ring = make_ring()
    ext = 3.0e-3
    n = 200
    x = np.linspace(-ext, ext, n)
    y = np.linspace(-ext, ext, n)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    Gmag = ring.grad_B(X, Y, Z)

    for i in range(ring.n_struts):
        dist = np.sqrt((X - ring.cx[i])**2 + (Y - ring.cy[i])**2)
        Gmag[dist < max(ring.w, ring.t) * 1.2] = np.nan

    fig, ax = plt.subplots(figsize=(8, 7))

    im = ax.pcolormesh(X * 1e3, Y * 1e3, Gmag,
                       norm=LogNorm(vmin=5, vmax=5000),
                       cmap=COLORS_CONCENTRATION, shading="auto")

    # Unified threshold colors (green → orange → red for low → medium → high)
    threshold_colors = {
        "40 T/m": COLORS_THRESHOLD,
        "100 T/m": COLORS_THRESHOLD_MEDIUM,
        "300 T/m": COLORS_THRESHOLD_HIGH,
    }

    Gmag_filled = np.nan_to_num(Gmag, nan=0)
    for lbl, val in THRESHOLDS.items():
        c = threshold_colors[lbl]
        cs = ax.contour(X * 1e3, Y * 1e3, Gmag_filled,
                        levels=[val], colors=[c], linewidths=2)
        if cs.allsegs[0]:
            ax.clabel(cs, fmt=f"{val} T/m", fontsize=9, colors=[c])

    ax.add_patch(Circle((0, 0), DEFAULTS["R"] * 1e3,
                         fill=False, color="w", ls="--", lw=1.5))
    for i in range(ring.n_struts):
        ax.plot(ring.cx[i] * 1e3, ring.cy[i] * 1e3, "ws", ms=5)

    fig.colorbar(im, ax=ax, label="|grad_B (T/m)", shrink=0.85)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(
        f"Gradient Cross-Section with Capture Thresholds — 3-D model\n"
        f"({DEFAULTS['n_struts']} struts, R = {DEFAULTS['R']*1e3:.1f} mm, "
        f"M = {DEFAULTS['M']/1e6:.1f} MA/m)"
    )
    ax.set_aspect("equal")
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 7: Gradient heatmap with contours…")
    fig = make_figure()
    fig.savefig(OUT / "fig7_gradient_contours.png")
    fig.savefig(OUT / "fig7_gradient_contours.pdf")
    plt.close(fig)
    print("  [OK] fig7_gradient_contours saved")


if __name__ == "__main__":
    main()
