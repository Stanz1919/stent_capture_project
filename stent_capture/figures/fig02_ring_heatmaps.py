"""
Fig 2 — Stent ring cross-section: |B| and |grad_B heatmaps (3-D at z = 0).

Run standalone::

    python -m stent_capture.figures.fig02_ring_heatmaps
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.figures.style import COLORS_CONCENTRATION


def make_figure():
    ring = make_ring()
    ext = 3.5e-3
    n = 200
    x = np.linspace(-ext, ext, n)
    y = np.linspace(-ext, ext, n)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    Bmag = ring.B_magnitude(X, Y, Z)
    Gmag = ring.grad_B(X, Y, Z)

    # Mask strut interiors
    for i in range(ring.n_struts):
        dist = np.sqrt((X - ring.cx[i])**2 + (Y - ring.cy[i])**2)
        mask = dist < max(ring.w, ring.t) * 1.2
        Bmag[mask] = np.nan
        Gmag[mask] = np.nan

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    im = ax.pcolormesh(X * 1e3, Y * 1e3, Bmag * 1000,
                       norm=LogNorm(vmin=0.01, vmax=100),
                       cmap=COLORS_CONCENTRATION, shading="auto")
    fig.colorbar(im, ax=ax, label="|B| (mT)", shrink=0.85)
    ax.add_patch(Circle((0, 0), DEFAULTS["R"] * 1e3, fill=False, color="w", ls="--", lw=1))
    for i in range(ring.n_struts):
        ax.plot(ring.cx[i] * 1e3, ring.cy[i] * 1e3, "ws", ms=4)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("(a) Magnetic flux density |B|")
    ax.set_aspect("equal")

    ax = axes[1]
    im = ax.pcolormesh(X * 1e3, Y * 1e3, Gmag,
                       norm=LogNorm(vmin=1, vmax=5000),
                       cmap=COLORS_CONCENTRATION, shading="auto")
    fig.colorbar(im, ax=ax, label="|grad_B (T/m)", shrink=0.85)
    ax.add_patch(Circle((0, 0), DEFAULTS["R"] * 1e3, fill=False, color="w", ls="--", lw=1))
    for i in range(ring.n_struts):
        ax.plot(ring.cx[i] * 1e3, ring.cy[i] * 1e3, "ws", ms=4)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("(b) Field gradient |grad_B")
    ax.set_aspect("equal")

    fig.suptitle(
        f"Stent Ring Cross-Section — 3-D model at z = 0\n"
        f"({DEFAULTS['n_struts']} struts, R = {DEFAULTS['R']*1e3:.1f} mm, "
        f"M = {DEFAULTS['M']/1e6:.1f} MA/m)",
        fontsize=13, y=0.98,
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 2: Ring cross-section heatmaps…")
    fig = make_figure()
    fig.savefig(OUT / "fig2_ring_heatmaps.png")
    fig.savefig(OUT / "fig2_ring_heatmaps.pdf")
    plt.close(fig)
    print("  [OK] fig2_ring_heatmaps saved")


if __name__ == "__main__":
    main()
