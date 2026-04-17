"""
Fig 11 — Gradient in the r-z (axial cross-section) plane.

Reveals the axial extent of the capture zone and end-effects at strut tips.
Only possible with the 3-D model.

Run standalone::

    python -m stent_capture.figures.fig11_rz_heatmap
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CONCENTRATION,
    COLORS_THRESHOLD, COLORS_THRESHOLD_MEDIUM, COLORS_THRESHOLD_HIGH
)


def make_figure():
    ring = make_ring()
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    L = DEFAULTS["L"]

    r_vals = np.linspace(R - t, R + 3e-3, 120)
    z_vals = np.linspace(-3e-3, 3e-3, 200)
    R_grid, Z_grid = np.meshgrid(r_vals, z_vals)

    # Evaluate along θ = 0 (through a strut)
    obs_x = R_grid
    obs_y = np.zeros_like(R_grid)
    obs_z = Z_grid

    Gmag = ring.grad_B(obs_x, obs_y, obs_z)

    # Mask strut interior
    inside = (np.abs(R_grid - R) < t / 2) & (np.abs(Z_grid) < L / 2)
    Gmag[inside] = np.nan

    fig, ax = plt.subplots(figsize=(10, 7))

    # Unified threshold colors (green → orange → red for low → medium → high)
    threshold_colors = {
        "40 T/m": COLORS_THRESHOLD,
        "100 T/m": COLORS_THRESHOLD_MEDIUM,
        "300 T/m": COLORS_THRESHOLD_HIGH,
    }

    im = ax.pcolormesh(
        (R_grid - R) * 1e3, Z_grid * 1e3, Gmag,
        norm=LogNorm(vmin=5, vmax=5000),
        cmap=COLORS_CONCENTRATION, shading="auto",
    )
    fig.colorbar(im, ax=ax, label="|grad_B (T/m)", shrink=0.85)

    Gmag_filled = np.nan_to_num(Gmag, nan=0)
    for lbl, val in THRESHOLDS.items():
        c = threshold_colors[lbl]
        cs = ax.contour(
            (R_grid - R) * 1e3, Z_grid * 1e3, Gmag_filled,
            levels=[val], colors=[c], linewidths=1.8,
        )
        if cs.allsegs[0]:
            ax.clabel(cs, fmt=f"{val} T/m", fontsize=9, colors=[c])

    strut = Rectangle(
        (-t / 2 * 1e3, -L / 2 * 1e3), t * 1e3, L * 1e3,
        linewidth=1.2, edgecolor="white", facecolor="none", linestyle="--",
    )
    ax.add_patch(strut)

    ax.set_xlabel("Radial distance from stent wall (mm)")
    ax.set_ylabel("Axial position z (mm)")
    ax.set_title(
        f"Gradient in r-z Plane — 3-D Akoun & Yonnet\n"
        f"({DEFAULTS['n_struts']} struts, L = {L*1e3:.1f} mm, "
        f"M = {DEFAULTS['M']/1e6:.1f} MA/m, through-strut slice)"
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 11: r-z heatmap (3-D)…")
    fig = make_figure()
    fig.savefig(OUT / "fig11_rz_heatmap.png")
    fig.savefig(OUT / "fig11_rz_heatmap.pdf")
    plt.close(fig)
    print("  [OK] fig11_rz_heatmap saved")


if __name__ == "__main__":
    main()
