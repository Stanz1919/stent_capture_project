"""
Fig 1 — Single magnetised strut: |B| and |grad_B vs radial distance.

Uses the 3-D Akoun & Yonnet model with a single strut placed at the origin,
magnetised in +x.  Observation along the radial (+x) direction at z = 0.
Previously the original script compared a 2-D surface-charge model against
a dipole approximation; this figure shows the exact 3-D result for the
default strut length (L = 500 um) alongside the long-strut limit (L = 5 mm),
demonstrating the effect of finite axial length.

Run standalone::

    python -m stent_capture.figures.fig01_single_strut
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT, COLORS_CODE_CALIBRATED,
    COLORS_THRESHOLD, COLORS_THRESHOLD_MEDIUM, COLORS_THRESHOLD_HIGH
)


def make_figure():
    # Unified threshold colors (green → orange → red for low → medium → high)
    threshold_colors = {
        "40 T/m": COLORS_THRESHOLD,
        "100 T/m": COLORS_THRESHOLD_MEDIUM,
        "300 T/m": COLORS_THRESHOLD_HIGH,
    }

    p = DEFAULTS

    # Single strut at the origin, magnetised along +x.
    # Use n_struts=1, R=0: the Akoun & Yonnet code places the strut at (0,0,0)
    # with local frame = global frame for angle=0, so M_local=[M,0,0].
    ring_def = make_ring(n_struts=1, R=0.0)                    # L = 500 um
    ring_long = make_ring(n_struts=1, R=0.0, L=5e-3)          # L = 5 mm (near-infinite)

    d = np.linspace(10e-6, 2e-3, 300)
    surface = p["t"] / 2
    obs_x = d + surface
    obs_y = np.zeros_like(obs_x)
    obs_z = np.zeros_like(obs_x)

    B_def  = ring_def.B_magnitude(obs_x, obs_y, obs_z)
    G_def  = ring_def.grad_B(obs_x, obs_y, obs_z)
    B_long = ring_long.B_magnitude(obs_x, obs_y, obs_z)
    G_long = ring_long.grad_B(obs_x, obs_y, obs_z)

    d_um = d * 1e6

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.semilogy(d_um, B_def  * 1000, color=COLORS_CODE_DEFAULT, lw=2.0,
               label=f"3-D model L = {p['L']*1e3:.1f} mm")
    ax.semilogy(d_um, B_long * 1000, color=COLORS_CODE_CALIBRATED, ls="--", lw=1.8,
               label="3-D model L = 5 mm (long-strut limit)")
    ax.set_xlabel("Distance from strut surface (um)")
    ax.set_ylabel("Magnetic flux density |B| (mT)")
    ax.set_title("(a) Magnetic Flux Density vs. Distance")
    ax.legend(fontsize=10, loc='upper right', framealpha=0.95)
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3, which='both')

    ax = axes[1]
    ax.semilogy(d_um, G_def, color=COLORS_CODE_DEFAULT, lw=2.0,
               label=f"3-D model L = {p['L']*1e3:.1f} mm")
    ax.semilogy(d_um, G_long, color=COLORS_CODE_CALIBRATED, ls="--", lw=1.8,
               label="3-D model L = 5 mm (long-strut limit)")
    for lbl, val in THRESHOLDS.items():
        c = threshold_colors[lbl]
        ax.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)
    ax.set_xlabel("Distance from strut surface (um)")
    ax.set_ylabel("Field gradient |grad_B| (T/m)")
    ax.set_title("(b) Magnetic Field Gradient vs. Distance")
    ax.legend(fontsize=9, loc='upper right', framealpha=0.95)
    ax.set_xlim(0, 2000)
    ax.set_ylim(0.1, None)
    ax.grid(True, alpha=0.3, which='both')

    fig.suptitle(
        f"Figure 1: Single Magnetised Strut - 3D Akoun & Yonnet\n"
        f"(w = {p['w']*1e6:.0f} um, t = {p['t']*1e6:.0f} um, M = {p['M']/1e6:.1f} MA/m, radial magnetisation)",
        fontsize=13, fontweight='bold', y=0.98,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def main():
    print("  Fig 1: Single strut radial profile (3-D)…")
    fig = make_figure()
    fig.savefig(OUT / "fig1_single_strut.png")
    fig.savefig(OUT / "fig1_single_strut.pdf")
    plt.close(fig)
    print("  [OK] fig1_single_strut saved")


if __name__ == "__main__":
    main()
