"""
Fig 1 — Single magnetised strut: |B| and |∇|B|| vs radial distance.

Uses the 3-D Akoun & Yonnet model with a single strut placed at the origin,
magnetised in +x.  Observation along the radial (+x) direction at z = 0.
Previously the original script compared a 2-D surface-charge model against
a dipole approximation; this figure shows the exact 3-D result for the
default strut length (L = 500 µm) alongside the long-strut limit (L = 5 mm),
demonstrating the effect of finite axial length.

Run standalone::

    python -m stent_capture.figures.fig01_single_strut
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring


def make_figure():
    p = DEFAULTS

    # Single strut at the origin, magnetised along +x.
    # Use n_struts=1, R=0: the Akoun & Yonnet code places the strut at (0,0,0)
    # with local frame = global frame for angle=0, so M_local=[M,0,0].
    ring_def = make_ring(n_struts=1, R=0.0)                    # L = 500 µm
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
    ax.semilogy(d_um, B_def  * 1000, "b-",  lw=2,   label=f"3-D model  L = {p['L']*1e3:.1f} mm")
    ax.semilogy(d_um, B_long * 1000, "r--", lw=1.5, label="3-D model  L = 5 mm (long-strut)")
    ax.set_xlabel("Distance from strut surface (µm)")
    ax.set_ylabel("|B| (mT)")
    ax.set_title("(a) Magnetic flux density")
    ax.legend()
    ax.set_xlim(0, 2000)

    ax = axes[1]
    ax.semilogy(d_um, G_def,  "b-",  lw=2,   label=f"3-D model  L = {p['L']*1e3:.1f} mm")
    ax.semilogy(d_um, G_long, "r--", lw=1.5, label="3-D model  L = 5 mm (long-strut)")
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)
    ax.set_xlabel("Distance from strut surface (µm)")
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title("(b) Field gradient")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 2000)
    ax.set_ylim(0.1, None)

    fig.suptitle(
        f"Single Magnetised Strut — 3-D Akoun & Yonnet\n"
        f"(w = {p['w']*1e6:.0f} µm, t = {p['t']*1e6:.0f} µm, "
        f"M = {p['M']/1e6:.1f} MA/m, radial magnetisation)",
        fontsize=13, y=1.02,
    )
    plt.tight_layout()
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
