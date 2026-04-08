"""
Fig 9 — Axial (z) gradient profile at the stent wall for several strut lengths.

Only possible with the 3-D model — shows finite-length end-effects.

Run standalone::

    python -m stent_capture.figures.fig09_axial_profile
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring


def make_figure():
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_obs = R + t / 2 + 200e-6  # 200 µm outside stent wall, through strut (θ=0)

    L_vals = [200e-6, 500e-6, 1e-3, 2e-3, 5e-3]
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(L_vals)))

    z = np.linspace(-4e-3, 4e-3, 300)
    obs_x = np.full_like(z, r_obs)
    obs_y = np.zeros_like(z)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for L, col in zip(L_vals, colors):
        ring = make_ring(L=L)
        G = ring.grad_B(obs_x, obs_y, z)
        lbl = f"L = {L*1e3:.1f} mm"
        axes[0].semilogy(z * 1e3, G, color=col, lw=1.8, label=lbl)
        axes[1].plot(z * 1e3, G, color=col, lw=1.8, label=lbl)

    for ax in axes:
        for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
            ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.6, label=lbl)
        ax.axvline(0, color="gray", ls="--", lw=1, alpha=0.5)
        ax.set_xlabel("Axial position z (mm)")
        ax.set_ylabel("|∇|B|| (T/m)")
        ax.set_xlim(-4, 4)
        ax.legend(fontsize=8)

    axes[0].set_title("(a) Log scale")
    axes[0].set_ylim(0.1, 1e4)
    axes[1].set_title("(b) Linear scale")
    axes[1].set_ylim(0, 1200)

    fig.suptitle(
        "Axial Gradient Profile — 3-D Akoun & Yonnet\n"
        "(obs. point: r = R + t/2 + 200 µm, through strut)",
        fontsize=13, y=1.02,
    )
    plt.tight_layout()
    return fig


def main():
    print("  Fig 9: Axial gradient profile (3-D)…")
    fig = make_figure()
    fig.savefig(OUT / "fig9_axial_profile.png")
    fig.savefig(OUT / "fig9_axial_profile.pdf")
    plt.close(fig)
    print("  [OK] fig9_axial_profile saved")


if __name__ == "__main__":
    main()
