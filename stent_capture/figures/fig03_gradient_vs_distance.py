"""
Fig 3 — Gradient vs distance: through-strut and between-struts profiles.

Run standalone::

    python -m stent_capture.figures.fig03_gradient_vs_distance
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring


def make_figure():
    ring = make_ring()
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    d = np.linspace(5e-6, 2e-3, 300)
    z = np.zeros_like(d)

    # Through a strut (θ = 0)
    G_through = ring.grad_B(d + r_outer, np.zeros_like(d), z)

    # Between struts (θ = π / n_struts)
    angle_bet = np.pi / DEFAULTS["n_struts"]
    r_bet = d + r_outer
    G_between = ring.grad_B(
        r_bet * np.cos(angle_bet),
        r_bet * np.sin(angle_bet),
        z,
    )

    d_um = d * 1e6

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for ax, yscale in zip(axes, ["linear", "log"]):
        plot_fn = ax.plot if yscale == "linear" else ax.semilogy
        plot_fn(d_um, G_through,  "b-",  lw=2,   label="Through strut")
        plot_fn(d_um, G_between,  "b--", lw=1.5, alpha=0.7, label="Between struts")
        for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
            ax.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)
        ax.set_xlabel("Distance from stent surface (µm)")
        ax.set_ylabel("|∇|B|| (T/m)")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 2000)

    axes[0].set_title("(a) Linear scale")
    axes[0].set_ylim(0, min(3000, float(np.nanmax(G_through)) * 1.1))
    axes[1].set_title("(b) Log scale")
    axes[1].set_ylim(0.1, 1e4)

    fig.suptitle(
        f"Gradient vs Distance from Stent Surface — 3-D model\n"
        f"({DEFAULTS['n_struts']} struts, M = {DEFAULTS['M']/1e6:.1f} MA/m)",
        fontsize=13, y=1.02,
    )
    plt.tight_layout()
    return fig


def main():
    print("  Fig 3: Gradient vs distance…")
    fig = make_figure()
    fig.savefig(OUT / "fig3_gradient_vs_distance.png")
    fig.savefig(OUT / "fig3_gradient_vs_distance.pdf")
    plt.close(fig)
    print("  [OK] fig3_gradient_vs_distance saved")


if __name__ == "__main__":
    main()
