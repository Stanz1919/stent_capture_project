"""
Fig 4 — Effect of magnetisation on gradient and capture distance.

Run standalone::

    python -m stent_capture.figures.fig04_magnetisation_sweep
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, THRESHOLDS, TH_COLORS, OUT, make_ring


def make_figure():
    M_vals = [0.2e6, 0.5e6, 0.8e6, 1.0e6, 1.2e6]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(M_vals)))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2
    d = np.linspace(5e-6, 1.5e-3, 250)
    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    obs_z = np.zeros_like(d)
    d_um = d * 1e6

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    cap_dists: dict[str, list[float]] = {lbl: [] for lbl in THRESHOLDS}

    for M, col in zip(M_vals, colors):
        ring = make_ring(M=M)
        G = ring.grad_B(obs_x, obs_y, obs_z)
        axes[0].semilogy(d_um, G, color=col, lw=1.8,
                         label=f"M = {M/1e6:.1f} MA/m")
        for lbl, thr in THRESHOLDS.items():
            above = np.where(G >= thr)[0]
            cap_dists[lbl].append(float(d_um[above[-1]]) if len(above) > 0 else 0.0)

    ax = axes[0]
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.5)
    ax.set_xlabel("Distance from stent surface (µm)")
    ax.set_ylabel("|∇|B|| (T/m)")
    ax.set_title("(a) Gradient profiles")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 1500)
    ax.set_ylim(0.1, 1e4)

    ax = axes[1]
    M_arr = np.array(M_vals) / 1e6
    for (lbl, dists), c in zip(cap_dists.items(), TH_COLORS):
        ax.plot(M_arr, dists, "o-", color=c, lw=2, ms=8, label=lbl)
    ax.set_xlabel("Magnetisation M (MA/m)")
    ax.set_ylabel("Capture distance (µm)")
    ax.set_title("(b) Capture distance vs magnetisation")
    ax.legend(fontsize=8)

    fig.suptitle("Effect of Magnetisation on Cell Capture Range — 3-D model",
                 fontsize=13, y=1.02)
    plt.tight_layout()
    return fig


def main():
    print("  Fig 4: Magnetisation sweep…")
    fig = make_figure()
    fig.savefig(OUT / "fig4_magnetisation_sweep.png")
    fig.savefig(OUT / "fig4_magnetisation_sweep.pdf")
    plt.close(fig)
    print("  [OK] fig4_magnetisation_sweep saved")


if __name__ == "__main__":
    main()
