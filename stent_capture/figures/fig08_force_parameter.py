"""
Fig 8 — B·grad_B force parameter (proportional to magnetic capture force).

Run standalone::

    python -m stent_capture.figures.fig08_force_parameter
"""

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.figures.style import COLORS_CODE_DEFAULT, COLORS_THRESHOLD_HIGH


def make_figure():
    ring = make_ring()
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    d = np.linspace(5e-6, 1.5e-3, 250)
    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    obs_z = np.zeros_like(d)

    Bmag = ring.B_magnitude(obs_x, obs_y, obs_z)
    Gmag = ring.grad_B(obs_x, obs_y, obs_z)
    F_param = Bmag * Gmag  # T²/m

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    ax.semilogy(d * 1e6, F_param, color=COLORS_CODE_DEFAULT, lw=2)
    ax.set_xlabel("Distance from stent surface (um)")
    ax.set_ylabel("B·grad_B (T²/m)")
    ax.set_title("(a) Force parameter vs distance")

    ax = axes[1]
    l1, = ax.semilogy(d * 1e6, Bmag * 1000, color=COLORS_CODE_DEFAULT, lw=2, label="|B| (mT)")
    ax2 = ax.twinx()
    l2, = ax2.semilogy(d * 1e6, Gmag, color=COLORS_THRESHOLD_HIGH, lw=2, label="|grad_B (T/m)")
    ax.set_xlabel("Distance from stent surface (um)")
    ax.set_ylabel("|B| (mT)", color=COLORS_CODE_DEFAULT)
    ax2.set_ylabel("|grad_B (T/m)", color=COLORS_THRESHOLD_HIGH)
    ax.set_title("(b) B and gradient decomposition")
    ax.legend([l1, l2], [l1.get_label(), l2.get_label()], fontsize=9)

    fig.suptitle("Magnetic Force Parameter for Cell Capture — 3-D model",
                 fontsize=13, y=0.98)
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 8: Force parameter…")
    fig = make_figure()
    fig.savefig(OUT / "fig8_force_parameter.png")
    fig.savefig(OUT / "fig8_force_parameter.pdf")
    plt.close(fig)
    print("  [OK] fig8_force_parameter saved")


if __name__ == "__main__":
    main()
