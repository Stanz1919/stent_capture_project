"""
Fig 25 — COMSOL multi-geometry gradient comparison.

Extends the single-geometry validation in ``fig04_comsol_gradient_validation``
to the full V-series (V1 = 6, V2 = 12, V3 = 18 cells per loop) plus a 2-D
linearity cross-check at two applied field strengths.

Panels
------
(a) ``G(d)`` profiles for V1 / V2 / V3 (log–log). Each geometry's COMSOL
    scatter is plotted alongside the Excel power-law fit; the analytical
    code prediction is overlaid for the V2 calibration geometry (the one
    the dissertation's M = 2.20 MA/m was calibrated against).
(b) Fitted exponent ``n(N)`` against strut count per loop. Compared with
    the empirical trend line ``n(N) = −1.094 − 2.916/N`` fitted in the
    source workbook.
(c) Fitted prefactor ``A(N)`` against strut count per loop, illustrating
    the saturation of ``A`` at ~23.5 for ``N ≥ 10``.
(d) Linearity cross-check: ratio ``G(d)[1.5 T] / G(d)[0.2433 T]`` for the
    2-D single-cell cut line, expected to equal ``1.5 / 0.2433 ≈ 6.17``
    under the COMSOL μ_r = 2 soft-ferromagnet constitutive law. Confirms
    both the label correction and the linearity of the FEM model below
    strut-material saturation.

Run standalone::

    python -m stent_capture.figures.fig25_comsol_multigeometry
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import (
    DEFAULTS, THRESHOLDS, TH_COLORS, OUT,
    make_ring, M_COMSOL_EFF_B15,
)
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.data.comsol_loader import load_dataset


# ---------------------------------------------------------------------------
# Analytical helper
# ---------------------------------------------------------------------------

def _code_gradient_profile(
    d_m: np.ndarray, *, n_struts: int, M: float, B0_T: float,
) -> np.ndarray:
    """Compute the code's analytical ``|∇|B||`` along the strut-aligned axis."""
    ring = make_ring(B0_magnitude=None, n_struts=n_struts, M=M)
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, B0_T]))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2
    z = np.zeros_like(d_m)
    return tf.grad_B(d_m + r_outer, np.zeros_like(d_m), z)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

_GEOM_COLOURS = {
    "V1": "#8e44ad",  # purple
    "V2": "#e67e22",  # orange — calibration geometry
    "V3": "#16a085",  # teal
}


def make_figure():
    v1 = load_dataset("V1")
    v2 = load_dataset("V2")
    v3 = load_dataset("V3")
    v5_new = load_dataset("V5_new")
    d2_15 = load_dataset("2D_1T")      # B_true = 1.5 T
    d2_024 = load_dataset("2D_015T")   # B_true = 0.2433 T

    fig = plt.figure(figsize=(14, 11))

    # -----------------------------------------------------------------------
    # Panel (a) — V1/V2/V3 gradient profiles + analytical overlay
    # -----------------------------------------------------------------------
    ax_a = fig.add_subplot(2, 2, 1)

    for ds in (v1, v2, v3):
        c = _GEOM_COLOURS[ds.key]
        ax_a.loglog(ds.d_mm, ds.grad_T_per_m, "o", ms=4.5, color=c,
                    mec="black", mew=0.3,
                    label=f"{ds.key}: N={ds.n_struts}  (COMSOL)", zorder=3)
        fit_d = np.geomspace(max(ds.d_mm.min(), 1e-3), ds.d_mm.max(), 200)
        ax_a.loglog(fit_d, ds.fit_predict(fit_d), "--", color=c, lw=1.2, alpha=0.75,
                    zorder=2)

    # Analytical code prediction at the V2 calibration (N=12, M=2.20 MA/m, B₀=1.5 T)
    d_code = np.geomspace(5e-6, 1.5e-3, 400)
    G_code = _code_gradient_profile(
        d_code, n_struts=12, M=M_COMSOL_EFF_B15, B0_T=1.5,
    )
    ax_a.loglog(d_code * 1e3, G_code, color="#2c3e50", lw=2.2,
                label="Code (N=12, M=2.20 MA/m, B₀=1.5 T)", zorder=4)

    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax_a.axhline(val, color=c, ls=":", lw=1.1, alpha=0.6)
        ax_a.text(1.05, val, lbl, color=c, fontsize=7.5, va="center",
                  bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=0.6))

    ax_a.set_xlabel("Distance from strut surface $d$ (mm)", fontsize=10.5, fontweight="bold")
    ax_a.set_ylabel(r"$|\nabla|B||$ (T/m)", fontsize=10.5, fontweight="bold")
    ax_a.set_title("(a) V-series gradient profiles (B₀ = 1.5 T)", fontsize=11, fontweight="bold")
    ax_a.legend(fontsize=8, loc="lower left")
    ax_a.grid(True, which="both", alpha=0.3, ls="--")
    ax_a.set_xlim(1e-2, 1.0)

    # -----------------------------------------------------------------------
    # Panel (b) — n(N) trend
    # -----------------------------------------------------------------------
    ax_b = fig.add_subplot(2, 2, 2)
    geom = [v1, v2, v3, v5_new]
    N_vals = np.array([g.n_struts for g in geom], dtype=float)
    n_vals = np.array([g.fit_n for g in geom], dtype=float)
    colours = [_GEOM_COLOURS.get(g.key, "#7f8c8d") for g in geom]

    for Nv, nv, lbl, col in zip(N_vals, n_vals, [g.key for g in geom], colours):
        ax_b.scatter(Nv, nv, s=90, color=col, edgecolor="black", lw=1.2, zorder=3,
                     label=f"{lbl} (N={int(Nv)})")

    # Empirical trend from workbook: n(N) = -1.094 - 2.916/N
    N_grid = np.linspace(4, 20, 200)
    n_trend = -1.094 - 2.916 / N_grid
    ax_b.plot(N_grid, n_trend, "--", color="#34495e", lw=1.5,
              label=r"Trend: $n(N) = -1.094 - 2.916/N$")

    ax_b.axhline(-1.094, color="#7f8c8d", ls=":", lw=1.0, alpha=0.7)
    ax_b.text(18.5, -1.08, r"Asymptote $n \to -1.094$", fontsize=8,
              color="#7f8c8d", ha="right", va="bottom")

    ax_b.set_xlabel("Struts per loop $N$", fontsize=10.5, fontweight="bold")
    ax_b.set_ylabel("Fitted exponent $n$", fontsize=10.5, fontweight="bold")
    ax_b.set_title("(b) Exponent vs strut count", fontsize=11, fontweight="bold")
    ax_b.grid(True, alpha=0.3, ls="--")
    ax_b.legend(fontsize=8, loc="lower right")
    ax_b.set_xlim(4, 20)

    # -----------------------------------------------------------------------
    # Panel (c) — A(N) trend
    # -----------------------------------------------------------------------
    ax_c = fig.add_subplot(2, 2, 3)
    A_vals = np.array([g.fit_A for g in geom], dtype=float)

    for Nv, Av, lbl, col in zip(N_vals, A_vals, [g.key for g in geom], colours):
        ax_c.scatter(Nv, Av, s=90, color=col, edgecolor="black", lw=1.2, zorder=3,
                     label=f"{lbl}: A={Av:.2f}")

    ax_c.axhline(23.5, color="#7f8c8d", ls=":", lw=1.0, alpha=0.7)
    ax_c.text(19, 23.5, "A → 23.5 (N ≥ 10)", fontsize=8, color="#7f8c8d",
              ha="right", va="bottom")

    ax_c.set_xlabel("Struts per loop $N$", fontsize=10.5, fontweight="bold")
    ax_c.set_ylabel("Fitted prefactor $A$ (T/m · mm$^{-n}$)", fontsize=10.5, fontweight="bold")
    ax_c.set_title("(c) Prefactor vs strut count", fontsize=11, fontweight="bold")
    ax_c.grid(True, alpha=0.3, ls="--")
    ax_c.legend(fontsize=8, loc="lower right")
    ax_c.set_xlim(4, 20)
    ax_c.set_ylim(0, 30)

    # -----------------------------------------------------------------------
    # Panel (d) — linearity cross-check
    # -----------------------------------------------------------------------
    ax_d = fig.add_subplot(2, 2, 4)
    B_ratio_expected = 1.5 / 0.2433

    # The two 2-D sheets share the same cut-line geometry → same d values.
    # In practice openpyxl may order rows slightly differently, so match on d.
    d_common, idx_15, idx_024 = np.intersect1d(
        np.round(d2_15.d_mm, 6), np.round(d2_024.d_mm, 6),
        return_indices=True,
    )
    g15 = d2_15.grad_T_per_m[idx_15]
    g024 = d2_024.grad_T_per_m[idx_024]
    ratio = g15 / g024

    ax_d.semilogx(d_common, ratio, "o", ms=5.5, color="#c0392b",
                  mec="black", mew=0.4, label="COMSOL ratio", zorder=3)
    ax_d.axhline(B_ratio_expected, color="#2980b9", ls="--", lw=2,
                 label=f"Linearity: 1.5/0.2433 = {B_ratio_expected:.2f}")

    mean_ratio = float(np.mean(ratio))
    ax_d.axhline(mean_ratio, color="#27ae60", ls=":", lw=1.5,
                 label=f"Observed mean: {mean_ratio:.2f}")

    ax_d.set_xlabel("Distance from strut surface $d$ (mm)", fontsize=10.5, fontweight="bold")
    ax_d.set_ylabel(r"$G(d)\,|_{B_0 = 1.5\,T}\,/\,G(d)\,|_{B_0 = 0.2433\,T}$",
                    fontsize=10.5, fontweight="bold")
    ax_d.set_title("(d) 2-D linearity check (μ_r = 2)", fontsize=11, fontweight="bold")
    ax_d.set_ylim(5.5, 7.0)
    ax_d.legend(fontsize=8, loc="lower right")
    ax_d.grid(True, which="both", alpha=0.3, ls="--")

    # -----------------------------------------------------------------------
    # Overall annotations
    # -----------------------------------------------------------------------
    note = (
        f"V2 calibration:  code M = {M_COMSOL_EFF_B15 / 1e6:.2f} MA/m "
        f"reproduces the V2 (N=12) COMSOL profile at B₀ = 1.5 T to <1% at the "
        "40–100 T/m capture thresholds. "
        f"Linearity:  mean 2-D ratio = {mean_ratio:.2f} vs expected "
        f"{B_ratio_expected:.2f} — confirms the μ_r = 2 material is in its "
        "linear regime and validates the label correction 1 T → 1.5 T, "
        "0.15 T → 0.2433 T."
    )
    fig.text(0.5, 0.01, note, ha="center", fontsize=8.5, style="italic",
             bbox=dict(boxstyle="round,pad=0.8", facecolor="lightyellow", alpha=0.85),
             wrap=True)

    fig.suptitle(
        "COMSOL multi-geometry validation — V-series gradient decay + linearity",
        fontsize=12.5, fontweight="bold", y=0.995,
    )
    plt.tight_layout(rect=[0, 0.055, 1, 0.96])
    return fig, mean_ratio


def main():
    print("  Fig 25: COMSOL multi-geometry validation...")
    fig, mean_ratio = make_figure()
    fig.savefig(OUT / "fig25_comsol_multigeometry.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "fig25_comsol_multigeometry.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] fig25_comsol_multigeometry saved  (mean 2-D ratio = {mean_ratio:.3f})")


if __name__ == "__main__":
    main()
