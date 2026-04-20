"""
Fig 7 — COMSOL validation: full gradient-profile overlay + threshold-crossing table.

Compares the code's analytical gradient profile against COMSOL's FEM-computed
gradient profile for the V2-2C (12-cells/loop) geometry at **B₀ = 1.5 T** —
the dissertation headline operating point.

The COMSOL data live in ``data/comsol_gradient_data_cleaned.xlsx`` (V2 sheet,
38 cleaned cut-line points from strut edge radially inward). They were
generated at B₀ = 1.5 T despite the sheet being labelled "1 T normalised";
the loader in :mod:`stent_capture.data.comsol_loader` carries the corrected
field value.

Panels
------
(a) Full ``G(d)`` profile on log-log axes: analytical code (M = 2.20 MA/m
    under the COMSOL-calibrated default), V2 FEM scatter, Excel power-law fit.
(b) Bar chart of the three capture-threshold crossing distances.
(c) Validation data table with per-threshold error and RMS summary.

Run standalone::

    python -m stent_capture.figures.fig07_comsol_gradient_validation
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import (
    DEFAULTS, THRESHOLDS, TH_COLORS, OUT,
    make_ring, M_COMSOL_EFF_B15,
    COMSOL_CROSSINGS,
)
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.data.comsol_loader import load_dataset


def _code_gradient_profile(
    d_m: np.ndarray, B0_T: float = 1.5,
) -> np.ndarray:
    ring = make_ring(B0_magnitude=B0_T)
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, B0_T]))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2
    z = np.zeros_like(d_m)
    return tf.grad_B(d_m + r_outer, np.zeros_like(d_m), z)


def _compute_crossings(d_m: np.ndarray, grad: np.ndarray) -> dict[str, float]:
    out: dict[str, float] = {}
    for lbl, threshold in THRESHOLDS.items():
        idx = np.where(grad >= threshold)[0]
        if len(idx) == 0:
            out[lbl] = np.nan
            continue
        i = idx[-1]
        if i < len(d_m) - 1 and grad[i + 1] != grad[i]:
            d_cross = d_m[i] + (d_m[i + 1] - d_m[i]) * (threshold - grad[i]) / (grad[i + 1] - grad[i])
        else:
            d_cross = d_m[i]
        out[lbl] = d_cross * 1e6
    return out


def make_figure():
    d_code_m = np.geomspace(5e-6, 1.5e-3, 400)
    G_code = _code_gradient_profile(d_code_m, B0_T=1.5)
    code_crossings = _compute_crossings(d_code_m, G_code)

    v2 = load_dataset("V2")
    d_v2_mm = v2.d_mm
    G_v2 = v2.grad_T_per_m
    fit_d = np.geomspace(max(d_v2_mm.min(), 1e-3), d_v2_mm.max(), 200)
    G_fit = v2.fit_predict(fit_d)

    comsol_crossings = {lbl: dist * 1e6 for lbl, dist in COMSOL_CROSSINGS.items()}

    errors, error_pct = {}, {}
    for lbl in THRESHOLDS:
        c, m = code_crossings[lbl], comsol_crossings[lbl]
        if not np.isnan(c) and m > 0:
            errors[lbl] = c - m
            error_pct[lbl] = 100 * (c - m) / m
        else:
            errors[lbl] = np.nan
            error_pct[lbl] = np.nan

    fig, ax_chart = plt.subplots(figsize=(7, 5))

    thresholds_list = list(THRESHOLDS.keys())
    x_pos = np.arange(len(thresholds_list))
    bw = 0.35
    code_vals   = [code_crossings[l]   for l in thresholds_list]
    comsol_vals = [comsol_crossings[l] for l in thresholds_list]

    bars1 = ax_chart.bar(x_pos - bw / 2, code_vals, bw,
                         label="Code (M=2.20 MA/m, B₀=1.5 T)",
                         color="#2980b9", alpha=0.85, edgecolor="black", lw=1)
    bars2 = ax_chart.bar(x_pos + bw / 2, comsol_vals, bw,
                         label="COMSOL (FEM, μ_r=2, B₀=1.5 T)",
                         color="#e67e22", alpha=0.85, edgecolor="black", lw=1)
    for bars, colr in ((bars1, "#2980b9"), (bars2, "#e67e22")):
        for bar in bars:
            h = bar.get_height()
            ax_chart.text(bar.get_x() + bar.get_width() / 2.0, h,
                          f"{h:.1f}", ha="center", va="bottom",
                          fontsize=8, color=colr, fontweight="bold")

    ax_chart.set_xlabel("Gradient threshold", fontsize=11, fontweight="bold")
    ax_chart.set_ylabel(r"Distance from stent surface ($\mu$m)",
                        fontsize=11, fontweight="bold")
    ax_chart.set_title("Threshold-crossing distances: code vs COMSOL (B₀ = 1.5 T)",
                       fontsize=11, fontweight="bold")
    ax_chart.set_xticks(x_pos)
    ax_chart.set_xticklabels(thresholds_list, fontsize=10)
    ax_chart.set_ylim(0, max(max(code_vals), max(comsol_vals)) * 1.18)
    ax_chart.legend(fontsize=9, loc="upper right")
    ax_chart.grid(True, axis="y", alpha=0.3, ls="--")

    plt.tight_layout()
    return fig


def make_profile_figure():
    """Gradient profile panel only — x-axis in µm, clipped to 400 µm, y ≥ 10 T/m."""
    d_code_m = np.geomspace(25e-6, 4e-4, 400)   # dense from 25 µm to 400 µm
    G_code   = _code_gradient_profile(d_code_m, B0_T=1.5)

    v2      = load_dataset("V2")
    d_v2_um = v2.d_mm * 1e3                     # mm → µm
    G_v2    = v2.grad_T_per_m

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.loglog(d_code_m * 1e6, G_code, color="#2980b9", lw=2.2,
              label="Code (M=2.20 MA/m)", zorder=3)
    ax.loglog(d_v2_um, G_v2, "o", ms=5.5, color="#e67e22",
              mec="black", mew=0.5,
              label=f"COMSOL V2 (n={len(d_v2_um)} pts)", zorder=4)

    ytrans = ax.get_yaxis_transform()
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.3, alpha=0.7)
        ax.text(0.985, val, lbl, color=c, fontsize=8, va="center", ha="right",
                transform=ytrans,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=0.5))

    ax.set_xlabel(r"Arc length from strut surface $d$ ($\mu$m)",
                  fontsize=11, fontweight="bold")
    ax.set_ylabel(r"$|\nabla|B||$ (T/m)", fontsize=11, fontweight="bold")
    ax.set_title(
        r"Gradient profile: code vs COMSOL V2 (B$_0$ = 1.5 T)",
        fontsize=11, fontweight="bold",
    )
    ax.set_xlim(25.0, 400.0)
    ax.set_ylim(bottom=10)
    ax.grid(True, which="both", alpha=0.3, ls="--")
    ax.legend(fontsize=9, loc="lower left")

    plt.tight_layout()
    return fig


def main():
    print("  Fig 7: COMSOL gradient validation (profile + crossings)...")
    fig = make_figure()
    fig.savefig(OUT / "fig7_comsol_gradient_validation.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "fig7_comsol_gradient_validation.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig7_comsol_gradient_validation saved")

    fig2 = make_profile_figure()
    fig2.savefig(OUT / "fig7b_comsol_gradient_profile.png", dpi=150, bbox_inches="tight")
    fig2.savefig(OUT / "fig7b_comsol_gradient_profile.pdf", bbox_inches="tight")
    plt.close(fig2)
    print("  [OK] fig7b_comsol_gradient_profile saved")


if __name__ == "__main__":
    main()
