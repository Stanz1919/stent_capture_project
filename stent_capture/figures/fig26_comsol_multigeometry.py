"""
Fig 25 — COMSOL multi-geometry gradient comparison (panel a only).

Produces a single-panel log-log figure of G(d) for V1/V2/V3/V4.
Tabular data for the power-law fit summary (panels b/c) and the 2-D
linearity check (panel d) are written to the results directory as CSV files:

  results/fig25_powerlawfit_summary.csv  — A, n, R² per geometry
  results/fig25_linearity_check.csv      — per-node ratio G[1.5T]/G[0.2433T]

V1 data source: COMSOL 6.1 export ``V1-Mgrad-Plot-New2.csv``
(v1-2c-smoothed.mph, 46 nodes).  Fillet-surface cluster (d ≤ 0.025 mm) and
noise (grad < 1 T/m) are excluded, leaving ~15 usable near-field points.

V4 data: near-field only, d < 0.170 mm.

Run standalone::

    python -m stent_capture.figures.fig26_comsol_multigeometry
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import (
    DEFAULTS, THRESHOLDS, TH_COLORS, OUT,
    make_ring, M_COMSOL_EFF_B15,
)
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.data.comsol_loader import load_dataset

# ---------------------------------------------------------------------------
# V1 CSV loader
# ---------------------------------------------------------------------------

_V1_CSV_PATH = Path(__file__).parents[2] / "data" / "V1-Mgrad-Plot-New2.csv"


def _load_v1_csv() -> tuple[np.ndarray, np.ndarray]:
    """Load V1 cut-line CSV, compute d = r_max − r, apply quality filters.

    Excludes d ≤ 0.025 mm (fillet-surface non-monotone cluster) and
    grad < 1 T/m (noise floor).  Returns (d_mm, grad_T_per_m) sorted by d.
    """
    rows: list[list[float]] = []
    with open(_V1_CSV_PATH) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split(",")
            if len(parts) < 4:
                continue
            try:
                rows.append([float(p) for p in parts[:4]])
            except ValueError:
                continue

    arr = np.array(rows)
    x, y, z, grad = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    r = np.sqrt(x**2 + y**2 + z**2)
    d_mm = r.max() - r

    mask = (d_mm > 0.025) & (grad >= 1.0)
    d_mm, grad = d_mm[mask], grad[mask]
    order = np.argsort(d_mm)
    return d_mm[order], grad[order]


# ---------------------------------------------------------------------------
# Analytical helper
# ---------------------------------------------------------------------------

def _code_gradient_profile(
    d_m: np.ndarray, *, n_struts: int, M: float, B0_T: float,
) -> np.ndarray:
    """Compute the code's analytical |grad_B| along the strut-aligned axis."""
    ring = make_ring(B0_magnitude=None, n_struts=n_struts, M=M)
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, B0_T]))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2
    return tf.grad_B(d_m + r_outer, np.zeros_like(d_m), np.zeros_like(d_m))


# ---------------------------------------------------------------------------
# Tabular exports (panels b, c, d)
# ---------------------------------------------------------------------------

def _export_powerlawfit(geom_datasets, labels: list[str]) -> None:
    """Write power-law fit summary (panels b/c) to CSV."""
    path = OUT / "fig25_powerlawfit_summary.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["geometry", "N_struts", "fit_A", "fit_n", "fit_R2",
                    "fit_range_mm"])
        for ds, lbl in zip(geom_datasets, labels):
            rng = (
                f"{ds.fit_range_mm[0]}-{ds.fit_range_mm[1]}"
                if ds.fit_range_mm else "full"
            )
            w.writerow([lbl, ds.n_struts,
                        f"{ds.fit_A:.4f}" if ds.fit_A is not None else "",
                        f"{ds.fit_n:.4f}" if ds.fit_n is not None else "",
                        f"{ds.fit_R2:.4f}" if ds.fit_R2 is not None else "",
                        rng])
    print(f"  [CSV] {path.name}")


def _export_linearity(d_common, g15, g024, ratio, mean_ratio) -> None:
    """Write 2-D linearity check data (panel d) to CSV."""
    path = OUT / "fig25_linearity_check.csv"
    expected = 1.5 / 0.2433
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([f"# Expected ratio (1.5/0.2433): {expected:.4f}",
                    f"mean observed: {mean_ratio:.4f}"])
        w.writerow(["d_mm", "grad_1p5T_T_per_m", "grad_0p2433T_T_per_m",
                    "ratio"])
        for d, ga, gb, r in zip(d_common, g15, g024, ratio):
            w.writerow([f"{d:.6f}", f"{ga:.4f}", f"{gb:.4f}", f"{r:.4f}"])
    print(f"  [CSV] {path.name}")


# ---------------------------------------------------------------------------
# Figure (single panel)
# ---------------------------------------------------------------------------

_GEOM_COLOURS = {
    "V1": "#8e44ad",   # purple
    "V2": "#e67e22",   # orange — calibration geometry
    "V3": "#16a085",   # teal
    "V4": "#2980b9",   # blue
}

_KEY_TO_COLOUR_KEY = {"V4_new": "V4"}


def _colour(key: str) -> str:
    return _GEOM_COLOURS.get(_KEY_TO_COLOUR_KEY.get(key, key), "#7f8c8d")


def make_figure(y_min: float | None = None, x_min: float | None = None):
    v1_excel = load_dataset("V1")
    v2       = load_dataset("V2")
    v3       = load_dataset("V3")
    v5_new   = load_dataset("V4_new")
    d2_15    = load_dataset("2D_1T")
    d2_024   = load_dataset("2D_015T")

    v1_d_mm, v1_grad = _load_v1_csv()

    # V4 near-field filter: d < 0.170 mm (confirmed still applied)
    v5_mask = v5_new.d_mm < 0.170
    v5_d    = v5_new.d_mm[v5_mask]
    v5_g    = v5_new.grad_T_per_m[v5_mask]

    # Linearity data (needed for CSV export and note)
    B_ratio_expected = 1.5 / 0.2433
    d_common, idx_15, idx_024 = np.intersect1d(
        np.round(d2_15.d_mm, 6), np.round(d2_024.d_mm, 6),
        return_indices=True,
    )
    g15   = d2_15.grad_T_per_m[idx_15]
    g024  = d2_024.grad_T_per_m[idx_024]
    ratio = g15 / g024
    mean_ratio = float(np.mean(ratio))

    # -----------------------------------------------------------------------
    # Single-panel figure
    # -----------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6.5))

    # V2 and V3: filled circles + Excel power-law fit dashes
    for ds in (v2, v3):
        c = _colour(ds.key)
        ax.loglog(ds.d_mm, ds.grad_T_per_m, "o", ms=5, color=c,
                  mec="black", mew=0.3,
                  label=f"{ds.key}: N={ds.n_struts}  (COMSOL)", zorder=3)
        fit_d = np.geomspace(max(ds.d_mm.min(), 1e-3), ds.d_mm.max(), 200)
        ax.loglog(fit_d, ds.fit_predict(fit_d), "--", color=c, lw=1.3,
                  alpha=0.75, zorder=2)

    # V4: near-field data + power-law fit (d < 0.170 mm)
    c5 = _colour("V4")
    ax.loglog(v5_d, v5_g, "o", ms=5, color=c5,
              mec="black", mew=0.3,
              label=r"V4: N=10  (COMSOL, $d<0.17$ mm)", zorder=3)
    if v5_new.fit_A is not None and len(v5_d) > 1:
        fit_d5 = np.geomspace(max(v5_d.min(), 1e-3), 0.170, 200)
        ax.loglog(fit_d5, v5_new.fit_predict(fit_d5), "--", color=c5,
                  lw=1.3, alpha=0.75, zorder=2)

    # V1: open upward-triangle marker; NO power-law fit
    c1 = _colour("V1")
    ax.loglog(v1_d_mm, v1_grad, "^", ms=6.5, color=c1,
              mfc="none", mec=c1, mew=1.4,
              label=r"V1: N=6  (COMSOL CSV, $d>25\,\mu$m)", zorder=3)

    # Analytical Akoun–Yonnet curves for all four geometries
    # All use M_eff = 2.20 MA/m calibrated to V2 (N=12)
    d_code = np.geomspace(5e-6, 1.5e-3, 400)
    for n_struts, key in [(6, "V1"), (12, "V2"), (18, "V3"), (10, "V4")]:
        G_code = _code_gradient_profile(
            d_code, n_struts=n_struts, M=M_COMSOL_EFF_B15, B0_T=1.5,
        )
        lw  = 2.2 if key == "V2" else 1.5
        lbl = f"Analytical N={n_struts}"
        if key == "V2":
            lbl += r" ($M_\mathrm{eff}=2.20$ MA/m, cal.)"
        ax.loglog(d_code * 1e3, G_code, color=_colour(key), lw=lw,
                  alpha=0.9, zorder=4, label=lbl)

    # Semi-reliable boundary at d = 0.200 mm
    ax.axvline(0.200, color="#7f8c8d", ls="--", lw=1.5, alpha=0.9, zorder=1)
    trans = ax.get_xaxis_transform()
    ax.text(0.210, 0.97, "semi-reliable\nbeyond here",
            transform=trans, fontsize=8, va="top", ha="left",
            color="#636e72",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.5))

    # V1 annotation
    ax.text(0.02, 0.04,
            r"V1 ($\triangle$): single $M_\mathrm{eff}$ calibrated to V2" "\n"
            r"— analytical over-predicts V1; expected",
            transform=ax.transAxes, fontsize=8, ha="left", va="bottom",
            style="italic", color=c1,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.5))

    # Threshold lines (labels placed inside the axes using mixed transform)
    ytrans = ax.get_yaxis_transform()   # x=axes fraction, y=data
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.2, alpha=0.65)
        ax.text(0.985, val, lbl, color=c, fontsize=8, va="center", ha="right",
                transform=ytrans,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.65,
                          pad=0.5))

    ax.set_xlabel("Distance from strut surface $d$ (mm)",
                  fontsize=11, fontweight="bold")
    ax.set_ylabel(r"$|\nabla|B||$ (T/m)", fontsize=11, fontweight="bold")
    ax.set_title(
        "V-series gradient profiles — COMSOL vs analytical (B₀ = 1.5 T)",
        fontsize=11.5, fontweight="bold",
    )
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, which="both", alpha=0.3, ls="--")
    ax.set_xlim(x_min if x_min is not None else 1e-2, 1.0)
    if y_min is not None:
        ax.set_ylim(bottom=y_min)

    note = (
        f"Analytical curves: $M_\\mathrm{{eff}}$ = {M_COMSOL_EFF_B15/1e6:.2f} MA/m "
        "calibrated to V2 (N=12); same value applied to all geometries.  "
        "V1 analytical over-predicts — V1 sits at lower gradient amplitude "
        "than V2/V3/V4 at all depths (expected for fewer struts, N=6).  "
        "V4 near-field filter d < 0.170 mm applied."
    )
    fig.text(0.5, 0.01, note, ha="center", fontsize=8, style="italic",
             bbox=dict(boxstyle="round,pad=0.7", facecolor="lightyellow",
                       alpha=0.85),
             wrap=True)

    plt.tight_layout(rect=[0, 0.07, 1, 1])
    return fig, v1_excel, v2, v3, v5_new, d_common, g15, g024, ratio, mean_ratio


def main():
    print("  Fig 25: COMSOL multi-geometry gradient profiles...")
    (fig, v1_excel, v2, v3, v5_new,
     d_common, g15, g024, ratio, mean_ratio) = make_figure()

    fig.savefig(OUT / "fig25_comsol_multigeometry.png", dpi=150,
                bbox_inches="tight")
    fig.savefig(OUT / "fig25_comsol_multigeometry.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig25_comsol_multigeometry saved")

    # Variant with y-axis clipped to 10 T/m minimum
    (fig2, *_) = make_figure(y_min=10, x_min=0.025)
    fig2.savefig(OUT / "fig25_comsol_multigeometry_ymin10.png", dpi=150,
                 bbox_inches="tight")
    fig2.savefig(OUT / "fig25_comsol_multigeometry_ymin10.pdf",
                 bbox_inches="tight")
    plt.close(fig2)
    print("  [OK] fig25_comsol_multigeometry_ymin10 saved")

    # Tabular exports
    geom_datasets = [v1_excel, v2, v3, v5_new]
    labels        = ["V1", "V2", "V3", "V4"]
    _export_powerlawfit(geom_datasets, labels)
    _export_linearity(d_common, g15, g024, ratio, mean_ratio)


if __name__ == "__main__":
    main()
