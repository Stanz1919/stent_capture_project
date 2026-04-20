"""
Fig 9 — SPION loading parameter sweep.

2x1 panel figure:

(a) Capture distance (um, from stent inner surface into lumen) vs SPION
    loading (1-300 pg), for mean blood velocities 0.05, 0.2, 0.5 m/s.
    B0 = 1.5 T axial (MRI).  Capture distance is defined as the outermost
    distance from the stent inner surface at which |F_mag| ≥ |F_drag|,
    swept along the through-strut radial line into the lumen.

(b) Force ratio |F_mag|/|F_drag| at a fixed point 5 um from the stent
    inner surface, vs SPION loading.  Log-log scale.  Horizontal line
    at ratio = 1 (capture threshold).

Overlays on (a):
  - Horizontal dashed line: inter-strut half-arc (π R / n_struts  ~  589 um)
  - Horizontal dashed line: lumen inner radius (R − t/2  ~  1460 um)
  - Vertical dotted lines at 10, 50, 200 pg (literature benchmarks)
  - Shaded band 30-100 pg (typical experimental range)

The expensive field computation (|B| and grad_B along the radial line) is
performed ONCE; subsequent sweeps over loading and velocity values just
scale the precomputed force-parameter array — total runtime  ~  3-5 s.

Run standalone::

    python -m stent_capture.figures.fig09_spion_loading_sweep
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT, COLORS_CODE_CALIBRATED, COLORS_THRESHOLD_HIGH,
    COLORS_THRESHOLD, COLORS_MARKER_REFERENCE
)
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, MU_0
from stent_capture.physics.hydrodynamics import BloodFlow, stokes_drag
from stent_capture.core.gradient import compute_gradient_vector

_B0_Z    = 1.5
_R_VES   = 1.54e-3
_V_CASES = [0.05, 0.2, 0.4, 0.6]
_COLORS  = [
    COLORS_CODE_DEFAULT,
    COLORS_CODE_CALIBRATED,
    COLORS_THRESHOLD_HIGH,
    "#8e44ad",
]
_CELL    = SPIONLabelledCell()

_LOADINGS_PG = np.logspace(0, np.log10(300), 50)
_LOADINGS_KG = _LOADINGS_PG * 1e-15

_BENCH_PG  = [200]
_BENCH_LABELS = {
    200: "Polyak experimental (200 pg)",
}


def _precompute():
    """
    Returns
    -------
    d_vals : (N,) array — distance from stent inner surface (m)
    force_per_kg : (N,) array — |F_mag| per kg of SPION mass, B0=0.5T axial
        i.e.  |F_mag|(m) = m_kg x force_per_kg   [N / kg]
    drag_profiles : dict  v_mean → (N,) |F_drag| in N
    """
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_inner = R - t / 2

    N_pts = 400
    d_vals = np.linspace(2e-6, r_inner * 0.999, N_pts)
    x_vals = r_inner - d_vals
    pts = np.column_stack([x_vals, np.zeros(N_pts), np.zeros(N_pts)])

    print("    Precomputing B-field and gradient along radial line...")
    ring = make_ring(B0_magnitude=_B0_Z)
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    B_vecs  = tf.field_at(pts)
    B_mag   = np.linalg.norm(B_vecs, axis=1)
    grad_v  = compute_gradient_vector(tf.field_at, pts)
    grad_mag = np.linalg.norm(grad_v, axis=1)

    chi   = _CELL.spion_susceptibility
    rho   = _CELL.spion_density
    force_per_kg = (chi / (rho * MU_0)) * B_mag * grad_mag

    print("    Precomputing Stokes drag for each velocity...")
    drag_profiles = {}
    for v in _V_CASES:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        F = stokes_drag(_CELL, flow, pts)
        drag_profiles[v] = np.linalg.norm(F, axis=1)

    return d_vals, force_per_kg, drag_profiles


def _capture_dist_curve(
    loadings_kg: np.ndarray,
    force_per_kg: np.ndarray,
    drag: np.ndarray,
    d_vals: np.ndarray,
) -> np.ndarray:
    """Return capture distance (um) for each loading value."""
    result = np.zeros(len(loadings_kg))
    for j, m in enumerate(loadings_kg):
        F_mag = m * force_per_kg
        captured = F_mag > drag
        if np.any(captured):
            result[j] = float(d_vals[captured][-1]) * 1e6
    return result


def make_figure():
    d_vals, force_per_kg, drag_profiles = _precompute()

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    n = DEFAULTS["n_struts"]

    inter_strut_half = np.pi * R / n * 1e6
    lumen_radius_um  = (R - t / 2) * 1e6

    idx_5um = int(np.argmin(np.abs(d_vals - 5e-6)))

    cap_dist: dict[float, np.ndarray] = {}
    ratio_5um: dict[float, np.ndarray] = {}

    for v, col in zip(_V_CASES, _COLORS):
        drag = drag_profiles[v]
        cap_dist[v] = _capture_dist_curve(_LOADINGS_KG, force_per_kg, drag, d_vals)

        F_drag_5um = drag[idx_5um]
        if F_drag_5um > 0:
            ratio_5um[v] = _LOADINGS_KG * force_per_kg[idx_5um] / F_drag_5um
        else:
            ratio_5um[v] = np.full_like(_LOADINGS_KG, np.nan)


    bench_kg = [m * 1e-15 for m in _BENCH_PG]
    print("\n    Capture distance table (um from stent inner surface):")
    header = f"    {'Loading':>12}  " + "  ".join(f"v={v:.2f} m/s{''!s:>4}" for v in _V_CASES)
    print(header)
    table_data: dict[int, dict[float, float]] = {}
    for m_pg, m_kg in zip(_BENCH_PG, bench_kg):
        row: dict[float, float] = {}
        parts = [f"    {m_pg:>10} pg"]
        for v in _V_CASES:
            drag = drag_profiles[v]
            F_mag = m_kg * force_per_kg
            captured = F_mag > drag
            d_cap = float(d_vals[captured][-1]) * 1e6 if np.any(captured) else 0.0
            row[v] = d_cap
            parts.append(f"  {d_cap:>10.1f} um")
        table_data[m_pg] = row
        print("".join(parts))


    print("\n    Loading for capture distance first > 100 um:")
    for v in _V_CASES:
        cap_arr = cap_dist[v]
        indices = np.where(cap_arr > 100)[0]
        if len(indices):
            m_thresh = _LOADINGS_PG[indices[0]]
            print(f"      v = {v:.2f} m/s : {m_thresh:.1f} pg")
        else:
            print(f"      v = {v:.2f} m/s : > {_LOADINGS_PG[-1]:.0f} pg (not reached in sweep)")


    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6))


    for v, col in zip(_V_CASES, _COLORS):
        ax_a.semilogx(_LOADINGS_PG, cap_dist[v],
                      color=col, lw=2.5, label=f"v_mean = {v:.2f} m/s")



    ax_a.axvspan(5, 50, alpha=0.12, color=COLORS_THRESHOLD,
                 label="Literature range (5–50 pg)")


    bench_styles = [":"]
    bench_colors = ["#16a085"]
    for m_pg, ls, bc in zip(_BENCH_PG, bench_styles, bench_colors):
        ax_a.axvline(m_pg, color=bc, ls=ls, lw=1.5,
                     label=_BENCH_LABELS[m_pg])

    ax_a.set_xlabel("SPION loading (pg iron oxide per cell)")
    ax_a.set_ylabel(r"Capture distance, $d$ ($\mu$m)")
    ax_a.set_title("(a) Capture distance vs SPION loading\n(B0 = 1.5 T axial, inward radial sweep)")
    ax_a.set_xlim(1, 300)
    ax_a.set_ylim(bottom=0)
    ax_a.legend(fontsize=7, loc="upper left")
    ax_a.grid(True, which="both", alpha=0.3)


    for v, col in zip(_V_CASES, _COLORS):
        ax_b.loglog(_LOADINGS_PG, ratio_5um[v],
                    color=col, lw=2.5, label=f"v_mean = {v:.2f} m/s")

    ax_b.axhline(1.0, color=COLORS_MARKER_REFERENCE, ls="-", lw=1.5,
                 label="Ratio = 1.0 (capture threshold)")
    ax_b.axvspan(5, 50, alpha=0.12, color=COLORS_THRESHOLD,
                 label="Literature range (5–50 pg)")

    for m_pg, ls, bc in zip(_BENCH_PG, bench_styles, bench_colors):
        ax_b.axvline(m_pg, color=bc, ls=ls, lw=1.5, label=_BENCH_LABELS[m_pg])

    ax_b.set_xlabel("SPION loading (pg iron oxide per cell)")
    ax_b.set_ylabel(r"$|F_\mathrm{mag}| / |F_\mathrm{drag}|$")
    ax_b.set_title("(b) Force ratio at 5 µm from stent inner surface\n(B0 = 1.5 T axial; ratio > 1 = capture)")
    ax_b.set_xlim(1, 300)
    ax_b.set_ylim(0.02, 1e2)
    ax_b.legend(fontsize=7, loc="upper left")
    ax_b.grid(True, which="both", alpha=0.3)

    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 9: SPION loading sweep...")
    fig = make_figure()
    fig.savefig(OUT / "fig9_spion_loading_sweep.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig9_spion_loading_sweep.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig9_spion_loading_sweep saved")


if __name__ == "__main__":
    main()
