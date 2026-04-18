"""
Fig 17 — SPION loading parameter sweep.

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

    python -m stent_capture.figures.fig17_spion_loading_sweep
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

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

_B0_Z    = 1.5               # T, axial — MRI-strength, matches COMSOL
_R_VES   = 1.54e-3           # m, vessel radius
_V_CASES = [0.05, 0.2, 0.4, 0.6]  # m/s  (0.6 = healthy MCA, Aaslid et al. 1982)
# Unified color palette: blue (low), orange (medium), red (high), purple (healthy MCA)
_COLORS  = [
    COLORS_CODE_DEFAULT,      # 0.05 m/s - low velocity (blue)
    COLORS_CODE_CALIBRATED,   # 0.2 m/s - medium velocity (dark orange)
    COLORS_THRESHOLD_HIGH,    # 0.4 m/s - high velocity (red)
    "#8e44ad",                # 0.6 m/s - healthy MCA (purple)
]
_CELL    = SPIONLabelledCell()   # default (10 pg); used only for drag

# Loading sweep: 1-1000 pg, 50 points log-spaced
_LOADINGS_PG = np.logspace(0, 3, 50)
_LOADINGS_KG = _LOADINGS_PG * 1e-15

# Reference loadings for the 3x3 table and annotations
_BENCH_PG  = [10, 50, 200]
_BENCH_LABELS = {
    10:  "Polyak 2008 default",
    50:  "Chorny 2007 typical",
    200: "Riegler 2011 upper",
}

# ---------------------------------------------------------------------------
# Precompute
# ---------------------------------------------------------------------------

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
    r_inner = R - t / 2   # stent inner surface,  ~  1.460 mm

    # Radial line: d = distance from stent inner surface into lumen
    N_pts = 400
    d_vals = np.linspace(2e-6, r_inner * 0.999, N_pts)
    x_vals = r_inner - d_vals
    pts = np.column_stack([x_vals, np.zeros(N_pts), np.zeros(N_pts)])

    print("    Precomputing B-field and gradient along radial line...")
    ring = make_ring(B0_magnitude=_B0_Z)  # Adaptive M for COMSOL calibration at 1.5 T
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    B_vecs  = tf.field_at(pts)                          # (N, 3)
    B_mag   = np.linalg.norm(B_vecs, axis=1)            # (N,)
    grad_v  = compute_gradient_vector(tf.field_at, pts) # (N, 3)
    grad_mag = np.linalg.norm(grad_v, axis=1)           # (N,)

    # |F_mag| = (m / rho * chi / mu0) * |B| * |grad_B
    #         = m * (chi / (rho * mu0)) * |B| * |grad_B
    chi   = _CELL.spion_susceptibility   # 2.0
    rho   = _CELL.spion_density          # 5170 kg/m³
    force_per_kg = (chi / (rho * MU_0)) * B_mag * grad_mag   # N/kg

    print("    Precomputing Stokes drag for each velocity...")
    drag_profiles = {}
    for v in _V_CASES:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        F = stokes_drag(_CELL, flow, pts)
        drag_profiles[v] = np.linalg.norm(F, axis=1)   # N

    return d_vals, force_per_kg, drag_profiles


# ---------------------------------------------------------------------------
# Compute capture distance curves
# ---------------------------------------------------------------------------

def _capture_dist_curve(
    loadings_kg: np.ndarray,
    force_per_kg: np.ndarray,
    drag: np.ndarray,
    d_vals: np.ndarray,
) -> np.ndarray:
    """Return capture distance (um) for each loading value."""
    result = np.zeros(len(loadings_kg))
    for j, m in enumerate(loadings_kg):
        F_mag = m * force_per_kg     # (N,)
        captured = F_mag > drag      # (N,) bool
        if np.any(captured):
            result[j] = float(d_vals[captured][-1]) * 1e6   # um
    return result


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    d_vals, force_per_kg, drag_profiles = _precompute()

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    n = DEFAULTS["n_struts"]

    inter_strut_half = np.pi * R / n * 1e6     # um   ~  589
    lumen_radius_um  = (R - t / 2) * 1e6       # um   ~  1460

    # Index for 5 um from stent inner surface
    idx_5um = int(np.argmin(np.abs(d_vals - 5e-6)))

    # -----------------------------------------------------------------------
    # Compute capture-distance curves and force-ratio at 5 um
    # -----------------------------------------------------------------------
    cap_dist: dict[float, np.ndarray] = {}
    ratio_5um: dict[float, np.ndarray] = {}

    for v, col in zip(_V_CASES, _COLORS):
        drag = drag_profiles[v]
        cap_dist[v] = _capture_dist_curve(_LOADINGS_KG, force_per_kg, drag, d_vals)
        # Force ratio at 5 um: F_mag / F_drag = (m * force_per_kg[idx]) / drag[idx]
        F_drag_5um = drag[idx_5um]
        if F_drag_5um > 0:
            ratio_5um[v] = _LOADINGS_KG * force_per_kg[idx_5um] / F_drag_5um
        else:
            ratio_5um[v] = np.full_like(_LOADINGS_KG, np.nan)

    # -----------------------------------------------------------------------
    # 3x3 table: capture distances at benchmark loadings x velocities
    # -----------------------------------------------------------------------
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

    # -----------------------------------------------------------------------
    # Threshold loading: first loading where capture distance > 100 um
    # -----------------------------------------------------------------------
    print("\n    Loading for capture distance first > 100 um:")
    for v in _V_CASES:
        cap_arr = cap_dist[v]
        indices = np.where(cap_arr > 100)[0]
        if len(indices):
            m_thresh = _LOADINGS_PG[indices[0]]
            print(f"      v = {v:.2f} m/s : {m_thresh:.1f} pg")
        else:
            print(f"      v = {v:.2f} m/s : > 1000 pg (not reached in sweep)")

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6))

    # --- Panel (a): capture distance vs loading ---
    for v, col in zip(_V_CASES, _COLORS):
        ax_a.semilogx(_LOADINGS_PG, cap_dist[v],
                      color=col, lw=2.5, label=f"v_mean = {v:.2f} m/s")

    # Overlays
    # Shaded band: typical experimental range 30-100 pg
    ax_a.axvspan(30, 100, alpha=0.10, color=COLORS_THRESHOLD,
                 label="Typical experimental\nrange (30-100 pg)")

    # Inter-strut half-distance annotated as text (589 um exceeds capture range at all loadings tested)
    ax_a.annotate(
        f"Inter-strut half-distance: {inter_strut_half:.0f} µm",
        xy=(1, 1), xycoords="axes fraction",
        xytext=(-5, -5), textcoords="offset points",
        ha="right", va="top", fontsize=8,
        color=COLORS_MARKER_REFERENCE,
    )

    # Vertical reference lines (literature benchmarks)
    bench_styles = [":", ":", ":"]
    bench_colors = ["#2c3e50", "#8e44ad", "#16a085"]
    for m_pg, ls, bc in zip(_BENCH_PG, bench_styles, bench_colors):
        ax_a.axvline(m_pg, color=bc, ls=ls, lw=1.5,
                     label=_BENCH_LABELS[m_pg])

    ax_a.set_xlabel("SPION loading (pg iron oxide per cell)")
    ax_a.set_ylabel("Capture distance from stent inner surface (um)")
    ax_a.set_title(
        "(a) Capture distance vs SPION loading\n"
        "(B0 = 1.5 T axial (MRI), inward radial sweep)"
    )
    ax_a.set_xlim(1, 1000)
    ax_a.set_ylim(bottom=0)
    ax_a.legend(fontsize=7, loc="upper left")
    ax_a.grid(True, which="both", alpha=0.3)

    # --- Panel (b): force ratio at 5 um from inner surface ---
    for v, col in zip(_V_CASES, _COLORS):
        ax_b.loglog(_LOADINGS_PG, ratio_5um[v],
                    color=col, lw=2.5, label=f"v_mean = {v:.2f} m/s")

    ax_b.axhline(1.0, color=COLORS_MARKER_REFERENCE, ls="-", lw=1.5,
                 label="Ratio = 1.0 (capture threshold)")
    ax_b.axvspan(30, 100, alpha=0.10, color=COLORS_THRESHOLD, label="Typical experimental range")

    for m_pg, ls, bc in zip(_BENCH_PG, bench_styles, bench_colors):
        ax_b.axvline(m_pg, color=bc, ls=ls, lw=1.5, label=_BENCH_LABELS[m_pg])

    ax_b.set_xlabel("SPION loading (pg iron oxide per cell)")
    ax_b.set_ylabel("|F_mag| / |F_drag|  (at 5 um from stent inner surface)")
    ax_b.set_title(
        "(b) Force ratio at 5 um from stent inner surface\n"
        "(B0 = 1.5 T axial (MRI); ratio > 1 = capture)"
    )
    ax_b.set_xlim(1, 1000)
    ax_b.set_ylim(1e-3, 1e2)
    ax_b.legend(fontsize=7, loc="upper left")
    ax_b.grid(True, which="both", alpha=0.3)

    fig.suptitle(
        "SPION loading parameter sweep  —  B0 = 1.5 T axial (MRI), 12-strut stent (V2-2C), R = 1.5 mm\n"
        "Static Furlani & Ng criterion: |F_mag| > |F_drag|  (Stage 3 will extend to trajectories)\n"
        "Cell: 10 um radius, χ = 2.0, magnetite (5170 kg/m³);  vessel R = 1.54 mm, η = 4 mPa·s\n"
        "v = 0.6 m/s: healthy MCA mean (Aaslid et al. 1982);  v = 0.2 m/s: distal/diseased vessel",
        fontsize=10, y=0.98,
    )
    plt.tight_layout(pad=1.0)
    return fig


def main():
    print("  Fig 17: SPION loading sweep...")
    fig = make_figure()
    fig.savefig(OUT / "fig17_spion_loading_sweep.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig17_spion_loading_sweep.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig17_spion_loading_sweep saved")


if __name__ == "__main__":
    main()
