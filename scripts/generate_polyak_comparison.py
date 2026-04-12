"""
Polyak (2008) Comparison Run — B₀ = 0.1 T
==========================================

Reproduces the physical conditions closest to Polyak et al. (2008) PNAS:
  - External field:  B₀ = 0.1 T  (1,000 G, as used experimentally by Polyak)
  - Cell:            EC, radius 10 µm (bovine aortic endothelial cells, BAECs)
  - SPION loading:   sweep [10, 50, 100, 200, 300 pg]; 200 pg is Polyak's chosen dose
  - SPION chi:       2.0 (same as all variants; appropriate at 0.1 T per Polyak Fig. 2g)
  - Blood flow:      MCA-representative Poiseuille (0.2 m/s) for primary runs;
                     velocity sweep at 200 pg for context

Key purpose: establish how much the code's standard B₀ = 0.5 T overestimates
capture distances relative to Polyak's experimental field, and present a
side-by-side numerical comparison.

Outputs to results-polyak/ (new folder).
"""

import sys
import time
from pathlib import Path

proj_root = Path(__file__).parent.parent
sys.path.insert(0, str(proj_root))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from stent_capture.figures.common import DEFAULTS, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force
from stent_capture.physics.hydrodynamics import BloodFlow, stokes_drag
from stent_capture.physics.capture_criterion import capture_distance
from stent_capture.simulation.trajectories import integrate_trajectory

# ── Output directory ─────────────────────────────────────────────────────────
OUT_DIR = proj_root / "results-polyak"
OUT_DIR.mkdir(exist_ok=True)

# ── Polyak experimental parameters ───────────────────────────────────────────
_B0_POLYAK = 0.1        # T  — Polyak 2008: 1,000 G homogeneous field
_B0_CODE   = 0.5        # T  — code standard for comparison

_R_VES     = 1.54e-3    # m  — vessel radius (same geometry)
_V_MCA     = 0.2        # m/s — MCA representative velocity
_CELL_R    = 10e-6      # m  — EC radius (BAECs, same as code default)

# Loading sweep: includes Polyak's 200 pg dose + surrounding range
_LOADINGS_PG  = np.array([10., 50., 100., 200., 300.])   # pg
_VELOCITIES   = np.array([0.02, 0.05, 0.10, 0.20, 0.50])  # m/s

_N_ITER = 10   # binary search depth (~1.4 µm resolution)
_TRAJ_KW = dict(z_end=2e-3, max_time=1.5, rtol=1e-5, atol=1e-8)


# ── Known 0.5 T results (from existing code runs) for side-by-side comparison ─
# Source: results-200pg/COMPARISON_50pg_vs_200pg.md
_CODE_05T_LOADING = {   # (static_um, trajectory_um) at v=0.2 m/s, B0=0.5 T
    10:  (0.0,   51.4),
    50:  (6.8,   96.9),
    100: (24.4, 119.7),
    200: (39.0, 142.5),
    300: (47.7, 153.9),
}
_CODE_05T_VELOCITY = {   # trajectory_um at 200 pg, B0=0.5 T
    0.02: (94.5,  256.5),
    0.05: (71.1,  199.5),
    0.10: (56.5,  165.3),
    0.20: (39.0,  142.5),
    0.50: (18.5,  108.3),
}


# ── Helper: cell factory ──────────────────────────────────────────────────────
def _make_ec_cell(mass_pg):
    return SPIONLabelledCell(
        radius=_CELL_R,
        spion_mass_per_cell=mass_pg * 1e-15,
    )


# ── Helper: trajectory binary search ─────────────────────────────────────────
def _find_traj_capture(cell, tf, flow, ring, n_iter=_N_ITER):
    r_inner = ring.R - ring.t / 2
    lo, hi = 0.0, r_inner * 0.95

    traj_lo = integrate_trajectory(
        cell, tf, flow, ring,
        np.array([r_inner - 5e-6, 0.0, -2e-3]),
        **_TRAJ_KW,
    )
    if traj_lo.status != "captured":
        return 0.0

    for _ in range(n_iter):
        mid = (lo + hi) / 2.0
        traj = integrate_trajectory(
            cell, tf, flow, ring,
            np.array([r_inner - mid, 0.0, -2e-3]),
            **_TRAJ_KW,
        )
        if traj.status == "captured":
            lo = mid
        else:
            hi = mid
    return lo   # metres


# ─────────────────────────────────────────────────────────────────────────────
# Fig A: Force vs Radial Distance — 0.1 T vs 0.5 T
# ─────────────────────────────────────────────────────────────────────────────
def fig_A_force_comparison():
    print("  Fig A — Force vs distance (0.1 T vs 0.5 T)...", flush=True)

    ring_p = make_ring(); ring_p.assume_saturation = True
    ring_c = make_ring(); ring_c.assume_saturation = True

    tf_polyak = TotalField(ring_p, UniformExternalField([0.0, 0.0, _B0_POLYAK]))
    tf_code   = TotalField(ring_c, UniformExternalField([0.0, 0.0, _B0_CODE]))

    # Sample along +x axis from near stent inner surface to 300 µm into lumen
    r_inner = DEFAULTS["R"] - DEFAULTS["t"] / 2   # ≈ 1.46 mm
    r_vals = np.linspace(r_inner - 300e-6, r_inner - 2e-6, 200)
    d_vals = (r_inner - r_vals) * 1e6   # µm from inner stent surface, increasing inward
    pts = np.column_stack([r_vals, np.zeros_like(r_vals), np.zeros_like(r_vals)])

    cell_200 = _make_ec_cell(200)
    cell_10  = _make_ec_cell(10)

    F_200_p = np.linalg.norm(magnetic_force(cell_200, tf_polyak, pts), axis=1) * 1e12
    F_200_c = np.linalg.norm(magnetic_force(cell_200, tf_code,   pts), axis=1) * 1e12
    F_10_p  = np.linalg.norm(magnetic_force(cell_10,  tf_polyak, pts), axis=1) * 1e12
    F_10_c  = np.linalg.norm(magnetic_force(cell_10,  tf_code,   pts), axis=1) * 1e12

    # Stokes drag at v_mean = 0.2 m/s at r ≈ r_inner (near-wall)
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
    F_drag_200 = np.linalg.norm(stokes_drag(cell_200, flow, pts), axis=1) * 1e12

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(d_vals, F_200_c, 'r-',  lw=2.5, label="200 pg, B₀ = 0.5 T (code standard)")
    ax.semilogy(d_vals, F_200_p, 'r--', lw=2.5, label="200 pg, B₀ = 0.1 T (Polyak 2008)", alpha=0.8)
    ax.semilogy(d_vals, F_10_c,  'b-',  lw=1.5, label="10 pg,  B₀ = 0.5 T (code default)", alpha=0.7)
    ax.semilogy(d_vals, F_10_p,  'b--', lw=1.5, label="10 pg,  B₀ = 0.1 T (Polyak field)",  alpha=0.7)
    ax.semilogy(d_vals, F_drag_200, 'k:', lw=2.0, label="Stokes drag (v = 0.2 m/s, 200 pg cell)")

    ax.set_xlabel("Distance from stent inner surface (µm)", fontsize=11)
    ax.set_ylabel("Force magnitude (pN)", fontsize=11)
    ax.set_title(
        "Fig A — Magnetic force vs distance: B₀ = 0.1 T (Polyak) vs 0.5 T (code)\n"
        "EC cell (r = 10 µm), z = 0, θ = 0 (aligned with strut)", fontsize=11
    )
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(0, 300)

    fig.savefig(OUT_DIR / "figA_force_comparison.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figA_force_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved figA", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig B: Loading Sweep — 0.1 T vs 0.5 T (static + trajectory)
# ─────────────────────────────────────────────────────────────────────────────
def fig_B_loading_sweep():
    print("  Fig B — Loading sweep at 0.1 T (static + trajectory)...", flush=True)

    ring = make_ring(); ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_POLYAK]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)

    d_static_01 = []
    d_traj_01   = []

    for m_pg in _LOADINGS_PG:
        cell = _make_ec_cell(m_pg)
        ds = capture_distance(cell, tf, flow, direction="inward") * 1e6
        dt = _find_traj_capture(cell, tf, flow, ring) * 1e6
        d_static_01.append(ds)
        d_traj_01.append(dt)
        print(f"    {m_pg:5.0f} pg: static={ds:6.1f} µm  traj={dt:6.1f} µm", flush=True)

    d_static_01 = np.array(d_static_01)
    d_traj_01   = np.array(d_traj_01)

    # Pull 0.5 T values for the same loadings (from stored results)
    d_static_05 = np.array([_CODE_05T_LOADING.get(m, (np.nan, np.nan))[0] for m in _LOADINGS_PG])
    d_traj_05   = np.array([_CODE_05T_LOADING.get(m, (np.nan, np.nan))[1] for m in _LOADINGS_PG])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    fig.suptitle(
        "Fig B — EC loading sweep: B₀ = 0.1 T (Polyak) vs 0.5 T (code)  |  v = 0.2 m/s",
        fontsize=12
    )

    # Panel (a): Static criterion
    ax = axes[0]
    ax.plot(_LOADINGS_PG, d_static_01, 'o-', color="navy",    lw=2, ms=7, label="0.1 T (Polyak field)")
    ax.plot(_LOADINGS_PG, d_static_05, 's--', color="tomato", lw=2, ms=7, label="0.5 T (code standard)")
    for m, d in zip(_LOADINGS_PG, d_static_01):
        if d <= 0.5:
            ax.annotate("0 µm", xy=(m, 0.5), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=8, color="navy")
    ax.set_xlabel("SPION loading (pg/cell)", fontsize=11)
    ax.set_ylabel("Static capture distance (µm)", fontsize=11)
    ax.set_title("(a) Static criterion (|F_mag| > |F_drag|)", fontsize=10)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    # Panel (b): Trajectory
    ax = axes[1]
    ax.plot(_LOADINGS_PG, d_traj_01, 'o-', color="navy",    lw=2, ms=7, label="0.1 T (Polyak field)")
    ax.plot(_LOADINGS_PG, d_traj_05, 's--', color="tomato", lw=2, ms=7, label="0.5 T (code standard)")
    for m, d01, d05 in zip(_LOADINGS_PG, d_traj_01, d_traj_05):
        if not np.isnan(d05):
            ax.annotate(f"×{d05/max(d01,0.1):.1f}", xy=(m, (d01+d05)/2),
                        xytext=(5, 0), textcoords="offset points",
                        fontsize=7, color="gray", va="center")
    ax.axvline(200, color="steelblue", ls=":", lw=1.5, alpha=0.7, label="200 pg (Polyak dose)")
    ax.set_xlabel("SPION loading (pg/cell)", fontsize=11)
    ax.set_ylabel("Trajectory capture distance (µm)", fontsize=11)
    ax.set_title("(b) Trajectory integration (binary search, 10 iter)", fontsize=10)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "figB_loading_sweep_polyak_vs_code.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figB_loading_sweep_polyak_vs_code.pdf", bbox_inches="tight")
    plt.close(fig)

    return d_static_01, d_traj_01


# ─────────────────────────────────────────────────────────────────────────────
# Fig C: Velocity Sweep at 200 pg — 0.1 T vs 0.5 T
# ─────────────────────────────────────────────────────────────────────────────
def fig_C_velocity_sweep():
    print("  Fig C — Velocity sweep at 200 pg, 0.1 T...", flush=True)

    ring = make_ring(); ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_POLYAK]))
    cell = _make_ec_cell(200)

    d_static_01 = []
    d_traj_01   = []

    for v in _VELOCITIES:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        ds = capture_distance(cell, tf, flow, direction="inward") * 1e6
        dt = _find_traj_capture(cell, tf, flow, ring) * 1e6
        d_static_01.append(ds)
        d_traj_01.append(dt)
        print(f"    v={v:.3f} m/s: static={ds:6.1f} µm  traj={dt:6.1f} µm", flush=True)

    d_static_01 = np.array(d_static_01)
    d_traj_01   = np.array(d_traj_01)

    d_static_05 = np.array([_CODE_05T_VELOCITY.get(v, (np.nan, np.nan))[0] for v in _VELOCITIES])
    d_traj_05   = np.array([_CODE_05T_VELOCITY.get(v, (np.nan, np.nan))[1] for v in _VELOCITIES])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        "Fig C — EC velocity sweep at 200 pg: B₀ = 0.1 T (Polyak) vs 0.5 T (code)",
        fontsize=12
    )

    labels = [f"{v:.2f}" for v in _VELOCITIES]

    for ax, d01, d05, title in [
        (axes[0], d_static_01, d_static_05, "(a) Static criterion"),
        (axes[1], d_traj_01,   d_traj_05,   "(b) Trajectory integration"),
    ]:
        ax.plot(_VELOCITIES, d01, 'o-', color="navy",    lw=2, ms=7, label="0.1 T (Polyak field)")
        ax.plot(_VELOCITIES, d05, 's--', color="tomato", lw=2, ms=7, label="0.5 T (code standard)")
        ax.axvline(_V_MCA, color="gray", ls=":", lw=1.2, alpha=0.6, label="v = 0.2 m/s (MCA)")
        ax.set_xscale("log")
        ax.set_xlabel("Mean blood velocity (m/s)", fontsize=11)
        ax.set_ylabel("Capture distance (µm)", fontsize=11)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=10)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "figC_velocity_sweep_polyak_vs_code.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figC_velocity_sweep_polyak_vs_code.pdf", bbox_inches="tight")
    plt.close(fig)

    return d_static_01, d_traj_01


# ─────────────────────────────────────────────────────────────────────────────
# Numerical summary report
# ─────────────────────────────────────────────────────────────────────────────
def write_report(loading_static, loading_traj, vel_static, vel_traj):
    lines = []
    lines.append("# Polyak (2008) Comparison Report — B₀ = 0.1 T")
    lines.append("")
    lines.append(f"**Generated:** 2026-04-12")
    lines.append(f"**Cell:** EC, radius = 10 µm (BAEC representative)")
    lines.append(f"**B₀ (this run):** 0.1 T  (1,000 G — Polyak 2008 experimental)")
    lines.append(f"**B₀ (code standard):** 0.5 T  (5× stronger)")
    lines.append(f"**Binary search depth:** {_N_ITER} iterations (~1.4 µm precision)")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Loading Sweep at v = 0.2 m/s")
    lines.append("")
    lines.append("| Loading (pg) | Static 0.1T (µm) | Static 0.5T (µm) | Traj 0.1T (µm) | Traj 0.5T (µm) | Traj ratio 0.5T/0.1T |")
    lines.append("|---|---|---|---|---|---|")

    for i, m in enumerate(_LOADINGS_PG):
        s01 = loading_static[i]
        t01 = loading_traj[i]
        s05, t05 = _CODE_05T_LOADING.get(int(m), (np.nan, np.nan))
        ratio = t05 / max(t01, 0.1) if not np.isnan(t05) and t01 > 0 else float('nan')
        ratio_str = f"{ratio:.2f}×" if not np.isnan(ratio) else "N/A"
        lines.append(f"| {m:.0f} | {s01:.1f} | {s05:.1f} | {t01:.1f} | {t05:.1f} | {ratio_str} |")

    lines.append("")
    lines.append("**Key finding:** Code (0.5 T) overestimates trajectory capture distance relative to Polyak's experimental field (0.1 T).")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Velocity Sweep at 200 pg (Polyak dose)")
    lines.append("")
    lines.append("| Velocity (m/s) | Static 0.1T (µm) | Static 0.5T (µm) | Traj 0.1T (µm) | Traj 0.5T (µm) | Traj ratio 0.5T/0.1T |")
    lines.append("|---|---|---|---|---|---|")

    for i, v in enumerate(_VELOCITIES):
        s01 = vel_static[i]
        t01 = vel_traj[i]
        s05, t05 = _CODE_05T_VELOCITY.get(v, (np.nan, np.nan))
        ratio = t05 / max(t01, 0.1) if not np.isnan(t05) and t01 > 0 else float('nan')
        ratio_str = f"{ratio:.2f}×" if not np.isnan(ratio) else "N/A"
        lines.append(f"| {v:.3f} | {s01:.1f} | {s05:.1f} | {t01:.1f} | {t05:.1f} | {ratio_str} |")

    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("### Why the results differ")
    lines.append("")
    lines.append("The magnetic force is:")
    lines.append("```")
    lines.append("F = (V_spion × χ / μ₀) × |B_total| × |∇|B_total||")
    lines.append("```")
    lines.append("")
    lines.append("At distances far from the strut, |B_total| ≈ B₀. The 5× increase in B₀")
    lines.append("(0.1 → 0.5 T) increases |B_total| by ~5×. The gradient |∇|B_total|| also")
    lines.append("changes because B_total direction rotates. The net effect on capture distance")
    lines.append("is visible in the tables above.")
    lines.append("")
    lines.append("### Stent saturation note")
    lines.append("")
    lines.append("Both runs use `assume_saturation = True` (M = 1.0×10⁶ A/m). At 0.1 T,")
    lines.append("Polyak reports ~80–90% saturation for 304 SS — so the stent gradient field")
    lines.append("is slightly overestimated (~10–20%) in the 0.1 T run. This means the 0.1 T")
    lines.append("results here are a slight overestimate of Polyak's true experimental forces.")
    lines.append("")
    lines.append("### SPION saturation note")
    lines.append("")
    lines.append("χ = 2.0 is used at both B₀ values. Polyak reports ~80–90% SPION saturation")
    lines.append("at 0.1 T. At 0.5 T (code standard), SPIONs are fully saturated, making the")
    lines.append("linear χ approximation an overestimate — so the code's 0.5 T force is")
    lines.append("doubly overestimated (higher B₀ AND fixed χ rather than saturated M_sat).")
    lines.append("")
    lines.append("### Capture efficiency (Polyak: 20%)")
    lines.append("")
    lines.append("Polyak measured 20% capture of 2.5×10⁶ cells at 0.1 T, 200 pg, 30 ml/min.")
    lines.append("This code computes capture DISTANCE (µm), not efficiency (%). These cannot")
    lines.append("be directly compared without the cell spatial distribution and stent area.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Files Generated")
    lines.append("")
    lines.append("- `figA_force_comparison.{png,pdf}` — Force vs distance at both field strengths")
    lines.append("- `figB_loading_sweep_polyak_vs_code.{png,pdf}` — Loading sweep comparison")
    lines.append("- `figC_velocity_sweep_polyak_vs_code.{png,pdf}` — Velocity sweep comparison")

    report_path = OUT_DIR / "POLYAK_COMPARISON_REPORT.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"    Report saved to {report_path}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    t0 = time.time()
    print("=" * 60)
    print("Polyak (2008) Comparison -- B0 = 0.1 T")
    print("=" * 60)

    print("\n[Fig A] Force vs distance...")
    fig_A_force_comparison()

    print("\n[Fig B] Loading sweep at 0.1 T (v = 0.2 m/s)...")
    load_static, load_traj = fig_B_loading_sweep()

    print("\n[Fig C] Velocity sweep at 200 pg, 0.1 T...")
    vel_static, vel_traj = fig_C_velocity_sweep()

    print("\n[Report] Writing summary...")
    write_report(load_static, load_traj, vel_static, vel_traj)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed/60:.1f} min -- results in results-polyak/")
    print("=" * 60)
