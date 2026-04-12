"""
Generate all figures with 200 pg SPION loading (vs default 10 pg).
Output to results-200pg/ for comparison.
"""

import sys
import time
from pathlib import Path

# Add project root to path
proj_root = Path(__file__).parent.parent
sys.path.insert(0, str(proj_root))

import numpy as np
import matplotlib.pyplot as plt

from stent_capture.figures.common import make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.physics.capture_criterion import capture_distance
from stent_capture.simulation.trajectories import integrate_trajectory

# Output directory
OUT_DIR = proj_root / "results-200pg"
OUT_DIR.mkdir(exist_ok=True)

# Parameters (same as defaults except SPION mass)
_B0_Z = 0.5
_R_VES = 1.54e-3
_V_MCA = 0.2
_M_REF_PG = 200.0  # CHANGED: was 50 pg, now 200 pg
_N_ITER = 7
_SPION_MASS_PG = 200.0  # Fixed SPION mass for all figs

_LOADINGS_PG = np.array([10., 30., 50., 100., 200., 300., 400.])
_VELOCITIES = np.array([0.02, 0.05, 0.10, 0.20, 0.50])

_TRAJ_KW = dict(z_end=2e-3, max_time=1.5, rtol=1e-5, atol=1e-8)

_COLOR_STATIC = "#2980b9"
_COLOR_TRAJ = "#c0392b"


# ─────────────────────────────────────────────────────────────────────────────
# Fig 16: Capture Map (Stage 2 — static criterion with 200 pg)
# ─────────────────────────────────────────────────────────────────────────────

def fig16_200pg():
    """Static capture map with 200 pg SPION loading."""
    print("  Generating fig16 (capture map, 200 pg)...")

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)

    cell = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)

    # Capture map on r-z plane (theta=0)
    r_pts = np.linspace(1.4e-3, 1.55e-3, 100)
    z_pts = np.linspace(-3e-3, 3e-3, 150)
    R, Z = np.meshgrid(r_pts, z_pts)
    X = R * 1.0  # theta=0
    Y = R * 0.0

    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    from stent_capture.physics.magnetic_force import magnetic_force
    from stent_capture.core.gradient import compute_gradient_magnitude

    F_mag_vec = magnetic_force(cell, tf, pts)
    F_mag = np.linalg.norm(F_mag_vec, axis=1)

    from stent_capture.physics.hydrodynamics import stokes_drag
    F_drag = np.linalg.norm(stokes_drag(cell, flow, pts), axis=1)

    captured = (F_mag > F_drag).astype(float).reshape(R.shape)

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.contourf(R*1e3, Z*1e3, captured, levels=[0, 0.5, 1],
                     colors=["#ecf0f1", "#27ae60"], alpha=0.8)
    ax.contour(R*1e3, Z*1e3, captured, levels=[0.5], colors="black", linewidths=1.5)

    ax.set_xlabel("Radial distance from axis (mm)", fontsize=11)
    ax.set_ylabel("Axial distance (mm)", fontsize=11)
    ax.set_title(f"Fig 16 (modified) — Static capture criterion (200 pg SPION)\n"
                f"v = {_V_MCA} m/s, B₀ = {_B0_Z} T, strut-aligned axis (θ=0)",
                fontsize=12)
    ax.grid(True, alpha=0.3)

    fig.savefig(OUT_DIR / "fig16_capture_map_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig16_capture_map_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig16 to {OUT_DIR}/fig16_capture_map_200pg.*")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 21: Static vs Trajectory (Stage 3 — headline result with 200 pg)
# ─────────────────────────────────────────────────────────────────────────────

def _find_range_trajectory_200pg(cell, tf, flow, ring, label="", **kw):
    """Binary search for trajectory capture range."""
    r_inner = ring.R - ring.t / 2
    n = [0]

    def _check(d):
        n[0] += 1
        r = r_inner - max(d, 1e-7)
        traj = integrate_trajectory(
            cell, tf, flow, ring,
            np.array([r, 0.0, -2e-3]),
            **kw,
        )
        print(
            f"      [{label}] eval {n[0]:2d}: d = {d*1e6:7.1f} um  ->  {traj.status}",
            flush=True,
        )
        return traj.status

    d_lo = 1e-7
    d_hi = r_inner * 0.999

    if _check(d_lo) != "captured":
        return 0.0
    if _check(d_hi) == "captured":
        return float(r_inner)

    for _ in range(_N_ITER):
        d_mid = (d_lo + d_hi) / 2
        if _check(d_mid) == "captured":
            d_lo = d_mid
        else:
            d_hi = d_mid

    return (d_lo + d_hi) / 2


def fig21_200pg():
    """Static vs trajectory comparison with 200 pg SPION (fixed, not swept)."""
    print("  Generating fig21 (static vs trajectory, 200 pg fixed)...")

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    # Panel (a): loading sweep at fixed velocity
    print("    Panel (a): loading sweep...")
    cell_base = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)

    static_a, traj_a = [], []
    for m_pg in _LOADINGS_PG:
        cell = SPIONLabelledCell(spion_mass_per_cell=m_pg * 1e-15)
        print(f"      loading = {m_pg:.0f} pg")
        d_s = capture_distance(cell, tf, flow, direction="inward") * 1e6
        static_a.append(d_s)
        print(f"        static d = {d_s:.1f} um")
        d_t = _find_range_trajectory_200pg(
            cell, tf, flow, ring, label=f"{m_pg:.0f}pg", **_TRAJ_KW,
        ) * 1e6
        traj_a.append(d_t)
        print(f"        traj range = {d_t:.1f} um")

    # Panel (b): velocity sweep at fixed 200 pg
    print("    Panel (b): velocity sweep at 200 pg...")
    static_b, traj_b = [], []
    cell_200pg = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)

    for v in _VELOCITIES:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        print(f"      v_mean = {v:.3f} m/s")
        d_s = capture_distance(cell_200pg, tf, flow, direction="inward") * 1e6
        static_b.append(d_s)
        print(f"        static d = {d_s:.1f} um")
        d_t = _find_range_trajectory_200pg(
            cell_200pg, tf, flow, ring, label=f"{v:.3f}m/s", **_TRAJ_KW,
        ) * 1e6
        traj_b.append(d_t)
        print(f"        traj range = {d_t:.1f} um")

    # Reference values (200 pg, 0.2 m/s)
    idx_ref = int(np.where(_LOADINGS_PG == 200.0)[0][0])
    val_static_ref = static_a[idx_ref]
    val_traj_ref = traj_a[idx_ref]
    ratio_str = (
        f"{val_traj_ref / val_static_ref:.1f}x"
        if val_static_ref > 0.5 else ">> static (~0 um)"
    )

    # Build figure
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.subplots_adjust(wspace=0.30, left=0.09, right=0.97, top=0.82, bottom=0.13)
    max_y = max(max(traj_a), max(traj_b), 50) * 1.15

    # Panel (a)
    ax_a.semilogx(_LOADINGS_PG, static_a, "-o",
                  color=_COLOR_STATIC, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Static  |F$_{mag}$| > |F$_{drag}$|", zorder=5)
    ax_a.semilogx(_LOADINGS_PG, traj_a, "--s",
                  color=_COLOR_TRAJ, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Trajectory  (strut-aligned, binary search)", zorder=5)
    ax_a.fill_between(_LOADINGS_PG, static_a, traj_a, alpha=0.10, color=_COLOR_TRAJ)

    for m, tr in zip(_LOADINGS_PG, traj_a):
        ax_a.annotate(f"{tr:.0f}", xy=(m, tr), xytext=(0, 6),
                      textcoords="offset points", ha="center", fontsize=7, color=_COLOR_TRAJ)

    ax_a.axvline(200, color="red", ls="--", lw=1.5, alpha=0.8,
                 label="200 pg (this run)")
    ax_a.axvline(50, color="gray", ls="--", lw=1.0, alpha=0.7,
                 label="50 pg (project default)")
    ax_a.set_xlabel("SPION loading per cell (pg)", fontsize=10)
    ax_a.set_ylabel("Effective capture range ($\\mu$m)", fontsize=10)
    ax_a.set_xlim(7, 500)
    ax_a.set_ylim(0, max_y)
    ax_a.legend(fontsize=7, loc="upper left")
    ax_a.grid(True, which="both", alpha=0.25)
    ax_a.set_title(
        f"(a) Capture range vs SPION loading\n"
        f"$\\bar{{v}}$ = {_V_MCA:.1f} m/s (MCA mean), B$_0$ = {_B0_Z} T",
        fontsize=10,
    )

    # Panel (b)
    ax_b.semilogx(_VELOCITIES, static_b, "-o",
                  color=_COLOR_STATIC, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Static  |F$_{mag}$| > |F$_{drag}$|", zorder=5)
    ax_b.semilogx(_VELOCITIES, traj_b, "--s",
                  color=_COLOR_TRAJ, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Trajectory  (strut-aligned, binary search)", zorder=5)
    ax_b.fill_between(_VELOCITIES, static_b, traj_b, alpha=0.10, color=_COLOR_TRAJ)

    for v, tr in zip(_VELOCITIES, traj_b):
        ax_b.annotate(f"{tr:.0f}", xy=(v, tr), xytext=(0, 6),
                      textcoords="offset points", ha="center", fontsize=7, color=_COLOR_TRAJ)

    ax_b.axvline(0.20, color="gray", ls="--", lw=1.0, alpha=0.7,
                 label="0.20 m/s (MCA mean)")
    ax_b.set_xlabel("Mean blood velocity $\\bar{v}$ (m/s)", fontsize=10)
    ax_b.set_ylabel("Effective capture range ($\\mu$m)", fontsize=10)
    ax_b.set_xlim(0.014, 0.65)
    ax_b.set_ylim(0, max_y)
    ax_b.legend(fontsize=7, loc="upper right")
    ax_b.grid(True, which="both", alpha=0.25)
    ax_b.set_title(
        f"(b) Capture range vs flow velocity\n"
        f"200 pg SPION loading, B$_0$ = {_B0_Z} T",
        fontsize=10,
    )

    # Suptitle
    fig.suptitle(
        "Fig 21 (200 pg variant) — Static vs trajectory capture predictions\n"
        "Solid blue: Furlani & Ng (2006) static criterion  |  "
        "Dashed red: trajectory effective range (binary search over injection radius)\n"
        f"At 200 pg loading with MCA velocity ({_V_MCA:.1f} m/s, {_B0_Z} T):  "
        f"trajectory = {val_traj_ref:.0f} µm  vs  static = {val_static_ref:.0f} µm  "
        f"({ratio_str} larger)\n"
        "Radial drift over 2 mm approach extends capture zone beyond instantaneous force balance",
        fontsize=8.5, y=0.995,
    )

    fig.savefig(OUT_DIR / "fig21_static_vs_trajectory_200pg.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig21_static_vs_trajectory_200pg.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"    Saved fig21 to {OUT_DIR}/fig21_static_vs_trajectory_200pg.*")

    # Print summary table
    print("\n  Panel (a) -- loading sweep, v = 0.2 m/s:")
    print(f"  {'loading (pg)':>14}  {'static (um)':>12}  {'traj (um)':>12}  {'ratio':>8}")
    for m, st, tr in zip(_LOADINGS_PG, static_a, traj_a):
        ratio = f"{tr/st:.1f}x" if st > 0.5 else "inf"
        print(f"  {m:>14.0f}  {st:>12.1f}  {tr:>12.1f}  {ratio:>8}")

    print("\n  Panel (b) -- velocity sweep, 200 pg:")
    print(f"  {'v_mean (m/s)':>14}  {'static (um)':>12}  {'traj (um)':>12}  {'ratio':>8}")
    for v, st, tr in zip(_VELOCITIES, static_b, traj_b):
        ratio = f"{tr/st:.1f}x" if st > 0.5 else "inf"
        print(f"  {v:>14.3f}  {st:>12.1f}  {tr:>12.1f}  {ratio:>8}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"\n{'='*80}")
    print("GENERATING RESULTS WITH 200 PG SPION LOADING")
    print(f"Output directory: {OUT_DIR}")
    print(f"{'='*80}\n")

    t_start = time.time()

    fig16_200pg()
    print()
    fig21_200pg()

    elapsed = time.time() - t_start
    print(f"\n{'='*80}")
    print(f"All figures generated successfully in {elapsed:.1f} s")
    print(f"Results saved to: {OUT_DIR}")
    print(f"{'='*80}\n")
