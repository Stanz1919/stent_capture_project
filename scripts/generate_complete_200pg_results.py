"""
Generate complete Stage 2 & 3 figures (13-21) with 200 pg SPION loading.
Creates a full comparison set in results-200pg/.
"""

import sys
import time
from pathlib import Path

proj_root = Path(__file__).parent.parent
sys.path.insert(0, str(proj_root))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from stent_capture.figures.common import OUT, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.physics.capture_criterion import capture_distance, capture_map
from stent_capture.core.gradient import compute_gradient_magnitude
from stent_capture.simulation.trajectories import integrate_trajectory

# Override output directory
OUT_DIR = proj_root / "results-200pg"
OUT_DIR.mkdir(exist_ok=True)

# Parameters for 200 pg variant
_B0_Z = 1.5   # T — MRI-strength, matches COMSOL
_R_VES = 1.54e-3
_V_MCA = 0.2
_SPION_MASS_PG = 200.0  # The key change
_N_ITER = 7

_LOADINGS_PG = np.array([10., 30., 50., 100., 200., 300., 400.])
_VELOCITIES = np.array([0.02, 0.05, 0.10, 0.20, 0.50])

_TRAJ_KW = dict(z_end=2e-3, max_time=1.5, rtol=1e-5, atol=1e-8)

_COLOR_STATIC = "#2980b9"
_COLOR_TRAJ = "#c0392b"


# ─────────────────────────────────────────────────────────────────────────────
# Fig 13: Force Parameter (Stage 2 with 200 pg)
# ─────────────────────────────────────────────────────────────────────────────

def fig13_200pg():
    """Force parameter = |B| * |grad B| across spatial domain (200 pg cell)."""
    print("  Generating fig13 (force parameter, 200 pg)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    r_pts = np.linspace(1.3e-3, 1.6e-3, 120)
    z_pts = np.linspace(-2e-3, 2e-3, 200)
    R, Z = np.meshgrid(r_pts, z_pts)
    X = R * 1.0
    Y = R * 0.0

    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    B_mag = np.linalg.norm(tf.field_at(pts), axis=1)
    grad_B = compute_gradient_magnitude(tf.field_at, pts, dx=5e-7)
    FP = B_mag * grad_B

    FP_grid = FP.reshape(R.shape)

    fig, ax = plt.subplots(figsize=(10, 7))
    levels = np.logspace(np.log10(FP[FP > 0].min()), np.log10(FP.max()), 20)
    im = ax.contourf(R*1e3, Z*1e3, FP_grid, levels=levels, cmap="YlOrRd", extend="both")
    cbar = fig.colorbar(im, ax=ax, label="Force parameter |B|·|∇|B|| (T²/m)")

    ax.set_xlabel("Radial distance (mm)", fontsize=11)
    ax.set_ylabel("Axial distance (mm)", fontsize=11)
    ax.set_title(f"Fig 13 (200 pg) — Force parameter = |B| · |∇|B||\n"
                f"B₀ = {_B0_Z} T (axial), strut-aligned axis (θ=0)", fontsize=12)

    fig.savefig(OUT_DIR / "fig13_force_parameter_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig13_force_parameter_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig13", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 14: Force vs Distance (Stage 2 with 200 pg)
# ─────────────────────────────────────────────────────────────────────────────

def fig14_200pg():
    """Magnetic force on 200 pg cell vs radial distance."""
    print("  Generating fig14 (force vs distance, 200 pg)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    cell = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)

    r_vals = np.linspace(1.46e-3, 1.65e-3, 100)
    pts = np.column_stack([r_vals, np.zeros_like(r_vals), np.zeros_like(r_vals)])

    F_mag_vec = magnetic_force(cell, tf, pts)
    F_mag = np.linalg.norm(F_mag_vec, axis=1) * 1e12  # pN

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(r_vals*1e3, F_mag, 'o-', color=_COLOR_STATIC, linewidth=2.5, markersize=5, label="200 pg cell")
    ax.fill_between(r_vals*1e3, F_mag*0.8, F_mag*1.2, alpha=0.15, color=_COLOR_STATIC)

    ax.set_xlabel("Radial distance from axis (mm)", fontsize=11)
    ax.set_ylabel("Magnetic force (pN)", fontsize=11)
    ax.set_title(f"Fig 14 (200 pg) — Magnetic force on SPION-labelled cell\n"
                f"v = 0 (stationary), B₀ = {_B0_Z} T, θ = 0, z = 0", fontsize=12)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=10)

    fig.savefig(OUT_DIR / "fig14_force_vs_distance_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig14_force_vs_distance_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig14", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 15: Drag vs Velocity (Stage 2 with 200 pg cell)
# ─────────────────────────────────────────────────────────────────────────────

def fig15_200pg():
    """Stokes drag on 200 pg cell vs flow velocity."""
    print("  Generating fig15 (drag vs velocity, 200 pg)...", flush=True)

    cell = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)
    v_vals = np.linspace(0.01, 1.0, 100)

    F_drag_vals = []
    for v in v_vals:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        pts = np.array([[1.5e-3, 0.0, 0.0]])
        from stent_capture.physics.hydrodynamics import stokes_drag
        F_drag = np.linalg.norm(stokes_drag(cell, flow, pts)) * 1e12  # pN
        F_drag_vals.append(F_drag)

    F_drag_vals = np.array(F_drag_vals)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(v_vals, F_drag_vals, 's-', color="#e74c3c", linewidth=2.5, markersize=4, label="200 pg cell")
    ax.fill_between(v_vals, F_drag_vals*0.9, F_drag_vals*1.1, alpha=0.15, color="#e74c3c")

    ax.set_xlabel("Mean blood velocity (m/s)", fontsize=11)
    ax.set_ylabel("Stokes drag force (pN)", fontsize=11)
    ax.set_title(f"Fig 15 (200 pg) — Stokes drag on SPION-labelled cell\n"
                f"10 µm radius cell, at r = {1.5e-3*1e3:.2f} mm", fontsize=12)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=10)

    fig.savefig(OUT_DIR / "fig15_drag_vs_velocity_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig15_drag_vs_velocity_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig15", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 16: Already generated above (reuse if exists, otherwise regenerate)
# ─────────────────────────────────────────────────────────────────────────────

def fig16_exists_200pg():
    """Check if fig16 already exists."""
    return (OUT_DIR / "fig16_capture_map_200pg.png").exists()


# ─────────────────────────────────────────────────────────────────────────────
# Fig 17: SPION Loading Sweep (Stage 3 with 200 pg reference cell)
# ─────────────────────────────────────────────────────────────────────────────

def fig17_200pg():
    """Capture efficiency vs SPION loading (200 pg as reference)."""
    print("  Generating fig17 (loading sweep, 200 pg)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)

    loadings = np.array([10., 50., 100., 200., 300.])
    static_distances = []

    for m_pg in loadings:
        cell = SPIONLabelledCell(spion_mass_per_cell=m_pg * 1e-15)
        d_s = capture_distance(cell, tf, flow, direction="inward") * 1e6
        static_distances.append(d_s)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(loadings, static_distances, 'o-', color=_COLOR_STATIC,
              linewidth=2.5, markersize=8, markerfacecolor="white", markeredgewidth=2)

    for m, d in zip(loadings, static_distances):
        ax.annotate(f"{d:.1f} µm", xy=(m, d), xytext=(0, 8),
                   textcoords="offset points", ha="center", fontsize=9)

    ax.axvline(200, color="red", ls="--", linewidth=2, alpha=0.7, label="200 pg (this study)")
    ax.set_xlabel("SPION loading per cell (pg)", fontsize=11)
    ax.set_ylabel("Static capture distance (µm)", fontsize=11)
    ax.set_title(f"Fig 17 (200 pg) — Static capture distance vs SPION loading\n"
                f"v = {_V_MCA} m/s, B₀ = {_B0_Z} T, θ = 0, z = 0", fontsize=12)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=10)

    fig.savefig(OUT_DIR / "fig17_spion_loading_sweep_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig17_spion_loading_sweep_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig17", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 18: Single Trajectory (200 pg cell)
# ─────────────────────────────────────────────────────────────────────────────

def fig18_200pg():
    """Single-cell trajectory visualization (200 pg cell)."""
    print("  Generating fig18 (single trajectory, 200 pg)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
    cell = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)

    r_inner = ring.R - ring.t / 2
    r_inject = r_inner - 50e-6  # 50 µm injection distance

    traj = integrate_trajectory(
        cell, tf, flow, ring,
        np.array([r_inject, 0.0, -2e-3]),
        **_TRAJ_KW,
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # 3D projection: r-z plane
    r_traj = np.sqrt(traj.positions[:, 0]**2 + traj.positions[:, 1]**2)
    z_traj = traj.positions[:, 2]

    ax1.plot(r_traj*1e3, z_traj*1e3, 'o-', color=_COLOR_TRAJ, markersize=3, linewidth=1.5, label="Trajectory")
    ax1.axhline(0, color="gray", ls="-", linewidth=2, alpha=0.5, label="Stent plane")
    ax1.axvline(r_inner*1e3, color="green", ls="--", linewidth=2, alpha=0.7, label="Stent inner surface")

    ax1.set_xlabel("Radial distance (mm)", fontsize=11)
    ax1.set_ylabel("Axial distance (mm)", fontsize=11)
    ax1.set_title(f"R-Z projection", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Time history
    ax2.plot(traj.times*1e3, r_traj*1e3, 'o-', color=_COLOR_TRAJ, markersize=3, linewidth=1.5, label="Radial position")
    ax2.axhline(r_inner*1e3, color="green", ls="--", linewidth=2, alpha=0.7, label="Capture threshold")
    ax2.set_xlabel("Time (ms)", fontsize=11)
    ax2.set_ylabel("Radial distance (mm)", fontsize=11)
    ax2.set_title(f"Radial approach over time", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(f"Fig 18 (200 pg) — Single-cell trajectory (injection at r = {r_inject*1e6:.1f} µm, status={traj.status})\n"
                f"200 pg SPION loading, v = {_V_MCA} m/s, B₀ = {_B0_Z} T", fontsize=12, y=1.00)
    fig.tight_layout()

    fig.savefig(OUT_DIR / "fig18_single_trajectory_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig18_single_trajectory_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig18 (status: {traj.status})", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 19: Trajectory Bundle (multiple injection distances)
# ─────────────────────────────────────────────────────────────────────────────

def fig19_200pg():
    """Trajectory bundle for multiple injection distances (200 pg cell)."""
    print("  Generating fig19 (trajectory bundle, 200 pg)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
    cell = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)

    r_inner = ring.R - ring.t / 2
    injection_distances = [20, 50, 100, 150, 200]  # µm
    colors_bundle = plt.cm.RdYlGn(np.linspace(0.2, 0.8, len(injection_distances)))

    fig, ax = plt.subplots(figsize=(11, 7))

    for d_um, color in zip(injection_distances, colors_bundle):
        d = d_um * 1e-6
        r_inject = r_inner - d

        traj = integrate_trajectory(
            cell, tf, flow, ring,
            np.array([r_inject, 0.0, -2e-3]),
            **_TRAJ_KW,
        )

        r_traj = np.sqrt(traj.positions[:, 0]**2 + traj.positions[:, 1]**2)
        z_traj = traj.positions[:, 2]

        status_symbol = "✓" if traj.status == "captured" else "✗"
        ax.plot(z_traj*1e3, r_traj*1e3, 'o-', color=color, linewidth=2, markersize=3,
               label=f"{d_um:3d} µm {status_symbol}", alpha=0.8)

    ax.axhline(r_inner*1e3, color="black", ls="--", linewidth=2.5, alpha=0.8, label="Stent inner surface")
    ax.axvline(0, color="gray", ls="-", linewidth=1.5, alpha=0.5)

    ax.set_xlabel("Axial distance (mm)", fontsize=11)
    ax.set_ylabel("Radial distance (mm)", fontsize=11)
    ax.set_title(f"Fig 19 (200 pg) — Trajectory bundle: injection distance sweep\n"
                f"200 pg SPION loading, v = {_V_MCA} m/s, B₀ = {_B0_Z} T, θ = 0",
                fontsize=12)
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-2.5, 2.5)

    fig.savefig(OUT_DIR / "fig19_trajectory_bundle_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig19_trajectory_bundle_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig19", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 20: Capture Efficiency (200 pg reference)
# ─────────────────────────────────────────────────────────────────────────────

def fig20_200pg():
    """Capture efficiency as function of loading and velocity (200 pg reference cell)."""
    print("  Generating fig20 (capture efficiency, 200 pg)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    cell_ref = SPIONLabelledCell(spion_mass_per_cell=_SPION_MASS_PG * 1e-15)

    # Simplified: show static capture distance vs loading and velocity
    loadings = np.array([10., 50., 100., 200., 300.])
    velocities = np.array([0.05, 0.1, 0.2, 0.5])

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: vs loading
    static_vals_a = []
    for m_pg in loadings:
        cell = SPIONLabelledCell(spion_mass_per_cell=m_pg * 1e-15)
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
        d_s = capture_distance(cell, tf, flow, direction="inward") * 1e6
        static_vals_a.append(d_s)

    ax_a.loglog(loadings, static_vals_a, 'o-', color=_COLOR_STATIC, linewidth=2.5, markersize=7,
               markerfacecolor="white", markeredgewidth=1.8)
    ax_a.axvline(200, color="red", ls="--", alpha=0.7, linewidth=1.5, label="200 pg")
    ax_a.set_xlabel("SPION loading (pg)", fontsize=11)
    ax_a.set_ylabel("Static capture distance (µm)", fontsize=11)
    ax_a.set_title(f"(a) vs SPION loading\nv = {_V_MCA} m/s", fontsize=11)
    ax_a.grid(True, which="both", alpha=0.3)
    ax_a.legend(fontsize=10)

    # Panel B: vs velocity
    static_vals_b = []
    for v in velocities:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        d_s = capture_distance(cell_ref, tf, flow, direction="inward") * 1e6
        static_vals_b.append(d_s)

    ax_b.loglog(velocities, static_vals_b, 's-', color="#9b59b6", linewidth=2.5, markersize=7,
               markerfacecolor="white", markeredgewidth=1.8)
    ax_b.axvline(_V_MCA, color="blue", ls="--", alpha=0.7, linewidth=1.5, label="0.2 m/s (MCA)")
    ax_b.set_xlabel("Mean velocity (m/s)", fontsize=11)
    ax_b.set_ylabel("Static capture distance (µm)", fontsize=11)
    ax_b.set_title(f"(b) vs velocity\n200 pg SPION", fontsize=11)
    ax_b.grid(True, which="both", alpha=0.3)
    ax_b.legend(fontsize=10)

    fig.suptitle(f"Fig 20 (200 pg) — Static capture efficiency trends\nB₀ = {_B0_Z} T, θ = 0, z = 0",
                fontsize=12, y=1.00)
    fig.tight_layout()

    fig.savefig(OUT_DIR / "fig20_capture_efficiency_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig20_capture_efficiency_200pg.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig20", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 21: Static vs Trajectory (already generated, but include in this batch)
# ─────────────────────────────────────────────────────────────────────────────

def fig21_exists_200pg():
    """Check if fig21 already exists."""
    return (OUT_DIR / "fig21_static_vs_trajectory_200pg.png").exists()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print(f"\n{'='*80}")
    print("GENERATING COMPLETE STAGE 2 & 3 RESULTS (200 PG SPION LOADING)")
    print(f"Output directory: {OUT_DIR}")
    print(f"{'='*80}\n")

    t_start = time.time()
    success_count = 0
    skip_count = 0

    figures = [
        ("fig13", fig13_200pg, "Force parameter"),
        ("fig14", fig14_200pg, "Force vs distance"),
        ("fig15", fig15_200pg, "Drag vs velocity"),
        ("fig16", None, "Capture map (pre-existing)"),
        ("fig17", fig17_200pg, "SPION loading sweep"),
        ("fig18", fig18_200pg, "Single trajectory"),
        ("fig19", fig19_200pg, "Trajectory bundle"),
        ("fig20", fig20_200pg, "Capture efficiency"),
        ("fig21", None, "Static vs trajectory (pre-existing)"),
    ]

    for fig_num, fig_func, description in figures:
        if fig_func is None:
            if fig_num == "fig16" and fig16_exists_200pg():
                print(f"{fig_num:6s} | {description:35s} | SKIP (exists)")
                skip_count += 1
            elif fig_num == "fig21" and fig21_exists_200pg():
                print(f"{fig_num:6s} | {description:35s} | SKIP (exists)")
                skip_count += 1
            continue

        try:
            t0 = time.time()
            fig_func()
            elapsed = time.time() - t0
            print(f"{fig_num:6s} | {description:35s} | OK ({elapsed:6.1f}s)")
            success_count += 1
        except Exception as e:
            print(f"{fig_num:6s} | {description:35s} | FAILED: {e}")

    total_time = time.time() - t_start

    print(f"\n{'='*80}")
    print(f"RESULTS: {success_count} OK  |  {skip_count} skipped  |  Total time: {total_time/60:.1f} min")
    print(f"Output saved to: {OUT_DIR}")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
