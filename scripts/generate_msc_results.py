"""
Generate Stage 2 & 3 figures (13-21) for SPION-labelled Mesenchymal Stem Cells (MSCs).

Cell parameters derived from literature:
  - Radius:          12.5 µm  (diameter ~25 µm; conservative culture estimate)
  - SPION loading:   25 pg    (centre of 10-30 pg range typical for standard labeling;
                               up to ~30 pg with poly-L-lysine or transfection agents)
  - Cell density:    1050 kg/m³ (not currently modelled — gravity/buoyancy excluded
                                 per low-Re assumption; documented for reference)
  - Susceptibility:  2.0 (SI)   (magnetite, same as endothelial cell variant)
  - SPION density:   5170 kg/m³  (magnetite Fe3O4, same as endothelial cell variant)

Key references for SPION labelling of MSCs:
  1. Arbab et al. (2003) Radiology — 30.1 ± 3.7 pg Fe/cell with ferumoxides + transfection agent
  2. Farrell et al. (2008) — 22 pg Fe/cell with poly-L-lysine, no cytotoxicity
  3. Saldanha et al. (2021) PubMed 33747602 — 5-10 pg Fe/cell passive incubation with CMF particles

Physics note:
  - Larger MSC radius (12.5 vs 10 µm) increases Stokes drag by 25% at equal velocity
  - Lower SPION loading (25 vs 200 pg) reduces magnetic force by ~8x
  - Net effect: MSCs are significantly harder to capture than 200 pg endothelial cells,
    but the trajectory advantage (drift over 2 mm approach) remains meaningful

All other parameters (stent geometry, B0, blood flow) identical to original project defaults.
"""

import sys
import time
from pathlib import Path

proj_root = Path(__file__).parent.parent
sys.path.insert(0, str(proj_root))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

from stent_capture.figures.common import DEFAULTS, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force
from stent_capture.physics.hydrodynamics import BloodFlow, stokes_drag
from stent_capture.physics.capture_criterion import capture_distance
from stent_capture.core.gradient import compute_gradient_magnitude
from stent_capture.simulation.trajectories import integrate_trajectory

# ── Output directory ──────────────────────────────────────────────────────────
OUT_DIR = proj_root / "results-msc"
OUT_DIR.mkdir(exist_ok=True)

# ── MSC cell parameters ───────────────────────────────────────────────────────
_CELL_RADIUS_M   = 12.5e-6   # m  — 12.5 µm radius (25 µm diameter)
_SPION_MASS_PG   = 25.0      # pg — default loading (centre of 10-30 pg range)
_CELL_DENSITY    = 1050      # kg/m³ — for documentation only (not used in physics)

# ── Shared physics parameters (identical to all other variants) ───────────────
_B0_Z   = 0.5       # T
_R_VES  = 1.54e-3   # m
_V_MCA  = 0.2       # m/s

# ── Sweep ranges appropriate for MSCs ─────────────────────────────────────────
_LOADINGS_PG  = np.array([5., 10., 15., 25., 50., 75., 100.])   # pg
_VELOCITIES   = np.array([0.02, 0.05, 0.10, 0.20, 0.50])        # m/s
_N_ITER       = 10  # binary search depth (~1.4 µm resolution; 7 gave ~11 µm, causing fig19/fig21 discrepancy)

# ── Trajectory integration settings ──────────────────────────────────────────
_TRAJ_KW = dict(z_end=2e-3, max_time=1.5, rtol=1e-5, atol=1e-8)

# ── Colours ───────────────────────────────────────────────────────────────────
_COLOR_STATIC  = "#2980b9"
_COLOR_TRAJ    = "#c0392b"
_COLOR_MSC     = "#27ae60"


def _make_msc_cell(spion_mass_pg=None):
    """Return SPIONLabelledCell configured for an MSC."""
    if spion_mass_pg is None:
        spion_mass_pg = _SPION_MASS_PG
    return SPIONLabelledCell(
        radius=_CELL_RADIUS_M,
        spion_mass_per_cell=spion_mass_pg * 1e-15,
    )


def _find_trajectory_capture_range(cell, tf, flow, ring, n_iter=None):
    """
    Binary search for the maximum injection distance (from stent inner surface)
    at which a cell is still captured.

    Returns distance in metres.
    """
    if n_iter is None:
        n_iter = _N_ITER

    r_inner = ring.R - ring.t / 2
    lo = 0.0           # always captured (right at wall)
    hi = r_inner * 0.95  # upper bound (near vessel axis)

    # Quick check: is any capture possible?
    traj_lo = integrate_trajectory(
        cell, tf, flow, ring,
        np.array([r_inner - 5e-6, 0.0, -2e-3]),
        **_TRAJ_KW,
    )
    if traj_lo.status != "captured":
        return 0.0  # no capture even from wall proximity

    for _ in range(n_iter):
        mid = (lo + hi) / 2.0
        r_inject = r_inner - mid
        traj = integrate_trajectory(
            cell, tf, flow, ring,
            np.array([r_inject, 0.0, -2e-3]),
            **_TRAJ_KW,
        )
        if traj.status == "captured":
            lo = mid   # can capture from at least this far — try farther
        else:
            hi = mid   # too far — reduce

    return lo  # maximum capture distance


# ─────────────────────────────────────────────────────────────────────────────
# Fig 13: Force Parameter (B-field only, same as other variants)
# ─────────────────────────────────────────────────────────────────────────────

def fig13_msc():
    """Force parameter = |B| * |grad B| — cell-independent, included for completeness."""
    print("  Generating fig13 (force parameter)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    r_pts = np.linspace(1.3e-3, 1.6e-3, 120)
    z_pts = np.linspace(-2e-3, 2e-3, 200)
    R, Z = np.meshgrid(r_pts, z_pts)
    pts = np.column_stack([R.ravel(), np.zeros(R.size), Z.ravel()])

    B_mag  = np.linalg.norm(tf.field_at(pts), axis=1)
    grad_B = compute_gradient_magnitude(tf.field_at, pts, dx=5e-7)
    FP     = (B_mag * grad_B).reshape(R.shape)

    fig, ax = plt.subplots(figsize=(10, 7))
    levels = np.logspace(np.log10(FP[FP > 0].min()), np.log10(FP.max()), 20)
    im = ax.contourf(R * 1e3, Z * 1e3, FP, levels=levels, cmap="YlOrRd", extend="both")
    fig.colorbar(im, ax=ax, label="Force parameter |B|·|∇|B|| (T²/m)")
    ax.set_xlabel("Radial distance (mm)", fontsize=11)
    ax.set_ylabel("Axial distance (mm)", fontsize=11)
    ax.set_title(
        f"Fig 13 (MSC) — Force parameter = |B| · |∇|B|| (cell-type-independent)\n"
        f"B₀ = {_B0_Z} T (axial), strut-aligned axis (θ=0) — identical to original fig13", fontsize=12
    )

    fig.savefig(OUT_DIR / "fig13_force_parameter_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig13_force_parameter_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved fig13", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 14: Magnetic Force vs Radial Distance
# ─────────────────────────────────────────────────────────────────────────────

def fig14_msc():
    """Magnetic force on MSC vs radial distance — compare to 10 pg endothelial cell."""
    print("  Generating fig14 (force vs distance)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    r_vals = np.linspace(1.46e-3, 1.65e-3, 100)
    pts = np.column_stack([r_vals, np.zeros_like(r_vals), np.zeros_like(r_vals)])

    cell_msc     = _make_msc_cell()
    cell_ref_10  = SPIONLabelledCell(spion_mass_per_cell=10e-15)    # endothelial lower regime
    cell_ref_200 = SPIONLabelledCell(spion_mass_per_cell=200e-15)   # Polyak 2008 working dose

    F_msc     = np.linalg.norm(magnetic_force(cell_msc,     tf, pts), axis=1) * 1e12
    F_ref_10  = np.linalg.norm(magnetic_force(cell_ref_10,  tf, pts), axis=1) * 1e12
    F_ref_200 = np.linalg.norm(magnetic_force(cell_ref_200, tf, pts), axis=1) * 1e12

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(r_vals * 1e3, F_msc,     'o-',  color=_COLOR_MSC,   lw=2.5, ms=5,
                label=f"MSC ({_SPION_MASS_PG:.0f} pg, r={_CELL_RADIUS_M*1e6:.1f} µm)")
    ax.semilogy(r_vals * 1e3, F_ref_10,  's--', color="gray",       lw=1.5, ms=4, alpha=0.7,
                label="Endothelial (10 pg, r=10 µm) — lower regime")
    ax.semilogy(r_vals * 1e3, F_ref_200, '^--', color="steelblue",  lw=1.5, ms=4, alpha=0.7,
                label="Endothelial (200 pg, r=10 µm) — Polyak 2008")

    ax.set_xlabel("Radial distance from axis (mm)", fontsize=11)
    ax.set_ylabel("Magnetic force (pN)", fontsize=11)
    ax.set_title(
        f"Fig 14 (MSC) — Magnetic force on SPION-labelled cell\n"
        f"v = 0 (stationary), B₀ = {_B0_Z} T, θ = 0, z = 0", fontsize=12
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=10)

    fig.savefig(OUT_DIR / "fig14_force_vs_distance_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig14_force_vs_distance_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved fig14", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 15: Stokes Drag vs Flow Velocity
# ─────────────────────────────────────────────────────────────────────────────

def fig15_msc():
    """Stokes drag on MSC vs flow velocity — larger radius increases drag by 25% vs 10 µm cell."""
    print("  Generating fig15 (drag vs velocity)...", flush=True)

    cell_msc    = _make_msc_cell()
    cell_ref_10 = SPIONLabelledCell(spion_mass_per_cell=10e-15)

    v_vals = np.linspace(0.01, 1.0, 100)
    r_probe = np.array([[1.5e-3, 0.0, 0.0]])

    F_msc_vals, F_ref_vals = [], []
    for v in v_vals:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        F_msc_vals.append(np.linalg.norm(stokes_drag(cell_msc,    flow, r_probe)) * 1e12)
        F_ref_vals.append(np.linalg.norm(stokes_drag(cell_ref_10, flow, r_probe)) * 1e12)

    F_msc_vals = np.array(F_msc_vals)
    F_ref_vals = np.array(F_ref_vals)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(v_vals, F_msc_vals, 'o-', color=_COLOR_MSC,  lw=2.5, ms=4,
                label=f"MSC (r = {_CELL_RADIUS_M*1e6:.1f} µm)")
    ax.semilogy(v_vals, F_ref_vals, 's--', color="gray", lw=1.5, ms=3, alpha=0.7,
                label="Endothelial (r = 10 µm)")

    ratio = F_msc_vals[0] / F_ref_vals[0]
    ax.text(0.05, 0.92,
            f"Drag ratio MSC/EPC = {ratio:.2f}× (∝ radius ratio {_CELL_RADIUS_M*1e6/10:.2f}×)",
            transform=ax.transAxes, fontsize=9, color="navy",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    ax.set_xlabel("Mean blood velocity (m/s)", fontsize=11)
    ax.set_ylabel("Stokes drag force (pN)", fontsize=11)
    ax.set_title(
        f"Fig 15 (MSC) — Stokes drag on SPION-labelled cell\n"
        f"at r = 1.50 mm (mid-lumen), η = 4 mPa·s", fontsize=12
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=10)

    fig.savefig(OUT_DIR / "fig15_drag_vs_velocity_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig15_drag_vs_velocity_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved fig15", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 16: Static Capture Map (3-panel velocity sweep)
# ─────────────────────────────────────────────────────────────────────────────

_V_CASES_16  = [0.05, 0.2, 0.5]   # m/s  (same panels as original fig16)
_GRID_N      = 55
_LOG_MIN, _LOG_MAX = -4.0, 1.0


def _build_capture_grid():
    lim  = _R_VES * 1.05
    x_v  = np.linspace(-lim, lim, _GRID_N)
    xx, yy = np.meshgrid(x_v, x_v)
    pts  = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(_GRID_N ** 2)])
    r_grid = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    r_lumen = DEFAULTS["R"] - DEFAULTS["t"] / 2
    inside  = r_grid < r_lumen
    return xx, yy, pts, inside


def _strut_markers(ax):
    R = DEFAULTS["R"];  t = DEFAULTS["t"];  w = DEFAULTS["w"];  n = DEFAULTS["n_struts"]
    mr = max(t, w) / 2 * 1e3
    for k in range(n):
        th = 2 * np.pi * k / n
        ax.add_patch(plt.Circle((R * np.cos(th) * 1e3, R * np.sin(th) * 1e3),
                                mr, color="white", zorder=4, linewidth=0))


def fig16_msc():
    """Static capture criterion map for MSC — 3-panel velocity sweep."""
    print("  Generating fig16 (capture map)...", flush=True)

    xx, yy, pts, inside = _build_capture_grid()

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    cell = _make_msc_cell()
    F_mag_all = np.linalg.norm(magnetic_force(cell, tf, pts), axis=1) * 1e12

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    max_ratio = 0.0

    for ax, v_mean in zip(axes, _V_CASES_16):
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v_mean)
        F_drag_all = np.linalg.norm(
            np.column_stack([stokes_drag(cell, flow, pts[:, :].copy())]), axis=1
        ) * 1e12

        # recompute drag properly row-by-row
        F_drag_vec = stokes_drag(cell, flow, pts)
        F_drag_all = np.linalg.norm(F_drag_vec, axis=1) * 1e12

        ratio = np.full(len(pts), np.nan)
        safe  = np.where(F_drag_all > 0, F_drag_all, np.nan)
        ratio[inside] = (F_mag_all / safe)[inside]

        valid = ratio[np.isfinite(ratio)]
        if len(valid):
            max_ratio = max(max_ratio, float(valid.max()))

        log_r = np.full_like(ratio, np.nan)
        pos   = inside & (ratio > 0)
        log_r[pos] = np.log10(ratio[pos])
        log_grid = log_r.reshape(xx.shape)
        clipped  = np.clip(log_grid, _LOG_MIN, _LOG_MAX)

        norm = TwoSlopeNorm(vmin=_LOG_MIN, vcenter=0.0, vmax=_LOG_MAX)
        im = ax.pcolormesh(xx * 1e3, yy * 1e3, clipped,
                           cmap=plt.get_cmap("RdYlGn"), norm=norm,
                           shading="auto", zorder=1)

        vmask = np.isfinite(log_grid)
        if np.any(log_grid[vmask] > 0) and np.any(log_grid[vmask] < 0):
            ax.contour(xx * 1e3, yy * 1e3, log_grid, levels=[0],
                       colors="black", linewidths=1.5, zorder=2)

        theta = np.linspace(0, 2 * np.pi, 300)
        R_lum = DEFAULTS["R"] - DEFAULTS["t"] / 2
        R_out = DEFAULTS["R"] + DEFAULTS["t"] / 2
        ax.plot(R_lum * 1e3 * np.cos(theta), R_lum * 1e3 * np.sin(theta),
                "w--", lw=1.5, alpha=0.9, zorder=3, label="Lumen boundary")
        ax.plot(R_out * 1e3 * np.cos(theta), R_out * 1e3 * np.sin(theta),
                "w-",  lw=0.5, alpha=0.5, zorder=2)
        _strut_markers(ax)

        ax.set_aspect("equal")
        ax.set_xlim(-_R_VES * 1.05e3, _R_VES * 1.05e3)
        ax.set_ylim(-_R_VES * 1.05e3, _R_VES * 1.05e3)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)" if ax is axes[0] else "")
        ax.set_title(f"v_mean = {v_mean:.2f} m/s")

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if ax is axes[-1]:
            cbar.set_label("log10(|F_mag| / |F_drag|)")
        cbar.set_ticks([_LOG_MIN, -3, -2, -1, 0, _LOG_MAX])
        cbar.set_ticklabels(
            [f"10^{int(t)}" if t != 0 else "1 (balance)"
             for t in [_LOG_MIN, -3, -2, -1, 0, _LOG_MAX]]
        )
        cbar.ax.axhline(0, color="black", lw=1.5)

    fig.suptitle(
        f"Fig 16 (MSC) — Force ratio |F_mag|/|F_drag| map (z=0 cross-section)\n"
        f"Green: |F_mag| > |F_drag| (capture possible) | "
        f"Red: drag dominates | B₀ = {_B0_Z} T axial\n"
        f"MSC: r = {_CELL_RADIUS_M*1e6:.1f} µm, {_SPION_MASS_PG:.0f} pg SPION loading. "
        f"Max ratio = {max_ratio:.2f}",
        fontsize=11, y=0.98
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    fig.savefig(OUT_DIR / "fig16_capture_map_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig16_capture_map_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig16 (max ratio = {max_ratio:.3f})", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 17: SPION Loading Sweep
# ─────────────────────────────────────────────────────────────────────────────

def fig17_msc():
    """Static capture distance vs SPION loading — MSC radius throughout."""
    print("  Generating fig17 (loading sweep)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)

    static_dists = []
    for m_pg in _LOADINGS_PG:
        cell = _make_msc_cell(m_pg)
        d_s  = capture_distance(cell, tf, flow, direction="inward") * 1e6
        static_dists.append(d_s)

    static_dists = np.array(static_dists)
    positive_dists = static_dists[static_dists > 0]
    y_zero_label = positive_dists.min() * 0.4 if len(positive_dists) else 0.1

    fig, ax = plt.subplots(figsize=(10, 6))
    plot_vals = np.where(static_dists > 0, static_dists, y_zero_label * 0.5)
    ax.loglog(_LOADINGS_PG, plot_vals, 'o-', color=_COLOR_MSC,
              lw=2.5, ms=8, markerfacecolor="white", markeredgewidth=2,
              label="MSC (r = 12.5 µm)")

    for m, d in zip(_LOADINGS_PG, static_dists):
        if d > 0:
            ax.annotate(f"{d:.1f} µm", xy=(m, d), xytext=(0, 8),
                       textcoords="offset points", ha="center", fontsize=9)
        else:
            ax.annotate("0 µm\n(no static\ncapture)", xy=(m, y_zero_label * 0.5),
                       ha="center", fontsize=7, color="red")

    ax.axvline(_SPION_MASS_PG, color="green", ls="--", lw=2, alpha=0.8,
               label=f"{_SPION_MASS_PG:.0f} pg (MSC default)")
    ax.set_xlabel("SPION loading per cell (pg)", fontsize=11)
    ax.set_ylabel("Static capture distance (µm)", fontsize=11)
    ax.set_title(
        f"Fig 17 (MSC) — Static capture distance vs SPION loading\n"
        f"v = {_V_MCA} m/s, B₀ = {_B0_Z} T, MSC radius = {_CELL_RADIUS_M*1e6:.1f} µm",
        fontsize=12
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=10)

    fig.savefig(OUT_DIR / "fig17_spion_loading_sweep_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig17_spion_loading_sweep_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved fig17", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 18: Single Trajectory
# ─────────────────────────────────────────────────────────────────────────────

def fig18_msc():
    """Single-cell trajectory for MSC at 50 µm injection distance."""
    print("  Generating fig18 (single trajectory)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
    cell = _make_msc_cell()

    r_inner  = ring.R - ring.t / 2
    r_inject = r_inner - 50e-6   # 50 µm from stent inner surface

    traj = integrate_trajectory(
        cell, tf, flow, ring,
        np.array([r_inject, 0.0, -2e-3]),
        **_TRAJ_KW,
    )

    r_traj = np.sqrt(traj.positions[:, 0] ** 2 + traj.positions[:, 1] ** 2)
    z_traj = traj.positions[:, 2]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.plot(r_traj * 1e3, z_traj * 1e3, 'o-', color=_COLOR_TRAJ, ms=3, lw=1.5, label="Trajectory")
    ax1.axhline(0, color="gray", ls="-", lw=2, alpha=0.5, label="Stent plane")
    ax1.axvline(r_inner * 1e3, color="green", ls="--", lw=2, alpha=0.7, label="Stent inner surface")
    ax1.set_xlabel("Radial distance (mm)", fontsize=11)
    ax1.set_ylabel("Axial distance (mm)", fontsize=11)
    ax1.set_title("R-Z projection", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    ax2.plot(traj.times * 1e3, r_traj * 1e3, 'o-', color=_COLOR_TRAJ, ms=3, lw=1.5, label="Radial position")
    ax2.axhline(r_inner * 1e3, color="green", ls="--", lw=2, alpha=0.7, label="Capture threshold")
    ax2.set_xlabel("Time (ms)", fontsize=11)
    ax2.set_ylabel("Radial distance (mm)", fontsize=11)
    ax2.set_title("Radial approach over time", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        f"Fig 18 (MSC) — Single-cell trajectory (injection at r = 50 µm, status = {traj.status})\n"
        f"MSC: r = {_CELL_RADIUS_M*1e6:.1f} µm, {_SPION_MASS_PG:.0f} pg SPION, "
        f"v = {_V_MCA} m/s, B₀ = {_B0_Z} T",
        fontsize=12, y=1.00
    )
    fig.tight_layout()

    fig.savefig(OUT_DIR / "fig18_single_trajectory_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig18_single_trajectory_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved fig18 (status: {traj.status})", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 19: Trajectory Bundle
# ─────────────────────────────────────────────────────────────────────────────

def fig19_msc():
    """Trajectory bundle for MSC — multiple injection distances."""
    print("  Generating fig19 (trajectory bundle)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
    cell = _make_msc_cell()

    r_inner  = ring.R - ring.t / 2
    distances_um = [10, 20, 40, 60, 80, 100]   # µm from stent inner surface
    colors = plt.cm.RdYlGn(np.linspace(0.15, 0.85, len(distances_um)))

    fig, ax = plt.subplots(figsize=(11, 7))

    for d_um, color in zip(distances_um, colors):
        r_inj = r_inner - d_um * 1e-6
        traj  = integrate_trajectory(
            cell, tf, flow, ring,
            np.array([r_inj, 0.0, -2e-3]),
            **_TRAJ_KW,
        )
        r_traj = np.sqrt(traj.positions[:, 0] ** 2 + traj.positions[:, 1] ** 2)
        z_traj = traj.positions[:, 2]

        sym = "OK" if traj.status == "captured" else "X"
        ax.plot(z_traj * 1e3, r_traj * 1e3, 'o-', color=color, lw=2, ms=3,
               label=f"{d_um:3d} µm ({sym})", alpha=0.85)

    ax.axhline(r_inner * 1e3, color="black", ls="--", lw=2.5, alpha=0.8, label="Stent inner surface")
    ax.axvline(0, color="gray", ls="-", lw=1.5, alpha=0.5)
    ax.set_xlabel("Axial distance (mm)", fontsize=11)
    ax.set_ylabel("Radial distance (mm)", fontsize=11)
    ax.set_title(
        f"Fig 19 (MSC) — Trajectory bundle: injection distance sweep\n"
        f"MSC: r = {_CELL_RADIUS_M*1e6:.1f} µm, {_SPION_MASS_PG:.0f} pg SPION, "
        f"v = {_V_MCA} m/s, B₀ = {_B0_Z} T",
        fontsize=12
    )
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-2.5, 2.5)

    fig.savefig(OUT_DIR / "fig19_trajectory_bundle_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig19_trajectory_bundle_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved fig19", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 20: Capture Efficiency
# ─────────────────────────────────────────────────────────────────────────────

def fig20_msc():
    """Static capture distance vs loading and velocity for MSC."""
    print("  Generating fig20 (capture efficiency)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: vs loading at v = 0.2 m/s
    flow_ref = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
    d_static_load = []
    for m_pg in _LOADINGS_PG:
        cell = _make_msc_cell(m_pg)
        d_static_load.append(capture_distance(cell, tf, flow_ref, direction="inward") * 1e6)

    ax_a.loglog(_LOADINGS_PG, np.maximum(d_static_load, 0.01), 'o-', color=_COLOR_MSC,
               lw=2.5, ms=7, markerfacecolor="white", markeredgewidth=1.8,
               label="MSC (r = 12.5 µm)")
    for m, d in zip(_LOADINGS_PG, d_static_load):
        if d <= 0:
            ax_a.annotate("0 µm\n(no static\ncapture)", xy=(m, 0.01),
                         xytext=(0, 6), textcoords="offset points",
                         ha="center", fontsize=7, color="red")
    ax_a.axvline(_SPION_MASS_PG, color="green", ls="--", alpha=0.8, lw=1.5,
                label=f"{_SPION_MASS_PG:.0f} pg default")
    ax_a.set_xlabel("SPION loading (pg)", fontsize=11)
    ax_a.set_ylabel("Static capture distance (µm)", fontsize=11)
    ax_a.set_title(f"(a) vs SPION loading\nv = {_V_MCA} m/s", fontsize=11)
    ax_a.grid(True, which="both", alpha=0.3)
    ax_a.legend(fontsize=10)

    # Panel B: vs velocity at default SPION loading
    cell_def = _make_msc_cell()
    d_static_vel = []
    for v in _VELOCITIES:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        d_static_vel.append(capture_distance(cell_def, tf, flow, direction="inward") * 1e6)

    ax_b.loglog(_VELOCITIES, np.maximum(d_static_vel, 0.01), 's-', color="#9b59b6",
               lw=2.5, ms=7, markerfacecolor="white", markeredgewidth=1.8,
               label=f"MSC {_SPION_MASS_PG:.0f} pg")
    for v, d in zip(_VELOCITIES, d_static_vel):
        if d <= 0:
            ax_b.annotate("0 µm", xy=(v, 0.01), xytext=(0, 6),
                         textcoords="offset points",
                         ha="center", fontsize=7, color="red")
    ax_b.axvline(_V_MCA, color="blue", ls="--", alpha=0.7, lw=1.5, label="0.2 m/s (MCA)")
    ax_b.set_xlabel("Mean velocity (m/s)", fontsize=11)
    ax_b.set_ylabel("Static capture distance (µm)", fontsize=11)
    ax_b.set_title(f"(b) vs velocity\n{_SPION_MASS_PG:.0f} pg MSC", fontsize=11)
    ax_b.grid(True, which="both", alpha=0.3)
    ax_b.legend(fontsize=10)

    fig.suptitle(
        f"Fig 20 (MSC) — Static capture efficiency trends\n"
        f"B₀ = {_B0_Z} T, r_cell = {_CELL_RADIUS_M*1e6:.1f} µm, θ = 0, z = 0",
        fontsize=12, y=1.00
    )
    fig.tight_layout()

    fig.savefig(OUT_DIR / "fig20_capture_efficiency_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig20_capture_efficiency_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved fig20", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Fig 21: Static vs Trajectory — headline comparison
# ─────────────────────────────────────────────────────────────────────────────

def fig21_msc():
    """
    Headline figure: static vs trajectory capture range for MSC.

    Panel (a): loading sweep at v = 0.2 m/s
    Panel (b): velocity sweep at 25 pg (MSC default)
    """
    print("  Generating fig21 (static vs trajectory — this takes a few minutes)...", flush=True)

    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    # ── Panel (a): loading sweep ──────────────────────────────────────────────
    print("    Panel (a): loading sweep...", flush=True)
    static_a, traj_a = [], []
    for m_pg in _LOADINGS_PG:
        cell = _make_msc_cell(m_pg)
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MCA)
        d_s  = capture_distance(cell, tf, flow, direction="inward") * 1e6
        d_t  = _find_trajectory_capture_range(cell, tf, flow, ring) * 1e6
        static_a.append(d_s)
        traj_a.append(d_t)
        print(f"      {m_pg:5.0f} pg: static={d_s:.1f} µm, traj={d_t:.1f} µm", flush=True)

    # ── Panel (b): velocity sweep ─────────────────────────────────────────────
    print("    Panel (b): velocity sweep...", flush=True)
    cell_def = _make_msc_cell()
    static_b, traj_b = [], []
    for v in _VELOCITIES:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        d_s  = capture_distance(cell_def, tf, flow, direction="inward") * 1e6
        d_t  = _find_trajectory_capture_range(cell_def, tf, flow, ring) * 1e6
        static_b.append(d_s)
        traj_b.append(d_t)
        print(f"      v={v:.3f} m/s: static={d_s:.1f} µm, traj={d_t:.1f} µm", flush=True)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a)
    ax_a.plot(_LOADINGS_PG, static_a, 'o-', color=_COLOR_STATIC, lw=2.5, ms=8,
             markerfacecolor="white", markeredgewidth=2, label="Static criterion")
    ax_a.plot(_LOADINGS_PG, traj_a,   's-', color=_COLOR_TRAJ,   lw=2.5, ms=8,
             markerfacecolor="white", markeredgewidth=2, label="Trajectory (RK45)")

    for m, s, t in zip(_LOADINGS_PG, static_a, traj_a):
        if t > 1:
            ax_a.annotate(f"{t:.0f}", xy=(m, t), xytext=(0, 8),
                         textcoords="offset points", ha="center", fontsize=8, color=_COLOR_TRAJ)
        if s > 1:
            ax_a.annotate(f"{s:.0f}", xy=(m, s), xytext=(0, -14),
                         textcoords="offset points", ha="center", fontsize=8, color=_COLOR_STATIC)

    ax_a.axvline(_SPION_MASS_PG, color="green", ls="--", lw=1.5, alpha=0.7,
                label=f"{_SPION_MASS_PG:.0f} pg default")
    ax_a.set_xlabel("SPION loading per MSC (pg)", fontsize=11)
    ax_a.set_ylabel("Capture range from stent surface (µm)", fontsize=11)
    ax_a.set_title(f"(a) Loading sweep at v = {_V_MCA} m/s\nMSC radius = {_CELL_RADIUS_M*1e6:.1f} µm", fontsize=11)
    ax_a.legend(fontsize=10)
    ax_a.grid(True, alpha=0.3)
    ax_a.set_xlim(left=0)
    ax_a.set_ylim(bottom=0)

    # Panel (b)
    ax_b.plot(_VELOCITIES, static_b, 'o-', color=_COLOR_STATIC, lw=2.5, ms=8,
             markerfacecolor="white", markeredgewidth=2, label="Static criterion")
    ax_b.plot(_VELOCITIES, traj_b,   's-', color=_COLOR_TRAJ,   lw=2.5, ms=8,
             markerfacecolor="white", markeredgewidth=2, label="Trajectory (RK45)")

    for v, s, t in zip(_VELOCITIES, static_b, traj_b):
        if t > 1:
            ax_b.annotate(f"{t:.0f}", xy=(v, t), xytext=(0, 8),
                         textcoords="offset points", ha="center", fontsize=8, color=_COLOR_TRAJ)
        if s > 1:
            ax_b.annotate(f"{s:.0f}", xy=(v, s), xytext=(0, -14),
                         textcoords="offset points", ha="center", fontsize=8, color=_COLOR_STATIC)

    ax_b.axvline(_V_MCA, color="blue", ls="--", lw=1.5, alpha=0.7, label="0.2 m/s (MCA)")
    ax_b.set_xlabel("Mean blood velocity (m/s)", fontsize=11)
    ax_b.set_ylabel("Capture range from stent surface (µm)", fontsize=11)
    ax_b.set_title(f"(b) Velocity sweep at {_SPION_MASS_PG:.0f} pg\nMSC radius = {_CELL_RADIUS_M*1e6:.1f} µm", fontsize=11)
    ax_b.legend(fontsize=10)
    ax_b.grid(True, alpha=0.3)
    ax_b.set_xlim(left=0)
    ax_b.set_ylim(bottom=0)

    # Reference line at MCA, default loading
    ref_static = static_b[_VELOCITIES.tolist().index(0.2)] if 0.2 in _VELOCITIES else None
    ref_traj   = traj_b[  _VELOCITIES.tolist().index(0.2)] if 0.2 in _VELOCITIES else None
    ratio_str  = f"{ref_traj/ref_static:.1f}×" if ref_static and ref_static > 0 else "∞ (static = 0)"

    fig.suptitle(
        f"Fig 21 (MSC) — Static vs trajectory capture range\n"
        f"MSC: r = {_CELL_RADIUS_M*1e6:.1f} µm, B₀ = {_B0_Z} T | "
        f"At {_SPION_MASS_PG:.0f} pg, v = {_V_MCA} m/s: trajectory extends static by {ratio_str}",
        fontsize=12, y=1.02
    )
    fig.tight_layout()

    fig.savefig(OUT_DIR / "fig21_static_vs_trajectory_msc.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig21_static_vs_trajectory_msc.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Saved fig21", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print(f"\n{'='*80}")
    print("GENERATING MSC RESULTS (MESENCHYMAL STEM CELLS)")
    print(f"  Cell radius : {_CELL_RADIUS_M*1e6:.1f} um")
    print(f"  SPION load  : {_SPION_MASS_PG:.0f} pg (default); sweep {_LOADINGS_PG[0]:.0f}-{_LOADINGS_PG[-1]:.0f} pg")
    print(f"  Cell density: {_CELL_DENSITY} kg/m3 (documented; not modelled in physics)")
    print(f"Output directory: {OUT_DIR}")
    print(f"{'='*80}\n")

    t_start = time.time()
    success_count = 0
    fail_count    = 0

    figures = [
        ("fig13", fig13_msc, "Force parameter"),
        ("fig14", fig14_msc, "Force vs distance"),
        ("fig15", fig15_msc, "Drag vs velocity"),
        ("fig16", fig16_msc, "Capture map"),
        ("fig17", fig17_msc, "SPION loading sweep"),
        ("fig18", fig18_msc, "Single trajectory"),
        ("fig19", fig19_msc, "Trajectory bundle"),
        ("fig20", fig20_msc, "Capture efficiency"),
        ("fig21", fig21_msc, "Static vs trajectory"),
    ]

    for fig_num, fig_func, description in figures:
        print(f"\n{fig_num:6s} | {description}", flush=True)
        try:
            t0 = time.time()
            fig_func()
            elapsed = time.time() - t0
            print(f"       OK ({elapsed:.1f}s)", flush=True)
            success_count += 1
        except Exception as e:
            import traceback
            print(f"       FAILED: {e}")
            traceback.print_exc()
            fail_count += 1

    total_time = time.time() - t_start
    print(f"\n{'='*80}")
    print(f"RESULTS: {success_count} OK  |  {fail_count} FAILED  |  Total: {total_time/60:.1f} min")
    print(f"Output saved to: {OUT_DIR}")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
