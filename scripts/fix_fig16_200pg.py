"""
Regenerate fig16 for 200 pg using the SAME capture_map implementation as the original.
"""

import sys
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

OUT_DIR = proj_root / "results-200pg"

# Parameters - same as original
_V_CASES = [0.05, 0.2, 0.5]   # m/s
_R_VES   = 1.54e-3             # m
_B0_Z    = 0.5                 # T
_CELL    = SPIONLabelledCell(spion_mass_per_cell=200e-15)  # 200 PG - THE KEY CHANGE
_GRID_N  = 55                  # grid points per axis

_LOG_MIN = -4.0
_LOG_MAX =  1.0


def _build_grid():
    """Return (xx, yy, pts, inside_lumen) for the 2D cross-section grid."""
    lim = _R_VES * 1.05
    x_vals = np.linspace(-lim, lim, _GRID_N)
    y_vals = np.linspace(-lim, lim, _GRID_N)
    xx, yy = np.meshgrid(x_vals, y_vals)
    pts = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(_GRID_N ** 2)])

    r_grid = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_lumen = R - t / 2
    inside_lumen = r_grid < r_lumen
    return xx, yy, pts, inside_lumen


def _compute_F_mag(pts: np.ndarray) -> np.ndarray:
    ring = make_ring()
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    F = magnetic_force(_CELL, tf, pts)
    return np.linalg.norm(F, axis=1) * 1e12   # pN


def _compute_F_drag(pts: np.ndarray, v_mean: float) -> np.ndarray:
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v_mean)
    F = stokes_drag(_CELL, flow, pts)
    return np.linalg.norm(F, axis=1) * 1e12   # pN


def _add_strut_markers(ax):
    """Draw filled white circles at each strut centre position (in mm)."""
    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    w = DEFAULTS["w"]
    n = DEFAULTS["n_struts"]
    marker_radius = max(t, w) / 2 * 1e3   # mm
    for k in range(n):
        th = 2 * np.pi * k / n
        cx = R * np.cos(th) * 1e3   # mm
        cy = R * np.sin(th) * 1e3   # mm
        circle = plt.Circle((cx, cy), marker_radius,
                             color="white", zorder=4, linewidth=0)
        ax.add_patch(circle)


def make_figure_200pg():
    """Generate fig16 with 200 pg SPION loading."""
    print("  Building grid and computing magnetic force (once)...")
    xx, yy, pts, inside_lumen = _build_grid()
    F_mag_all = _compute_F_mag(pts)

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    max_ratio_overall = 0.0

    for ax, v_mean in zip(axes, _V_CASES):
        print(f"    Computing drag for v_mean = {v_mean} m/s...")
        F_drag_all = _compute_F_drag(pts, v_mean)

        # Compute force ratio; NaN outside lumen
        ratio = np.full(len(pts), np.nan)
        drag_safe = np.where(F_drag_all > 0, F_drag_all, np.nan)
        ratio[inside_lumen] = (F_mag_all / drag_safe)[inside_lumen]

        # Track maximum ratio
        valid_ratio = ratio[np.isfinite(ratio)]
        if len(valid_ratio):
            max_ratio_overall = max(max_ratio_overall, float(valid_ratio.max()))

        # Log10 of ratio
        log_ratio = np.full_like(ratio, np.nan)
        pos = inside_lumen & (ratio > 0)
        log_ratio[pos] = np.log10(ratio[pos])
        log_ratio_grid = log_ratio.reshape(xx.shape)

        # Clip to display range
        log_ratio_clipped = np.clip(log_ratio_grid, _LOG_MIN, _LOG_MAX)

        # Colormap
        norm = TwoSlopeNorm(vmin=_LOG_MIN, vcenter=0.0, vmax=_LOG_MAX)
        cmap = plt.get_cmap("RdYlGn")
        im = ax.pcolormesh(xx * 1e3, yy * 1e3, log_ratio_clipped,
                           cmap=cmap, norm=norm, shading="auto", zorder=1)

        # Capture boundary contour at ratio = 1
        valid_mask = np.isfinite(log_ratio_grid)
        if (np.any(log_ratio_grid[valid_mask] > 0) and
                np.any(log_ratio_grid[valid_mask] < 0)):
            ax.contour(xx * 1e3, yy * 1e3, log_ratio_grid, levels=[0],
                       colors="black", linewidths=1.5, zorder=2)

        # Lumen inner boundary
        theta = np.linspace(0, 2 * np.pi, 300)
        R_lumen = DEFAULTS["R"] - DEFAULTS["t"] / 2
        ax.plot(R_lumen * 1e3 * np.cos(theta),
                R_lumen * 1e3 * np.sin(theta),
                "w--", lw=1.5, alpha=0.9, zorder=3, label="Lumen boundary")

        # Stent outer circle
        R_outer = DEFAULTS["R"] + DEFAULTS["t"] / 2
        ax.plot(R_outer * 1e3 * np.cos(theta),
                R_outer * 1e3 * np.sin(theta),
                "w-", lw=0.5, alpha=0.5, zorder=2)

        # Strut markers
        _add_strut_markers(ax)

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

    print(f"    Maximum force ratio in lumen: {max_ratio_overall:.4f}")

    fig.suptitle(
        f"Fig 16 (200 pg) - Force ratio |F_mag|/|F_drag| map (z=0 cross-section)\n"
        f"Green: |F_mag| > |F_drag| (capture possible) | "
        f"Red: drag dominates | B0 = {_B0_Z} T axial\n"
        f"200 pg SPION loading per cell (vs 10 pg default). "
        f"Max ratio = {max_ratio_overall:.2f}",
        fontsize=11, y=0.98
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    return fig


def main():
    print("\nRegenerating fig16 for 200 pg (matching original design)...")
    fig = make_figure_200pg()

    # Remove old fig16_capture_map_200pg files if they exist
    old_fig16 = OUT_DIR / "fig16_capture_map_200pg.png"
    if old_fig16.exists():
        old_fig16.unlink()
        print(f"  Removed old {old_fig16.name}")

    old_fig16_pdf = OUT_DIR / "fig16_capture_map_200pg.pdf"
    if old_fig16_pdf.exists():
        old_fig16_pdf.unlink()
        print(f"  Removed old {old_fig16_pdf.name}")

    fig.savefig(OUT_DIR / "fig16_capture_map_200pg.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig16_capture_map_200pg.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"  Saved proper fig16 to {OUT_DIR}/fig16_capture_map_200pg.*")
    print("\nDone!")


if __name__ == "__main__":
    main()
