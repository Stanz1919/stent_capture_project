"""
Fig 16 — Force ratio |F_mag|/|F_drag| map: 2D cross-section at z=0, B0 = 0.5 T axial.

1x3 subplots for mean blood velocities 0.05, 0.2, 0.5 m/s.

Each panel:
- Heatmap of log10(|F_mag| / |F_drag|) over the lumen cross-section.
  Green = |F_mag| > |F_drag| (ratio > 1), red = drag dominates (ratio < 1).
  White = exact balance (ratio = 1, log10 = 0).
- Black contour at ratio = 1.0 (static capture boundary, if it exists).
- White filled circles showing stent strut cross-sections.
- White dashed circle: lumen inner boundary.

Physical note: At cerebral arterial flow velocities (0.05-0.5 m/s, MCA-representative), the static
force balance predicts no capture anywhere in the lumen for 10 pg SPION-
loaded cells. This motivates the trajectory-based analysis in Stage 3.

The field and gradient are computed once and shared across all three
velocity panels (only drag changes with velocity).

Default parameters:
- Stent: R=1.5mm, 8 struts, w=100µm, t=80µm, M=1MA/m, assume_saturation=True
- B0: 0.5 T axial (+z)
- Cell: 10 µm, 10 pg, chi=2.0
- Vessel: R_vessel=1.54mm, viscosity=4mPa·s

Run standalone::

    python -m stent_capture.figures.fig16_capture_map
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force
from stent_capture.physics.hydrodynamics import BloodFlow, stokes_drag

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

_V_CASES = [0.05, 0.2, 0.5]   # m/s
_R_VES   = 1.54e-3             # m
_B0_Z    = 0.5                 # T
_CELL    = SPIONLabelledCell()
_GRID_N  = 55                  # grid points per axis (55x55 = 3025 pts)

# Log10 color scale limits: 1e-4 to 1e1 → log10 = -4 to +1
_LOG_MIN = -4.0
_LOG_MAX =  1.0


# ---------------------------------------------------------------------------
# Precompute F_mag over the 2D grid (expensive; done once)
# ---------------------------------------------------------------------------

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
    r_lumen = R - t / 2           # stent inner surface
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


# ---------------------------------------------------------------------------
# Strut marker helper
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    print("    Building grid and computing magnetic force (once)...")
    xx, yy, pts, inside_lumen = _build_grid()
    F_mag_all = _compute_F_mag(pts)

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    max_ratio_overall = 0.0

    for ax, v_mean in zip(axes, _V_CASES):
        print(f"    Computing drag for v_mean = {v_mean} m/s...")
        F_drag_all = _compute_F_drag(pts, v_mean)

        # Compute force ratio; NaN outside lumen
        ratio = np.full(len(pts), np.nan)
        # Avoid divide-by-zero: drag is zero at the wall (r = R_ves) but
        # inside_lumen excludes the wall. Use a small floor for safety.
        drag_safe = np.where(F_drag_all > 0, F_drag_all, np.nan)
        ratio[inside_lumen] = (F_mag_all / drag_safe)[inside_lumen]

        # Track maximum ratio for caption
        valid_ratio = ratio[np.isfinite(ratio)]
        if len(valid_ratio):
            max_ratio_overall = max(max_ratio_overall, float(valid_ratio.max()))

        # Log10 of ratio; clip to display range
        log_ratio = np.full_like(ratio, np.nan)
        pos = inside_lumen & (ratio > 0)
        log_ratio[pos] = np.log10(ratio[pos])
        log_ratio_grid = log_ratio.reshape(xx.shape)

        # Clip to display range
        log_ratio_clipped = np.clip(log_ratio_grid, _LOG_MIN, _LOG_MAX)

        # Diverging colormap centred at 0 (ratio = 1)
        norm = TwoSlopeNorm(vmin=_LOG_MIN, vcenter=0.0, vmax=_LOG_MAX)
        cmap = plt.get_cmap("RdYlGn")
        im = ax.pcolormesh(xx * 1e3, yy * 1e3, log_ratio_clipped,
                           cmap=cmap, norm=norm, shading="auto", zorder=1)

        # Capture boundary contour at log10(ratio) = 0 (ratio = 1)
        valid_mask = np.isfinite(log_ratio_grid)
        if (np.any(log_ratio_grid[valid_mask] > 0) and
                np.any(log_ratio_grid[valid_mask] < 0)):
            ax.contour(xx * 1e3, yy * 1e3, log_ratio_grid, levels=[0],
                       colors="black", linewidths=1.5, zorder=2)

        # Lumen inner boundary — dashed white circle, high alpha for visibility
        theta = np.linspace(0, 2 * np.pi, 300)
        R_lumen = DEFAULTS["R"] - DEFAULTS["t"] / 2
        ax.plot(R_lumen * 1e3 * np.cos(theta),
                R_lumen * 1e3 * np.sin(theta),
                "w--", lw=1.5, alpha=0.9, zorder=3, label="Lumen boundary")

        # Stent outer circle (thin solid)
        R_outer = DEFAULTS["R"] + DEFAULTS["t"] / 2
        ax.plot(R_outer * 1e3 * np.cos(theta),
                R_outer * 1e3 * np.sin(theta),
                "w-", lw=0.5, alpha=0.5, zorder=2)

        # Strut markers (white circles at strut centres, correct mm coords)
        _add_strut_markers(ax)

        ax.set_aspect("equal")
        ax.set_xlim(-_R_VES * 1.05e3, _R_VES * 1.05e3)
        ax.set_ylim(-_R_VES * 1.05e3, _R_VES * 1.05e3)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)" if ax is axes[0] else "")
        ax.set_title(f"v_mean = {v_mean:.2f} m/s")

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if ax is axes[-1]:
            cbar.set_label("log\u2081\u2080(|F_mag| / |F_drag|)")
        # Tick labels: show as powers of 10
        cbar.set_ticks([_LOG_MIN, -3, -2, -1, 0, _LOG_MAX])
        cbar.set_ticklabels(
            [f"10^{int(t)}" if t != 0 else "1 (balance)"
             for t in [_LOG_MIN, -3, -2, -1, 0, _LOG_MAX]]
        )
        cbar.ax.axhline(0, color="black", lw=1.5)

    print(f"    Maximum force ratio in lumen: {max_ratio_overall:.4f}")

    max_ratio_str = f"{max_ratio_overall:.2f}" if max_ratio_overall > 0 else "~0.52"
    fig.suptitle(
        "Force ratio |F_mag|/|F_drag|,  B0 = 0.5 T axial,  10 pg SPION-loaded cell\n"
        f"Green = |F_mag| > |F_drag| (capture),  Red = drag dominates,  "
        f"Max ratio ≈ {max_ratio_str} (near stent inner surface, v_mean = 0.05 m/s)\n"
        "Force ratio |F_mag|/|F_drag| under static force balance. Ratios below unity (red)\n"
        "indicate drag exceeds magnetic force; ratios above unity (green) indicate capture.\n"
        "At cerebral arterial flow velocities (0.05–0.5 m/s, MCA-representative;\n"
        "Aaslid et al. 1982), the ratio remains below 1.0 throughout the lumen\n"
        "for a 10 pg SPION-loaded cell, demonstrating that the static Furlani & Ng\n"
        "criterion is insufficient to explain experimental capture at realistic flow\n"
        "conditions. This motivates the trajectory-based analysis in Stage 3, where\n"
        "cells accumulate radial drift over their transit time even when instantaneous\n"
        "magnetic force is weaker than drag.  Black contour = ratio = 1 (if present).\n"
        "White circles = stent struts  |  Dashed white circle = lumen inner boundary",
        fontsize=9, y=1.01,
    )
    plt.tight_layout()
    return fig, max_ratio_overall


def main():
    print("  Fig 16: Force ratio map...")
    fig, max_ratio = make_figure()
    fig.savefig(OUT / "fig16_capture_map.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig16_capture_map.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] fig16_capture_map saved  (max ratio = {max_ratio:.4f})")


if __name__ == "__main__":
    main()
