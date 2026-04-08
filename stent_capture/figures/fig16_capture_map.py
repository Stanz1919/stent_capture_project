"""
Fig 16 — Capture map: 2D cross-section at z=0, B0 = 0.5 T axial.

1x3 subplots for mean blood velocities 0.05, 0.2, 0.5 m/s.

Each panel:
- Heatmap of capture margin |F_mag| - |F_drag| (pN).
  Green = captured (margin > 0), red = washed away (margin < 0).
- Black contour at margin = 0 (capture envelope).
- White filled rectangles showing stent strut cross-sections.
- Vessel lumen boundary circle.

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
# Strut patch helper
# ---------------------------------------------------------------------------

def _strut_patches():
    """Return matplotlib Rectangle patches for each strut cross-section."""
    R = DEFAULTS["R"]
    w = DEFAULTS["w"]
    t = DEFAULTS["t"]
    n = DEFAULTS["n_struts"]
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    patches = []
    for th in angles:
        # Strut centre in global frame
        cx = R * np.cos(th)
        cy = R * np.sin(th)
        # Local frame: radial = (cos th, sin th), circumferential = (-sin th, cos th)
        # Width along circumferential (w), thickness along radial (t)
        # Bottom-left corner in rotated frame: (-w/2, -t/2) in local → global
        # Build a rotated rectangle: use a Transform
        corner_local = np.array([-w / 2, -t / 2])
        cos_th, sin_th = np.cos(th), np.sin(th)
        # Anchor for Rectangle (lower-left before rotation)
        rect = mpatches.Rectangle(
            (cx - t / 2 * cos_th + w / 2 * sin_th,
             cy - t / 2 * sin_th - w / 2 * cos_th),
            width=t, height=w,
            angle=np.degrees(th),
            rotation_point="xy",
            linewidth=0.5, edgecolor="white", facecolor="white", zorder=3,
        )
        patches.append(rect)
    return patches


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    print("    Building grid and computing magnetic force (once)...")
    xx, yy, pts, inside_lumen = _build_grid()
    F_mag_all = _compute_F_mag(pts)

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    for ax, v_mean in zip(axes, _V_CASES):
        print(f"    Computing drag for v_mean = {v_mean} m/s...")
        F_drag_all = _compute_F_drag(pts, v_mean)

        margin = np.full(len(pts), np.nan)
        margin[inside_lumen] = (F_mag_all - F_drag_all)[inside_lumen]
        margin_grid = margin.reshape(xx.shape)

        # Clip for display: focus on ±5 pN range; larger values compressed
        vmax_plot = 5.0
        margin_clipped = np.clip(margin_grid, -vmax_plot, vmax_plot)

        norm = TwoSlopeNorm(vmin=-vmax_plot, vcenter=0.0, vmax=vmax_plot)
        cmap = plt.get_cmap("RdYlGn")
        im = ax.pcolormesh(xx * 1e3, yy * 1e3, margin_clipped,
                           cmap=cmap, norm=norm, shading="auto", zorder=1)

        # Capture boundary contour
        valid_mask = ~np.isnan(margin_grid)
        if np.any(margin_grid[valid_mask] > 0) and np.any(margin_grid[valid_mask] < 0):
            ax.contour(xx * 1e3, yy * 1e3, margin_grid, levels=[0],
                       colors="black", linewidths=1.5, zorder=2)

        # Vessel lumen boundary
        theta = np.linspace(0, 2 * np.pi, 300)
        R_lumen = (DEFAULTS["R"] - DEFAULTS["t"] / 2)
        ax.plot(R_lumen * 1e3 * np.cos(theta),
                R_lumen * 1e3 * np.sin(theta),
                "w--", lw=1.0, alpha=0.7, zorder=2, label="Lumen boundary")

        # Stent outer circle
        R_outer = (DEFAULTS["R"] + DEFAULTS["t"] / 2)
        ax.plot(R_outer * 1e3 * np.cos(theta),
                R_outer * 1e3 * np.sin(theta),
                "w-", lw=0.5, alpha=0.4, zorder=2)

        # Strut rectangles
        for patch in _strut_patches():
            ax.add_patch(patch)
            patch.set_transform(
                ax.transData +
                plt.matplotlib.transforms.Affine2D().scale(1e3)
            )

        # Actually redraw struts more simply: filled circles at strut centres
        R_s = DEFAULTS["R"]
        t_s = DEFAULTS["t"]
        w_s = DEFAULTS["w"]
        n_s = DEFAULTS["n_struts"]
        for k in range(n_s):
            th_k = 2 * np.pi * k / n_s
            cx = R_s * np.cos(th_k) * 1e3
            cy = R_s * np.sin(th_k) * 1e3
            # Draw as a small white rectangle (approximate as circle for simplicity)
            strut_circle = plt.Circle((cx, cy),
                                       max(t_s, w_s) / 2 * 1e3,
                                       color="white", zorder=4, lw=0)
            ax.add_patch(strut_circle)

        ax.set_aspect("equal")
        ax.set_xlim(-_R_VES * 1.05e3, _R_VES * 1.05e3)
        ax.set_ylim(-_R_VES * 1.05e3, _R_VES * 1.05e3)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)" if ax is axes[0] else "")
        ax.set_title(f"v_mean = {v_mean:.2f} m/s")

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("|F_mag| - |F_drag| (pN)" if ax is axes[-1] else "")
        cbar.ax.axhline(0, color="black", lw=1.0)

    fig.suptitle(
        "Capture map: |F_mag| - |F_drag| (pN),  B0 = 0.5 T axial\n"
        "Green = captured (|F_mag| > |F_drag|),  Red = washed away\n"
        "Black contour = capture boundary (Furlani & Ng 2006 criterion)\n"
        "White circles = stent struts  |  Dashed circle = lumen boundary",
        fontsize=11, y=1.02,
    )
    plt.tight_layout()
    return fig


def main():
    print("  Fig 16: Capture map...")
    fig = make_figure()
    fig.savefig(OUT / "fig16_capture_map.png", dpi=200)
    fig.savefig(OUT / "fig16_capture_map.pdf")
    plt.close(fig)
    print("  [OK] fig16_capture_map saved")


if __name__ == "__main__":
    main()
