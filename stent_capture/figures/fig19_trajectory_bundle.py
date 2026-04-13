"""
Fig 19 — Trajectory bundle: capture efficiency vs mean blood velocity.

Three 3-D panels (one per velocity) showing 20 cell trajectories injected
along a radial line from r = 0.1 mm to r = 1.45 mm at z = −2 mm.  All cells
carry 200 pg SPION loading under B0 = 0.5 T axial field.

Velocities
----------
(a) v̄ = 0.02 m/s  — slow flow, most cells captured
(b) v̄ = 0.10 m/s  — intermediate flow
(c) v̄ = 0.50 m/s  — fast MCA-range flow, few/none captured

Colour coding
-------------
Blue  — captured
Red   — escaped
Green — injection markers

Parameters
----------
- Cell radius: 10 µm, SPION: 200 pg, χ = 2.0
- Stent: 8 struts, R = 1.5 mm, t = 80 µm, w = 100 µm, L = 500 µm, M = 1 MA/m
- B0: 0.5 T axial (+z)
- Vessel radius: 1.54 mm, η = 4 mPa·s
- Injection: 20 points along (0.1 mm, 0, −2 mm) → (1.45 mm, 0, −2 mm)

Run standalone::

    python -m stent_capture.figures.fig19_trajectory_bundle
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numpy import pi

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.simulation.capture_efficiency import sweep_injection_line

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

_B0_Z      = 1.5   # T — MRI-strength, matches COMSOL
_R_VES     = 1.54e-3
_N_CELLS   = 20
_VELOCITIES = [0.02, 0.10, 0.50]   # m/s

_LINE_START = np.array([0.10e-3, 0.0, -2e-3])
_LINE_END   = np.array([1.45e-3, 0.0, -2e-3])   # just inside lumen wall

_COLOR_CAPTURED = "#2980b9"   # blue
_COLOR_ESCAPED  = "#c0392b"   # red
_COLOR_INJECT   = "#27ae60"   # green

# ---------------------------------------------------------------------------
# Strut geometry helpers (local copy — also in fig18)
# ---------------------------------------------------------------------------

_FACE_IDX = [
    [0, 2, 3, 1],
    [4, 5, 7, 6],
    [0, 1, 5, 4],
    [2, 6, 7, 3],
    [0, 4, 6, 2],
    [1, 3, 7, 5],
]


def _strut_faces_mm(ring) -> list:
    """Return list of quad-face arrays (mm) for Poly3DCollection."""
    a = ring.t / 2
    b = ring.w / 2
    c = ring.L / 2
    corners = np.array([
        [sa * a, sb * b, sc * c]
        for sa in (-1, 1) for sb in (-1, 1) for sc in (-1, 1)
    ])   # (8, 3) m
    all_faces = []
    for i in range(ring.n_struts):
        rot    = ring.rot[i]
        centre = np.array([ring.cx[i], ring.cy[i], 0.0])
        verts  = (rot @ corners.T).T + centre    # (8, 3) m
        verts_mm = verts * 1e3
        for fidx in _FACE_IDX:
            all_faces.append(verts_mm[fidx])
    return all_faces


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure():
    ring = make_ring(B0_magnitude=_B0_Z)  # Adaptive M for COMSOL calibration at 1.5 T
    ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    cell = SPIONLabelledCell(spion_mass_per_cell=200e-15)

    R     = DEFAULTS["R"]
    t_str = DEFAULTS["t"]
    L     = DEFAULTS["L"]
    r_lumen_mm = (R - t_str / 2) * 1e3

    # -----------------------------------------------------------------------
    # Run trajectory bundles
    # -----------------------------------------------------------------------
    results = []   # list of (v_mean, trajectories, summary)
    for v in _VELOCITIES:
        print(f"  v_mean = {v:.2f} m/s …", flush=True)
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        trajs, summary = sweep_injection_line(
            cell, tf, flow, ring,
            _LINE_START, _LINE_END,
            n_points=_N_CELLS,
            z_end=2e-3, max_time=1.0,
        )
        print(f"    -> {summary['n_captured']}/{summary['n_total']} captured "
              f"(efficiency = {summary['efficiency']:.2f})", flush=True)
        results.append((v, trajs, summary))

    # -----------------------------------------------------------------------
    # Figure layout — 1 row × 3 columns of 3D axes
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 6))
    panel_labels = ('a', 'b', 'c')
    axes = [fig.add_subplot(1, 3, k + 1, projection='3d') for k in range(3)]

    # Precompute strut faces (same for all panels)
    faces = _strut_faces_mm(ring)
    theta = np.linspace(0, 2 * pi, 80)

    for ax, label, (v, trajs, summary) in zip(axes, panel_labels, results):
        # Strut prisms
        poly = Poly3DCollection(faces, alpha=0.22, facecolor='lightgray',
                                edgecolor='darkgray', linewidth=0.4)
        ax.add_collection3d(poly)

        # Lumen boundary rings
        for z_circ in [-L / 2, 0.0, L / 2]:
            ax.plot(r_lumen_mm * np.cos(theta),
                    r_lumen_mm * np.sin(theta),
                    np.full_like(theta, z_circ * 1e3),
                    'k--', lw=0.5, alpha=0.25)

        # Trajectories
        n_captured = 0
        n_escaped  = 0
        for traj in trajs:
            pos_mm = traj.positions * 1e3
            if traj.status == 'captured':
                ax.plot(pos_mm[:, 0], pos_mm[:, 1], pos_mm[:, 2],
                        color=_COLOR_CAPTURED, lw=1.2, alpha=0.75)
                n_captured += 1
            elif traj.status == 'escaped':
                ax.plot(pos_mm[:, 0], pos_mm[:, 1], pos_mm[:, 2],
                        color=_COLOR_ESCAPED, lw=0.9, alpha=0.35)
                n_escaped += 1
            else:   # 'error' (timeout on very slow near-wall cells)
                ax.plot(pos_mm[:, 0], pos_mm[:, 1], pos_mm[:, 2],
                        color='gray', lw=0.7, alpha=0.2)

        # Injection markers (start of each trajectory)
        inj = summary['injection_points'] * 1e3
        ax.scatter(inj[:, 0], inj[:, 1], inj[:, 2],
                   c=_COLOR_INJECT, s=15, zorder=6, depthshade=False)

        # Axes labels and limits
        ax.set_xlabel('x (mm)', labelpad=3)
        ax.set_ylabel('y (mm)', labelpad=3)
        ax.set_zlabel('z (mm)', labelpad=3)
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.6, 1.6)
        ax.view_init(elev=22, azim=-55)

        eff_pct = summary['efficiency'] * 100
        ax.set_title(
            f'({label}) v̄ = {v:.2f} m/s\n'
            f'{n_captured}/{summary["n_total"]} captured  ({eff_pct:.0f}%)',
            fontsize=10,
        )

    # -----------------------------------------------------------------------
    # Legend patches (shared across panels)
    # -----------------------------------------------------------------------
    import matplotlib.patches as mpatches
    legend_elements = [
        mpatches.Patch(facecolor=_COLOR_CAPTURED, label='Captured (blue)'),
        mpatches.Patch(facecolor=_COLOR_ESCAPED,  label='Escaped (red)'),
        mpatches.Patch(facecolor=_COLOR_INJECT,   label='Injection point'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=9, bbox_to_anchor=(0.5, 0.01))

    fig.suptitle(
        "Fig 19 — Trajectory bundles: 200 pg SPION-labelled cells, B0 = 1.5 T axial (MRI)\n"
        "20 cells injected along a wide radial line (r = 0.1–1.45 mm, z = \u22122 mm) "
        "for visualisation of trajectory shapes across the full lumen.  "
        "Most cells start far from the magnetically active near-wall region — see Fig 20 "
        "for quantitative efficiency in the r = 1.20–1.45 mm near-wall zone.",
        fontsize=9, y=1.01,
    )
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    return fig


def main():
    print("  Fig 19: Trajectory bundle (3 velocities)…")
    fig = make_figure()
    fig.savefig(OUT / "fig19_trajectory_bundle.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig19_trajectory_bundle.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig19_trajectory_bundle saved")


if __name__ == "__main__":
    main()
