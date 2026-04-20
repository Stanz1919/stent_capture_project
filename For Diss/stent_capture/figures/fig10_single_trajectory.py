"""
Fig 10 — Single-cell trajectory comparison.

Two trajectories under identical flow (v_mean = 0.05 m/s, B0 = 0.5 T axial)
differing only in SPION loading:
  - Primary:   200 pg (blue)  — captured
  - Secondary:  10 pg (red)   — escaped

Three-panel figure:

(a) 3D view: stent struts as grey semi-transparent prisms (Poly3DCollection),
    primary trajectory (blue line + time markers), secondary (red line, no
    markers).  Injection point green; capture point red star.

(b) Side view (x-z projection): same trajectories, stent band shaded grey,
    lumen inner boundary as dashed line, flow direction arrows.

(c) Time series of primary trajectory only (3 stacked axes, shared t-axis):
    - radial position r(t) in um
    - |v_cell| in mm/s
    - |F_mag| in pN

Parameters
----------
- Injection: (1.3 mm, 0, −2 mm)
- Cell radius: 10 um
- SPION susceptibility: 2.0 (magnetite)
- Stent: 8 struts, R = 1.5 mm, t = 80 um, w = 100 um, L = 500 um, M = 1 MA/m
- B0: 0.5 T axial (+z)
- Vessel radius: 1.54 mm, η = 4 mPa·s

Run standalone::

    python -m stent_capture.figures.fig10_single_trajectory
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numpy import pi

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.figures.style import (
    COLORS_CODE_DEFAULT, COLORS_THRESHOLD_HIGH, COLORS_THRESHOLD,
    COLORS_MARKER_SOURCE
)
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.simulation.trajectories import integrate_trajectory


_B0_Z     = 1.5
_R_VES    = 1.54e-3
_V_MEAN   = 0.05
_R_INJECT = np.array([1.3e-3, 0.0, -2e-3])
_Z_END    = 2e-3

_COLOR_PRIMARY   = COLORS_CODE_DEFAULT
_COLOR_SECONDARY = COLORS_THRESHOLD_HIGH
_COLOR_INJECT    = COLORS_THRESHOLD
_COLOR_CAPTURE   = COLORS_MARKER_SOURCE



_FACE_IDX = [
    [0, 2, 3, 1],
    [4, 5, 7, 6],
    [0, 1, 5, 4],
    [2, 6, 7, 3],
    [0, 4, 6, 2],
    [1, 3, 7, 5],
]


def _strut_faces_mm(ring) -> list:
    """
    Return a list of quad-face arrays (mm) for all struts, ready for
    Poly3DCollection.  Each face is shape (4, 3) in mm.
    """
    a = ring.t / 2
    b = ring.w / 2
    c = ring.L / 2
    corners_local = np.array([
        [sa * a, sb * b, sc * c]
        for sa in (-1, 1) for sb in (-1, 1) for sc in (-1, 1)
    ])
    all_faces = []
    for i in range(ring.n_struts):
        rot    = ring.rot[i]
        centre = np.array([ring.cx[i], ring.cy[i], 0.0])
        verts  = (rot @ corners_local.T).T + centre
        verts_mm = verts * 1e3
        for fidx in _FACE_IDX:
            all_faces.append(verts_mm[fidx])
    return all_faces



def _run_trajectories():
    ring = make_ring(B0_magnitude=_B0_Z)
    ring.assume_saturation = True
    tf   = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MEAN)

    cell_200 = SPIONLabelledCell(spion_mass_per_cell=200e-15)
    cell_10  = SPIONLabelledCell(spion_mass_per_cell=10e-15)

    print("    Integrating 200 pg trajectory (primary — expect capture)...")
    traj_p = integrate_trajectory(
        cell_200, tf, flow, ring, _R_INJECT,
        z_end=_Z_END, max_time=5.0,
    )
    print(f"    -> {traj_p}")

    print("    Integrating 10 pg trajectory (secondary — expect escape)...")
    traj_s = integrate_trajectory(
        cell_10, tf, flow, ring, _R_INJECT,
        z_end=_Z_END, max_time=5.0,
    )
    print(f"    -> {traj_s}")

    if traj_p.status != 'captured':
        print(f"    WARNING: primary trajectory status = {traj_p.status!r} "
              "(expected 'captured')")
    if traj_s.status != 'escaped':
        print(f"    WARNING: secondary trajectory status = {traj_s.status!r} "
              "(expected 'escaped')")

    return ring, tf, flow, cell_200, traj_p, traj_s



def make_figure():
    ring, tf, flow, cell_200, traj_p, traj_s = _run_trajectories()

    R     = DEFAULTS["R"]
    t_str = DEFAULTS["t"]
    L     = DEFAULTS["L"]
    r_lumen = (R - t_str / 2) * 1e3

    fig = plt.figure(figsize=(13, 7))
    gs_outer = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0],
                                wspace=0.30, left=0.07, right=0.97,
                                top=0.88, bottom=0.10)
    ax_side = fig.add_subplot(gs_outer[0])
    gs_ts   = gs_outer[1].subgridspec(2, 1, hspace=0.08)
    ax_r    = fig.add_subplot(gs_ts[0])
    ax_f    = fig.add_subplot(gs_ts[1], sharex=ax_r)

    ps = traj_s.positions * 1e3
    pp = traj_p.positions * 1e3

    r_outer_mm = (R + t_str / 2) * 1e3
    L_mm       = L * 1e3
    stent_band = mpatches.Rectangle(
        (r_lumen, -L_mm / 2), r_outer_mm - r_lumen, L_mm,
        linewidth=0, facecolor='lightgray', alpha=0.6, zorder=1,
        label='Stent wall',
    )
    ax_side.add_patch(stent_band)

    ax_side.axhline(0, color='k', lw=0.3, alpha=0.2)
    ax_side.axvline(r_lumen, color='k', ls='--', lw=1.0, alpha=0.7,
                    label=f'Lumen boundary ({r_lumen:.2f} mm)')

    ax_side.plot(ps[:, 0], ps[:, 2],
                 color=_COLOR_SECONDARY, lw=1.5, alpha=0.85,
                 label='10 pg — escaped')

    ax_side.plot(pp[:, 0], pp[:, 2],
                 color=_COLOR_PRIMARY, lw=2.5, label='200 pg — captured')

    ax_side.scatter([pp[0, 0]], [pp[0, 2]], c=_COLOR_INJECT,
                    s=60, marker='o', zorder=5)
    if traj_p.status == 'captured':
        ax_side.scatter([pp[-1, 0]], [pp[-1, 2]], c=_COLOR_CAPTURE,
                        s=80, marker='o', zorder=5, label='Capture')

    for z_arr in [-1.5, -1.0, -0.5]:
        ax_side.annotate('', xy=(0.85, z_arr + 0.2), xytext=(0.85, z_arr),
                         arrowprops=dict(arrowstyle='->', color='gray',
                                         lw=1.0, mutation_scale=10))

    ax_side.set_xlabel('x (mm)')
    ax_side.set_ylabel('z (mm)')
    ax_side.set_xlim(0.8, r_outer_mm + 0.1)
    z_min = min(traj_p.positions[:, 2].min(), traj_s.positions[:, 2].min()) * 1e3
    z_max = max(traj_p.positions[:, 2].max(), traj_s.positions[:, 2].max()) * 1e3
    ax_side.set_ylim(z_min - 0.1, z_max + 0.1)
    ax_side.legend(fontsize=7)
    ax_side.set_title('(a) Side view (x-z projection)', fontsize=10)
    ax_side.grid(True, alpha=0.25)

    t_ms = traj_p.times * 1e3

    r_vals = np.sqrt(traj_p.positions[:, 0]**2 +
                     traj_p.positions[:, 1]**2) * 1e6
    ax_r.plot(t_ms, r_vals, color=_COLOR_PRIMARY, lw=2)
    ax_r.axhline(r_lumen * 1e3, color='k', ls='--', lw=0.8, alpha=0.6,
                 label=f'Lumen boundary ({r_lumen*1e3:.0f} um)')
    ax_r.set_ylabel('r (um)')
    ax_r.legend(fontsize=6)
    ax_r.set_title('(b) Primary trajectory time series (200 pg, captured)',
                   fontsize=10)
    ax_r.tick_params(labelbottom=False)
    ax_r.grid(True, alpha=0.25)

    print("    Computing |F_mag| time series for primary trajectory...")
    F_mag_vals = np.linalg.norm(
        magnetic_force(cell_200, tf, traj_p.positions), axis=1
    ) * 1e12
    ax_f.semilogy(t_ms, np.maximum(F_mag_vals, 1e-4),
                  color=_COLOR_PRIMARY, lw=2)
    ax_f.set_ylabel('|F_mag| (pN)')
    ax_f.set_xlabel('Time (ms)')
    ax_f.grid(True, which='both', alpha=0.25)

    t_capture = traj_p.capture_time * 1e3 if traj_p.status == 'captured' else None
    for ax in (ax_r, ax_f):
        if t_capture is not None:
            ax.axvline(t_capture, color=_COLOR_CAPTURE, ls=':', lw=1.2,
                       alpha=0.8, label='Capture' if ax is ax_r else None)
        ax.set_xlim(0, t_ms[-1] * 1.02)

    cap = traj_p.capture_position
    cap_str = (f"x={cap[0]*1e3:.3f} mm, z={cap[2]*1e3:.3f} mm"
               if cap is not None else "not captured")
    t_cap_str = (f"{traj_p.capture_time*1e3:.1f} ms"
                 if traj_p.capture_time is not None else "—")

    fig.suptitle(
        "Single-cell trajectory comparison — B0 = 1.5 T axial (MRI), v_mean = 0.05 m/s "
        "(MCA-representative), injection at (1.3 mm, 0, −2 mm)\n"
        "Two cells, identical conditions, differing only in SPION loading.  "
        f"200 pg cell (blue): captured at {cap_str} after {t_cap_str}.  "
        "10 pg cell (red, Polyak 2008 default): escapes downstream.\n"
        "The 200 pg cell accumulates enough radial drift over its ~60 ms transit to reach "
        "a strut; |F_mag| rises sharply in the final ~100 µs as it approaches the surface.",
        fontsize=9, y=0.995,
    )
    return fig


def main():
    print("  Fig 10: Single-cell trajectory comparison...")
    fig = make_figure()
    fig.savefig(OUT / "fig10_single_trajectory.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig10_single_trajectory.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  [OK] fig10_single_trajectory saved")


if __name__ == "__main__":
    main()
