"""
Fig 11 - Static force-balance vs trajectory-integration capture prediction.

Headline comparison figure for Stage 3 / full project.

Both panels show two curves:
  Solid blue  : static capture distance (|F_mag| > |F_drag|, Furlani & Ng 2006).
  Dashed red  : trajectory effective capture range -- outermost injection
                distance d such that a strut-aligned cell injected at
                (r_inner-d, 0, -2 mm) is captured by the stent.
                Found by 7-step bisection (resolution ~11 um).

Panel (a) -- loading sweep : fixed v = 0.2 m/s (MCA mean, Aaslid 1982)
Panel (b) -- velocity sweep: fixed loading set by _M_REF_PG (default 200 pg, Polyak 2008)

Both metrics evaluated along strut-aligned axis (theta=0) for direct
comparison with the static criterion (also evaluated at theta=0, z=0).

Run standalone::

    python -m stent_capture.figures.fig11_static_vs_trajectory
"""

from __future__ import annotations

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from stent_capture.figures.common import DEFAULTS, OUT, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.physics.capture_criterion import capture_distance
from stent_capture.simulation.trajectories import integrate_trajectory



_B0_Z     = 1.5
_R_VES    = 1.54e-3
_V_MCA    = 0.2
_M_REF_PG = 200.0
_N_ITER   = 7

_LOADINGS_PG = np.array([10., 30., 50., 100., 200., 300.])
_VELOCITIES  = np.array([0.02, 0.05, 0.10, 0.20, 0.50, 0.60])

_TRAJ_KW = dict(z_end=2e-3, max_time=1.5, rtol=1e-5, atol=1e-8)

_COLOR_STATIC = "#2980b9"
_COLOR_TRAJ   = "#c0392b"



def _find_range_trajectory(cell, tf, flow, ring, label="", **kw):
    """
    Outermost injection distance d from stent inner surface at which a
    strut-aligned cell (theta=0, injected at (r_inner-d, 0, -2mm)) is
    captured.  Uses 7-step bisection.

    Returns
    -------
    float : effective capture range (m).  0.0 if wall-contact cell escapes;
            r_inner if centreline cell is captured.
    """
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



def _sweep_loading(ring, tf, v_mean):
    """Panel (a): sweep over _LOADINGS_PG at fixed v_mean."""
    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=float(v_mean))
    static_dists, traj_ranges = [], []
    for m_pg in _LOADINGS_PG:
        cell = SPIONLabelledCell(spion_mass_per_cell=float(m_pg) * 1e-15)
        print(f"    loading = {m_pg:.0f} pg", flush=True)
        d_s = capture_distance(cell, tf, flow, direction="inward") * 1e6
        static_dists.append(d_s)
        print(f"      static d = {d_s:.1f} um", flush=True)
        d_t = _find_range_trajectory(
            cell, tf, flow, ring, label=f"{m_pg:.0f}pg", **_TRAJ_KW,
        ) * 1e6
        traj_ranges.append(d_t)
        print(f"      traj range = {d_t:.1f} um", flush=True)
    return static_dists, traj_ranges


def _sweep_velocity(ring, tf, loading_kg):
    """Panel (b): sweep over _VELOCITIES at fixed loading_kg."""
    cell = SPIONLabelledCell(spion_mass_per_cell=float(loading_kg))
    static_dists, traj_ranges = [], []
    for v in _VELOCITIES:
        flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=float(v))
        print(f"    v_mean = {v:.3f} m/s", flush=True)
        d_s = capture_distance(cell, tf, flow, direction="inward") * 1e6
        static_dists.append(d_s)
        print(f"      static d = {d_s:.1f} um", flush=True)
        d_t = _find_range_trajectory(
            cell, tf, flow, ring, label=f"{v:.3f}m/s", **_TRAJ_KW,
        ) * 1e6
        traj_ranges.append(d_t)
        print(f"      traj range = {d_t:.1f} um", flush=True)
    return static_dists, traj_ranges



def make_figure():
    t_start = time.time()

    ring = make_ring(B0_magnitude=_B0_Z)
    ring.assume_saturation = True
    tf = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    print("  Panel (a): loading sweep at v_mean = 0.2 m/s ...", flush=True)
    t0 = time.time()
    static_a, traj_a = _sweep_loading(ring, tf, _V_MCA)
    print(f"  Panel (a) done in {time.time()-t0:.1f} s", flush=True)

    print(f"  Panel (b): velocity sweep at {_M_REF_PG:.0f} pg ...", flush=True)
    t0 = time.time()
    static_b, traj_b = _sweep_velocity(ring, tf, _M_REF_PG * 1e-15)
    print(f"  Panel (b) done in {time.time()-t0:.1f} s", flush=True)

    print(f"  Total computation: {time.time()-t_start:.1f} s", flush=True)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.subplots_adjust(wspace=0.30, left=0.09, right=0.97, top=0.82, bottom=0.13)
    max_y = 200

    ax_a.semilogx(_LOADINGS_PG, static_a, "-o",
                  color=_COLOR_STATIC, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Static  |F$_{mag}$| > |F$_{drag}$|", zorder=5)
    ax_a.semilogx(_LOADINGS_PG, traj_a, "--s",
                  color=_COLOR_TRAJ, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Trajectory  (strut-aligned, binary search)", zorder=5)
    ax_a.fill_between(_LOADINGS_PG, static_a, traj_a,
                      alpha=0.10, color=_COLOR_TRAJ)
    for m, tr in zip(_LOADINGS_PG, traj_a):
        ax_a.annotate(f"{tr:.0f}", xy=(m, tr), xytext=(0, 6),
                      textcoords="offset points", ha="center",
                      fontsize=7, color=_COLOR_TRAJ)
    ax_a.axvline(200, color="gray", ls="--", lw=1.2, alpha=0.8,
                 label="200 pg (Polyak et al. 2008)")
    ax_a.set_xlabel("SPION loading per cell (pg)", fontsize=10)
    ax_a.set_ylabel("Effective capture range ($\\mu$m)", fontsize=10)
    ax_a.set_xlim(7, 380)
    ax_a.set_ylim(0, 200)
    ax_a.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.4g"))
    ax_a.legend(fontsize=7, loc="upper left")
    ax_a.grid(True, which="both", alpha=0.25)
    ax_a.set_title(
        f"(a) Capture range vs SPION loading\n"
        f"$\\bar{{v}}$ = {_V_MCA:.1f} m/s (MCA mean), B$_0$ = {_B0_Z} T",
        fontsize=10,
    )

    ax_b.semilogx(_VELOCITIES, static_b, "-o",
                  color=_COLOR_STATIC, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Static  |F$_{mag}$| > |F$_{drag}$|", zorder=5)
    ax_b.semilogx(_VELOCITIES, traj_b, "--s",
                  color=_COLOR_TRAJ, lw=2.0, ms=7,
                  markerfacecolor="white", markeredgewidth=1.8,
                  label="Trajectory  (strut-aligned, binary search)", zorder=5)
    ax_b.fill_between(_VELOCITIES, static_b, traj_b,
                      alpha=0.10, color=_COLOR_TRAJ)
    for v, tr in zip(_VELOCITIES, traj_b):
        ax_b.annotate(f"{tr:.0f}", xy=(v, tr), xytext=(0, 6),
                      textcoords="offset points", ha="center",
                      fontsize=7, color=_COLOR_TRAJ)
    ax_b.axvline(0.20, color="gray", ls="--", lw=1.2, alpha=0.8,
                 label="0.20 m/s (distal MCA mean, Aaslid et al. 1982)")
    ax_b.set_xlabel("Mean blood velocity $\\bar{v}$ (m/s)", fontsize=10)
    ax_b.set_ylabel("Effective capture range ($\\mu$m)", fontsize=10)
    ax_b.set_xlim(0.014, 0.75)
    ax_b.set_ylim(0, 300)
    ax_b.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3g"))
    ax_b.legend(fontsize=7, loc="upper right")
    ax_b.grid(True, which="both", alpha=0.25)
    ax_b.set_title(
        f"(b) Capture range vs flow velocity\n"
        f"{_M_REF_PG:.0f} pg SPION loading, B$_0$ = {_B0_Z} T",
        fontsize=10,
    )

    return fig, static_a, traj_a, static_b, traj_b



def main():
    print("  Fig 11: Static vs trajectory capture comparison ...")
    t0_total = time.time()
    fig, static_a, traj_a, static_b, traj_b = make_figure()

    fig.savefig(OUT / "fig11_static_vs_trajectory.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig11_static_vs_trajectory.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] fig11_static_vs_trajectory saved  ({time.time()-t0_total:.0f} s total)")

    print(f"\n  Panel (a) -- loading sweep, v = {_V_MCA:.1f} m/s:")
    print(f"  {'loading (pg)':>14}  {'static (um)':>12}  {'traj (um)':>12}  {'ratio':>8}")
    for m, st, tr in zip(_LOADINGS_PG, static_a, traj_a):
        ratio = f"{tr/st:.1f}x" if st > 0.5 else "inf"
        print(f"  {m:>14.0f}  {st:>12.1f}  {tr:>12.1f}  {ratio:>8}")

    print(f"\n  Panel (b) -- velocity sweep, {_M_REF_PG:.0f} pg:")
    print(f"  {'v_mean (m/s)':>14}  {'static (um)':>12}  {'traj (um)':>12}  {'ratio':>8}")
    for v, st, tr in zip(_VELOCITIES, static_b, traj_b):
        ratio = f"{tr/st:.1f}x" if st > 0.5 else "inf"
        print(f"  {v:>14.3f}  {st:>12.1f}  {tr:>12.1f}  {ratio:>8}")

    return static_a, traj_a, static_b, traj_b


if __name__ == "__main__":
    main()
