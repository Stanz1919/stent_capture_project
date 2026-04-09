"""
simulation.capture_efficiency
==============================
Population-level capture efficiency for SPION-labelled cell bundles.

Three public functions:

sweep_injection_line
    Integrate N trajectories with injection points distributed along a 3-D line
    segment; optionally parallelised via ``multiprocessing.Pool``.

capture_efficiency_vs_velocity
    Sweep over blood-flow mean velocities and return the capture fraction at each.

capture_efficiency_vs_loading
    Sweep over SPION loading per cell and return the capture fraction at each.

Parallelisation notes
---------------------
On Windows, ``multiprocessing`` uses the ``spawn`` start method.  Worker
functions must therefore be importable at module level — ``_run_one`` satisfies
this requirement.  A single ``Pool`` is created per outer sweep call to avoid
repeated process-spawn overhead.
"""

from __future__ import annotations

import multiprocessing
import os
from typing import Callable, Sequence

import numpy as np

from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.simulation.trajectories import CellTrajectory, integrate_trajectory


# ---------------------------------------------------------------------------
# Top-level worker — MUST be at module level for multiprocessing pickling
# ---------------------------------------------------------------------------

def _run_one(args: tuple) -> CellTrajectory:
    """Picklable worker: unpack args and call integrate_trajectory."""
    cell, tf, flow, ring, pos, kw = args
    return integrate_trajectory(cell, tf, flow, ring, pos, **kw)


# ---------------------------------------------------------------------------
# Internal helper: build injection points along a line
# ---------------------------------------------------------------------------

def _injection_points(line_start: np.ndarray, line_end: np.ndarray,
                      n_points: int) -> np.ndarray:
    """Return (n_points, 3) array of evenly-spaced positions on the segment."""
    ts = np.linspace(0.0, 1.0, n_points)
    return line_start + ts[:, None] * (line_end - line_start)


def _make_summary(trajectories: list[CellTrajectory],
                  injection_points: np.ndarray) -> dict:
    statuses   = [t.status for t in trajectories]
    n_total    = len(trajectories)
    n_captured = statuses.count('captured')
    n_escaped  = statuses.count('escaped')
    n_error    = statuses.count('error')
    return {
        'n_total':          n_total,
        'n_captured':       n_captured,
        'n_escaped':        n_escaped,
        'n_error':          n_error,
        'efficiency':       n_captured / n_total if n_total > 0 else 0.0,
        'injection_points': injection_points,
    }


# ---------------------------------------------------------------------------
# sweep_injection_line
# ---------------------------------------------------------------------------

def sweep_injection_line(
    cell,
    total_field,
    blood_flow:   BloodFlow,
    stent_ring,
    line_start:   np.ndarray,
    line_end:     np.ndarray,
    n_points:     int = 50,
    n_workers:    int | None = None,
    **trajectory_kwargs,
) -> tuple[list[CellTrajectory], dict]:
    """
    Integrate *n_points* cell trajectories injected along a line segment.

    Parameters
    ----------
    cell : SPIONLabelledCell
    total_field : TotalField
    blood_flow : BloodFlow
    stent_ring : StentRing
    line_start, line_end : array-like, shape (3,)
        Start and end injection points (m).
    n_points : int
        Number of evenly-spaced injection positions.  Default 50.
    n_workers : int or None
        Worker processes for ``multiprocessing.Pool``.
        ``None`` → ``os.cpu_count()``.  ``1`` → serial (no Pool created).
    **trajectory_kwargs
        Forwarded verbatim to :func:`integrate_trajectory`
        (e.g. ``z_end``, ``max_time``, ``rtol``).

    Returns
    -------
    trajectories : list of CellTrajectory
        One entry per injection point, ordered start → end.
    summary : dict
        Keys: ``n_total``, ``n_captured``, ``n_escaped``, ``n_error``,
        ``efficiency`` (float ∈ [0, 1]),
        ``injection_points`` (ndarray, shape (N, 3)).
    """
    line_start = np.asarray(line_start, dtype=float)
    line_end   = np.asarray(line_end,   dtype=float)
    pts = _injection_points(line_start, line_end, n_points)

    args_list = [
        (cell, total_field, blood_flow, stent_ring, pts[i], trajectory_kwargs)
        for i in range(n_points)
    ]

    workers  = n_workers if n_workers is not None else os.cpu_count()
    use_pool = (workers > 1) and (n_points > 4)

    if use_pool:
        try:
            with multiprocessing.Pool(processes=min(workers, n_points)) as pool:
                trajectories = pool.map(_run_one, args_list)
        except Exception:
            trajectories = [_run_one(a) for a in args_list]
    else:
        trajectories = [_run_one(a) for a in args_list]

    return trajectories, _make_summary(trajectories, pts)


# ---------------------------------------------------------------------------
# capture_efficiency_vs_velocity
# ---------------------------------------------------------------------------

def capture_efficiency_vs_velocity(
    cell,
    total_field,
    stent_ring,
    velocities:           Sequence[float],
    line_start:           np.ndarray,
    line_end:             np.ndarray,
    n_cells_per_velocity: int   = 30,
    vessel_radius:        float = 1.54e-3,
    n_workers:            int | None = None,
    **trajectory_kwargs,
) -> dict:
    """
    Sweep over mean blood-flow velocities and compute capture efficiency.

    A fresh :class:`BloodFlow` is created for each velocity using
    *vessel_radius*.  Injection positions are constant across all velocities
    (same *line_start* → *line_end* with *n_cells_per_velocity* points).

    All trajectories for a single velocity are dispatched to a single
    ``multiprocessing.Pool`` to minimise process-spawn overhead.

    Parameters
    ----------
    cell : SPIONLabelledCell
    total_field : TotalField
    stent_ring : StentRing
    velocities : sequence of float
        Mean velocities to test (m/s).
    line_start, line_end : array-like, shape (3,)
        Injection line endpoints (m).
    n_cells_per_velocity : int
        Injection points per velocity.  Default 30.
    vessel_radius : float
        Vessel radius used when constructing each ``BloodFlow`` (m).
        Default 1.54 mm.
    n_workers : int or None
        Passed to the inner :func:`sweep_injection_line` call.
    **trajectory_kwargs
        Forwarded to :func:`integrate_trajectory`.

    Returns
    -------
    dict
        ``velocities``   — ndarray (n_v,) m/s
        ``efficiencies`` — ndarray (n_v,) capture fraction ∈ [0, 1]
        ``n_captured``   — ndarray (n_v,) int
        ``n_total``      — int (same for every velocity)
        ``trajectories`` — list of lists ``[n_v][n_cells]`` CellTrajectory
    """
    velocities = np.asarray(velocities, dtype=float)
    efficiencies     = np.empty(len(velocities))
    n_captured_arr   = np.empty(len(velocities), dtype=int)
    all_trajectories: list[list[CellTrajectory]] = []

    line_start = np.asarray(line_start, dtype=float)
    line_end   = np.asarray(line_end,   dtype=float)

    for i, v in enumerate(velocities):
        print(f"    v_mean = {v:.4f} m/s …", flush=True)
        flow = BloodFlow(vessel_radius=vessel_radius, mean_velocity=float(v))
        trajs, summary = sweep_injection_line(
            cell, total_field, flow, stent_ring,
            line_start, line_end,
            n_points=n_cells_per_velocity,
            n_workers=n_workers,
            **trajectory_kwargs,
        )
        efficiencies[i]   = summary['efficiency']
        n_captured_arr[i] = summary['n_captured']
        all_trajectories.append(trajs)
        print(f"      -> {summary['n_captured']}/{summary['n_total']} captured "
              f"(efficiency = {summary['efficiency']:.3f})", flush=True)

    return {
        'velocities':   velocities,
        'efficiencies': efficiencies,
        'n_captured':   n_captured_arr,
        'n_total':      n_cells_per_velocity,
        'trajectories': all_trajectories,
    }


# ---------------------------------------------------------------------------
# capture_efficiency_vs_loading
# ---------------------------------------------------------------------------

def capture_efficiency_vs_loading(
    cell_factory:        Callable[[float], SPIONLabelledCell],
    total_field,
    blood_flow:          BloodFlow,
    stent_ring,
    loadings_kg:         Sequence[float],
    line_start:          np.ndarray,
    line_end:            np.ndarray,
    n_cells_per_loading: int = 30,
    n_workers:           int | None = None,
    **trajectory_kwargs,
) -> dict:
    """
    Sweep over SPION loading per cell and compute capture efficiency.

    Parameters
    ----------
    cell_factory : callable
        ``cell_factory(loading_kg)`` → :class:`SPIONLabelledCell`.
    total_field : TotalField
    blood_flow : BloodFlow
    stent_ring : StentRing
    loadings_kg : sequence of float
        SPION masses per cell to test (kg).
    line_start, line_end : array-like, shape (3,)
        Injection line endpoints (m).
    n_cells_per_loading : int
        Injection points per loading value.  Default 30.
    n_workers : int or None
        Passed to the inner :func:`sweep_injection_line` call.
    **trajectory_kwargs
        Forwarded to :func:`integrate_trajectory`.

    Returns
    -------
    dict
        ``loadings_kg``  — ndarray (n_L,) kg
        ``loadings_pg``  — ndarray (n_L,) pg (convenience)
        ``efficiencies`` — ndarray (n_L,) capture fraction ∈ [0, 1]
        ``n_captured``   — ndarray (n_L,) int
        ``n_total``      — int (same for every loading)
        ``trajectories`` — list of lists ``[n_L][n_cells]`` CellTrajectory
    """
    loadings_kg = np.asarray(loadings_kg, dtype=float)
    efficiencies     = np.empty(len(loadings_kg))
    n_captured_arr   = np.empty(len(loadings_kg), dtype=int)
    all_trajectories: list[list[CellTrajectory]] = []

    line_start = np.asarray(line_start, dtype=float)
    line_end   = np.asarray(line_end,   dtype=float)

    for i, m in enumerate(loadings_kg):
        print(f"    loading = {m * 1e15:.1f} pg …", flush=True)
        cell = cell_factory(float(m))
        trajs, summary = sweep_injection_line(
            cell, total_field, blood_flow, stent_ring,
            line_start, line_end,
            n_points=n_cells_per_loading,
            n_workers=n_workers,
            **trajectory_kwargs,
        )
        efficiencies[i]   = summary['efficiency']
        n_captured_arr[i] = summary['n_captured']
        all_trajectories.append(trajs)
        print(f"      -> {summary['n_captured']}/{summary['n_total']} captured "
              f"(efficiency = {summary['efficiency']:.3f})", flush=True)

    return {
        'loadings_kg':  loadings_kg,
        'loadings_pg':  loadings_kg * 1e15,
        'efficiencies': efficiencies,
        'n_captured':   n_captured_arr,
        'n_total':      n_cells_per_loading,
        'trajectories': all_trajectories,
    }
