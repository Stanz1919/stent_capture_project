"""
simulation.trajectories
=======================
Single-cell trajectory integration for SPION-labelled MSCs
flowing through a magnetised stent under an applied field.

Physics
-------
For a 10 µm cell at physiological flow velocities, the particle Reynolds
number Re = ρ v d / η ≈ 1.1 (blood, v = 0.2 m/s).  Inertia is negligible:
the cell instantly reaches terminal velocity where viscous (Stokes) drag
balances the magnetic force.  The trajectory satisfies a first-order ODE:

    dr/dt = v_blood(r) + F_mag(r) / (6π η R_cell)

where the first term is the local Poiseuille blood velocity and the second
is the magnetic drift velocity.  Integration uses scipy.integrate.solve_ivp
with method='RK45', adaptive step control, and event-based termination.

Termination events
------------------
1. **Escape**: cell reaches z = z_end (downstream boundary).
2. **Strut proximity**: 3-D distance from cell to any strut centre < strut
   surface radius + capture_tolerance.  Uses a cylindrical approximation
   (radius = max(w, t) / 2) — acceptable for Stage 3a; a point-to-rectangle
   metric may be substituted in a later refinement.
3. **Wall**: radial distance r_xy reaches the stent inner-surface radius
   (R − t/2).  This is the primary capture event for cells inside the lumen
   approaching the stent wall radially.

References
----------
Furlani, E.P. & Ng, K.C. (2006). Analytical model of magnetic nanoparticle
    transport and capture in the microvasculature. Physical Review E, 73, 061919.
Tefft, B.J. et al. (2014). Magnetizable stent-grafts enable endothelial cell
    capture. IEEE Transactions on Magnetics, 50(11), 1-4.
"""

from __future__ import annotations

from numpy import pi
import numpy as np
from scipy.integrate import solve_ivp

from stent_capture.physics.magnetic_force import magnetic_force, SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow


# ---------------------------------------------------------------------------
# Event factories
# ---------------------------------------------------------------------------

def _make_escape_event(z_end: float):
    """
    Terminal event: cell reaches the downstream plane z = z_end.

    The event function z_end − y[2] starts positive and decreases to zero
    as the cell travels in +z.  direction = −1 ensures it fires only on the
    decreasing crossing (i.e., cell moving forward, not backward).
    """
    def event(t, y):
        return z_end - y[2]
    event.terminal  = True
    event.direction = -1
    return event


def _make_strut_proximity_event(stent_ring, capture_tolerance: float):
    """
    Terminal event: cell reaches within capture_tolerance of any strut surface.

    Cylindrical approximation: strut "surface radius" = max(w, t) / 2 = 50 µm
    (default geometry).  For cells inside the lumen approaching the stent wall,
    the wall event (_make_wall_event) typically fires first; this event provides
    a secondary catch for cells approaching struts from the radially outer side.

    Note: the cylindrical approximation may over-estimate the strut cross-section
    along the circumferential direction.  A point-to-rectangle distance metric
    can replace this in a later refinement (flagged in CHANGELOG).
    """
    strut_r = max(stent_ring.w, stent_ring.t) / 2
    cx = stent_ring.cx.copy()
    cy = stent_ring.cy.copy()
    n  = stent_ring.n_struts

    def event(t, y):
        min_dist = np.inf
        for i in range(n):
            dx = y[0] - cx[i]
            dy = y[1] - cy[i]
            dz = y[2]           # strut centres are all at z = 0
            surface_dist = np.sqrt(dx*dx + dy*dy + dz*dz) - strut_r
            if surface_dist < min_dist:
                min_dist = surface_dist
        return min_dist - capture_tolerance
    event.terminal  = True
    event.direction = -1
    return event


def _make_wall_event(lumen_radius: float):
    """
    Terminal event: cell reaches the stent inner-surface radius.

    lumen_radius = R − t/2 is the inner radius of the stent ring.  This is
    the primary capture mechanism for cells flowing inside the lumen: when the
    cell's radial position r_xy = sqrt(x² + y²) reaches lumen_radius, it has
    contacted the stent inner surface and is captured.
    """
    def event(t, y):
        r_xy = np.sqrt(y[0] * y[0] + y[1] * y[1])
        return lumen_radius - r_xy   # positive inside lumen, zero at capture
    event.terminal  = True
    event.direction = -1             # fires when r_xy increases through lumen_radius
    return event


# ---------------------------------------------------------------------------
# CellTrajectory — result object
# ---------------------------------------------------------------------------

class CellTrajectory:
    """
    Integrated path of a single SPION-labelled cell through the stent region.

    Produced by :func:`integrate_trajectory`.

    Attributes
    ----------
    positions : ndarray, shape (N, 3)
        Cell position [x, y, z] at each recorded time step (m).
    times : ndarray, shape (N,)
        Time at each recorded step (s).
    status : str
        Integration outcome:
        ``'captured'``  — cell reached a strut surface or the lumen wall.
        ``'escaped'``   — cell reached z = z_end without capture.
        ``'error'``     — max_time exceeded (timeout) or solver failure.
    cell : SPIONLabelledCell
    total_field : TotalField
    blood_flow : BloodFlow
    """

    def __init__(
        self,
        positions:   np.ndarray,
        times:       np.ndarray,
        status:      str,
        cell:        SPIONLabelledCell,
        total_field,
        blood_flow:  BloodFlow,
    ) -> None:
        self.positions   = positions
        self.times       = times
        self.status      = status
        self.cell        = cell
        self.total_field = total_field
        self.blood_flow  = blood_flow

    @property
    def capture_position(self) -> np.ndarray | None:
        """Position at capture (m), or None if not captured."""
        return self.positions[-1].copy() if self.status == 'captured' else None

    @property
    def capture_time(self) -> float | None:
        """Time at capture (s), or None if not captured."""
        return float(self.times[-1]) if self.status == 'captured' else None

    @property
    def velocities(self) -> np.ndarray:
        """
        Instantaneous cell velocity at each recorded position (m/s).

        v_cell = v_blood(r) + F_mag(r) / (6π η R_cell)

        Computed in one batched pass over all stored positions.

        Returns
        -------
        v : ndarray, shape (N, 3)
        """
        v_blood    = self.blood_flow.velocity_at(self.positions)          # (N, 3)
        F_mag_vec  = magnetic_force(self.cell, self.total_field,
                                    self.positions)                        # (N, 3)
        drag_coeff = 6.0 * pi * self.blood_flow.viscosity * self.cell.radius
        return v_blood + F_mag_vec / drag_coeff                           # (N, 3)

    def __repr__(self) -> str:
        t_ms = self.times[-1] * 1e3 if len(self.times) else 0.0
        return (
            f"CellTrajectory(status={self.status!r}, "
            f"n_steps={len(self.times)}, "
            f"t_final={t_ms:.2f} ms)"
        )


# ---------------------------------------------------------------------------
# integrate_trajectory — public entry point
# ---------------------------------------------------------------------------

def integrate_trajectory(
    cell:              SPIONLabelledCell,
    total_field,
    blood_flow:        BloodFlow,
    stent_ring,
    r_inject:          np.ndarray,
    z_end:             float = 2e-3,
    max_time:          float = 1.0,
    capture_tolerance: float = 5e-6,
    rtol:              float = 1e-6,
    atol:              float = 1e-9,
    max_step:          float = np.inf,
) -> CellTrajectory:
    """
    Integrate a single SPION-labelled cell trajectory through the stent region.

    Uses the terminal-velocity (low-Re) approximation:

        dr/dt = v_blood(r) + F_mag(r) / (6π η R_cell)

    Integration is performed by ``scipy.integrate.solve_ivp`` with
    method='RK45' and event-based termination.

    Parameters
    ----------
    cell : SPIONLabelledCell
        SPION loading, cell radius, susceptibility.
    total_field : TotalField
        Combined stent + external field (``field_at`` interface).
    blood_flow : BloodFlow
        Poiseuille flow model (provides ``velocity_at`` and viscosity).
    stent_ring : StentRing
        Provides strut geometry (cx, cy, R, t, w) for event functions and
        lumen radius R − t/2.
    r_inject : array-like, shape (3,)
        Initial position [x, y, z] in metres.  Must be inside the lumen
        (sqrt(x²+y²) < R − t/2) and upstream of z_end.
    z_end : float
        Downstream escape plane (m).  Default 2 mm.
    max_time : float
        Safety timeout (s).  Integration stops with ``status='error'`` if
        no terminal event fires within this time.  Default 1 s.
    capture_tolerance : float
        Strut-proximity threshold for the strut event (m).  Default 5 µm.
    rtol, atol : float
        Relative / absolute tolerances for RK45.  Defaults 1e-6 / 1e-9.
    max_step : float
        Maximum allowed step size (s).  Default np.inf (fully adaptive).
        Set to a finite value (e.g. 1e-3 s) to force finer output for
        smooth trajectories where RK45 would otherwise take large steps.

    Returns
    -------
    CellTrajectory
        Trajectory with ``status``, ``positions``, ``times``, and references
        to cell / field / flow objects for downstream analysis.
    """
    y0           = np.asarray(r_inject, dtype=float).ravel()
    lumen_radius = stent_ring.R - stent_ring.t / 2
    drag_coeff   = 6.0 * pi * blood_flow.viscosity * cell.radius

    # Build RHS closure (called hundreds of times by RK45)
    def rhs(t, y):
        pos        = y.reshape(1, 3)
        v_blood    = blood_flow.velocity_at(pos)[0]             # (3,)
        F_mag_vec  = magnetic_force(cell, total_field, pos)[0]  # (3,)
        return v_blood + F_mag_vec / drag_coeff

    events = [
        _make_escape_event(z_end),
        _make_strut_proximity_event(stent_ring, capture_tolerance),
        _make_wall_event(lumen_radius),
    ]

    sol = solve_ivp(
        rhs,
        t_span=(0.0, max_time),
        y0=y0,
        method='RK45',
        events=events,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )

    # Determine outcome
    if sol.status == -1:
        # Solver failure
        status = 'error'
    elif sol.status == 0:
        # Reached max_time without any terminal event → timeout
        status = 'error'
    else:
        # sol.status == 1: a terminal event fired
        if len(sol.t_events[0]) > 0:
            status = 'escaped'
        elif len(sol.t_events[1]) > 0 or len(sol.t_events[2]) > 0:
            status = 'captured'
        else:
            status = 'error'   # should not happen

    return CellTrajectory(
        positions=sol.y.T,   # (N, 3)
        times=sol.t,         # (N,)
        status=status,
        cell=cell,
        total_field=total_field,
        blood_flow=blood_flow,
    )
