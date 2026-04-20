"""
physics.capture_criterion
=========================
Capture criterion: |F_mag| > |F_drag| (Furlani & Ng 2006 scalar comparison).

A cell at position r is considered captured if the magnitude of the magnetic
force toward the strut exceeds the magnitude of the hydrodynamic drag:

    captured  iff  |F_mag(r)| > |F_drag(r)|

This is the standard conservative criterion from the magnetic targeting
literature (Furlani & Ng 2006, Eq. 14).  It is conservative because the
forces are not coaxial in general: the magnetic force is primarily radial
(toward the nearest strut) while the Stokes drag is axial (along the flow
direction).  A directional criterion would allow some capture even when
|F_mag| < |F_drag| because the radial component alone only needs to exceed
zero to pull the cell in.  The scalar criterion used here therefore
under-estimates the true capture zone; Stage 3 will resolve this with full
trajectory integration.

References
----------
Furlani, E.P. & Ng, K.C. (2006). Analytical model of magnetic nanoparticle
    transport and capture in the microvasculature. Physical Review E, 73, 061919.
Tefft, B.J. et al. (2014). Magnetizable stent-grafts enable endothelial cell
    capture. IEEE Transactions on Magnetics, 50(11), 1-4.
"""

from __future__ import annotations

import numpy as np

from stent_capture.physics.magnetic_force import magnetic_force, SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow, stokes_drag


def capture_map(
    cell: SPIONLabelledCell,
    total_field,
    blood_flow: BloodFlow,
    points: np.ndarray,
    dx: float = 5e-7,
) -> dict:
    """
    Evaluate the capture criterion at an array of observation points.

    Parameters
    ----------
    cell : SPIONLabelledCell
    total_field : TotalField or StentRing
        Field source with ``field_at(points)`` interface.
    blood_flow : BloodFlow
    points : ndarray, shape (N, 3)
        Observation coordinates in metres.
    dx : float, optional
        FD step for magnetic force gradient (m).

    Returns
    -------
    result : dict with keys:
        ``F_mag_vec``   (N, 3) N — magnetic force vectors
        ``F_drag_vec``  (N, 3) N — Stokes drag force vectors
        ``F_mag``       (N,)   N — |F_mag|
        ``F_drag``      (N,)   N — |F_drag|
        ``margin``      (N,)   N — |F_mag| - |F_drag|, positive = captured
        ``captured``    (N,)   bool — True where |F_mag| > |F_drag|
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        pts = pts[np.newaxis, :]

    F_mag_vec  = magnetic_force(cell, total_field, pts, dx=dx)
    F_drag_vec = stokes_drag(cell, blood_flow, pts)

    F_mag_mag  = np.linalg.norm(F_mag_vec,  axis=1)
    F_drag_mag = np.linalg.norm(F_drag_vec, axis=1)

    margin   = F_mag_mag - F_drag_mag
    captured = margin > 0.0

    return {
        "F_mag_vec":  F_mag_vec,
        "F_drag_vec": F_drag_vec,
        "F_mag":      F_mag_mag,
        "F_drag":     F_drag_mag,
        "margin":     margin,
        "captured":   captured,
    }


def capture_distance(
    cell: SPIONLabelledCell,
    total_field,
    blood_flow: BloodFlow,
    direction: str = "radial",
    n_points: int = 500,
    dx: float = 5e-7,
) -> float:
    """
    Outermost distance from the stent surface at which |F_mag| >= |F_drag|.

    Sweeps radially along the through-strut direction (+x) and returns the
    largest distance at which the cell is still captured.

    Parameters
    ----------
    cell : SPIONLabelledCell
    total_field : TotalField or StentRing
    blood_flow : BloodFlow
    direction : str
        ``'radial'`` — legacy outward sweep from stent outer surface along +x
        (note: since r_outer = R_vessel, these points are outside the vessel
        and experience zero drag; use ``'inward'`` for lumen analysis).
        ``'inward'`` — sweeps from the stent inner surface (R - t/2) toward
        the vessel centre, covering the lumen region where cells actually flow.
        This is the physically correct sweep for capture-distance analysis.
    n_points : int
        Number of sample points along the sweep (default 500).
    dx : float, optional
        FD step for gradient (m).

    Returns
    -------
    d_capture : float
        Distance (m) from the reference surface at which the capture criterion
        is last satisfied.  For ``'inward'``, this is the distance from the
        stent inner surface into the lumen.  Returns 0.0 if no capture.
    """
    from stent_capture.figures.common import DEFAULTS

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]

    if direction == "inward":
        r_inner = R - t / 2
        d = np.linspace(1e-6, r_inner * 0.999, n_points)
        x_vals = r_inner - d
        pts = np.column_stack([x_vals, np.zeros(n_points), np.zeros(n_points)])
    else:
        r_outer = R + t / 2
        d = np.linspace(1e-6, 1.5e-3, n_points)
        pts = np.column_stack([d + r_outer, np.zeros_like(d), np.zeros_like(d)])

    result = capture_map(cell, total_field, blood_flow, pts, dx=dx)
    mask = result["captured"]

    if not np.any(mask):
        return 0.0
    return float(d[mask][-1])
