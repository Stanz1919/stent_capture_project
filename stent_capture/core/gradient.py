"""
core.gradient
=============
Finite-difference gradient utilities for scalar magnetic field magnitudes.

The central function ``compute_gradient_magnitude`` is deliberately
field-model-agnostic: it accepts any callable that maps an (N, 3) array of
observation points to an (N, 3) B-field array.  This means it works equally
well with :class:`~stent_capture.core.field_model.StentRing`,
:class:`~stent_capture.physics.external_field.TotalField`, or any future
field source.

Physics note
------------
∇|B| is computed with 6 evaluations of |B| per point (central differences in
x, y, z).  When an external uniform field B0 is present:

    |∇|B_total|| ≠ |∇|B_stent||

even though ∇B0 = 0, because the gradient of the *magnitude* is nonlinear:

    ∂|B_total|/∂x = (B_total · ∂B_total/∂x) / |B_total|
                  = (B_total · ∂B_stent/∂x) / |B_total|

The component of ∂B_stent/∂x along the (non-zero) B_total direction can be
larger or smaller than along B_stent alone, depending on the orientation of
B0 relative to the stent field.
"""

from __future__ import annotations

from typing import Callable

import numpy as np


def compute_gradient_magnitude(
    field_func: Callable[[np.ndarray], np.ndarray],
    points: np.ndarray,
    dx: float = 5e-7,
) -> np.ndarray:
    """
    Compute |∇|B(r)|| via 3-D central finite differences.

    Parameters
    ----------
    field_func : callable
        ``f(points)`` where *points* has shape ``(N, 3)`` and returns a
        ``(N, 3)`` array of [Bx, By, Bz] values in Tesla.
    points : ndarray, shape (N, 3)
        Observation coordinates [x, y, z] in metres.
    dx : float, optional
        Step size for finite differences (m).  Default 500 nm gives a good
        balance between truncation error and round-off for stent-scale fields.

    Returns
    -------
    grad_mag : ndarray, shape (N,)
        |∇|B|| in T/m at each point.

    Notes
    -----
    The three partial derivatives are computed as:

        ∂|B|/∂x ≈ (|B(r + dx·x̂)| − |B(r − dx·x̂)|) / (2·dx)

    and combined as  |∇|B|| = sqrt((∂|B|/∂x)² + (∂|B|/∂y)² + (∂|B|/∂z)²).

    Six field evaluations per call (one +/− pair per spatial dimension).
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        pts = pts[np.newaxis, :]

    offsets = dx * np.eye(3, dtype=float)  # shape (3, 3): one row per axis

    grad_sq = np.zeros(len(pts))
    for dim in range(3):
        B_plus  = field_func(pts + offsets[dim])   # (N, 3)
        B_minus = field_func(pts - offsets[dim])   # (N, 3)
        dBmag = (
            np.linalg.norm(B_plus,  axis=1) -
            np.linalg.norm(B_minus, axis=1)
        ) / (2.0 * dx)
        grad_sq += dBmag**2

    return np.sqrt(grad_sq)
