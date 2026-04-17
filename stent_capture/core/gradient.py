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
B0 relative to the stent field.  When B0 is applied axially (perpendicular to
the radial stent magnetisation), B_total is rotated towards z and the
projection of ∂B_stent/∂x onto B_total is small — this *correctly* reduces
|∇|B_total||.  However, the cell-capture force scales as |B_total|·|∇|B_total||,
not just |∇|B_total|| alone (see physics.external_field.TotalField).

Design note: analytical gradient
---------------------------------
The Akoun & Yonnet kernel could in principle be analytically differentiated
(the arctan and log terms have elementary derivatives).  This was considered
but ruled out for the following reasons:

1. **Numerical sufficiency**: a diagnostic sweep across dx = 1e-8 to 5e-6 m
   with B0 = 0.5 T axial (the hardest cancellation case) shows max/min ratio
   of 1.000 — float64 central differences are well-conditioned here because
   |B_stent| / |B0| ~ 0.06, which still leaves ~10 significant figures in
   the magnitude difference.  See ``TestFDStability`` in tests/test_external_field.py.

2. **Complexity**: the analytical Jacobian of the Akoun & Yonnet kernel
   requires differentiating corner sums of arctan(UV/WR) and ln(U+R) terms
   through the chain rule, producing expressions at least as long as the
   kernel itself, with additional special-case handling near R -> 0.

3. **Future path**: if an analytical gradient is ever needed (e.g. for
   force optimisation in Stage 4), it should be implemented as a separate
   ``_akoun_yonnet_local_gradient`` kernel in core/field_model.py alongside
   the field kernel, not here.
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
    return np.linalg.norm(compute_gradient_vector(field_func, points, dx=dx), axis=1)


def compute_gradient_vector(
    field_func: Callable[[np.ndarray], np.ndarray],
    points: np.ndarray,
    dx: float = 5e-7,
) -> np.ndarray:
    """
    Compute the vector gradient ∇|B(r)| via 3-D central finite differences.

    Parameters
    ----------
    field_func : callable
        ``f(points)`` where *points* has shape ``(N, 3)`` and returns a
        ``(N, 3)`` array of [Bx, By, Bz] values in Tesla.
    points : ndarray, shape (N, 3)
        Observation coordinates [x, y, z] in metres.
    dx : float, optional
        Step size for finite differences (m).  Default 500 nm.

    Returns
    -------
    grad_vec : ndarray, shape (N, 3)
        (∂|B|/∂x, ∂|B|/∂y, ∂|B|/∂z) in T/m at each point.
        ``np.linalg.norm(grad_vec, axis=1)`` recovers ``compute_gradient_magnitude``.

    Notes
    -----
    Required by :func:`~stent_capture.physics.magnetic_force.magnetic_force`
    to obtain both the magnitude and direction of the magnetic force on a
    SPION-labelled cell.  Six field evaluations per call, same cost as
    :func:`compute_gradient_magnitude`.
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        pts = pts[np.newaxis, :]

    offsets = dx * np.eye(3, dtype=float)
    grad = np.zeros_like(pts)   # (N, 3), T/m

    for dim in range(3):
        B_plus  = field_func(pts + offsets[dim])   # (N, 3)
        B_minus = field_func(pts - offsets[dim])   # (N, 3)
        grad[:, dim] = (
            np.linalg.norm(B_plus,  axis=1) -
            np.linalg.norm(B_minus, axis=1)
        ) / (2.0 * dx)

    return grad   # (N, 3), T/m
