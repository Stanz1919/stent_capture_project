"""
physics.external_field
======================
Uniform external magnetic field and field superposition.

Classes
-------
UniformExternalField
    A spatially uniform B0 vector field (∇B0 = 0 everywhere).

TotalField
    Composes a :class:`~stent_capture.core.field_model.StentRing` with an
    optional :class:`UniformExternalField` to produce B_total = B_stent + B0.

Physical background
-------------------
When a uniform external field B0 is applied to the stent geometry:

1. **Stent magnetisation** — For ferromagnetic stents (304 SS, 430 SS) at
   typical clinical field strengths (0.1–1 T), the material is driven close to
   magnetic saturation.  Rather than modelling the nonlinear M(H) hysteresis
   curve, the caller supplies the saturation magnetisation M directly and sets
   ``assume_saturation=True`` on :class:`~stent_capture.core.field_model.StentRing`.
   This flag signals the intent; TotalField respects it but does not alter M.

2. **Superposition** — The Akoun & Yonnet model gives the stent's contribution
   B_stent(r) analytically.  The total field is simply:

       B_total(r) = B0 + B_stent(r)

   where B0 is spatially uniform.

3. **Gradient of the magnitude** — Even though ∇·B0 = 0 and ∇B0 = 0, the
   gradient of the *magnitude* satisfies:

       ∂|B_total|/∂x_i = (B_total · ∂B_stent/∂x_i) / |B_total|

   This is generally ≠ ∂|B_stent|/∂x_i / |B_stent| because |B_total| ≠ |B_stent|
   and B_total is not parallel to B_stent away from the magnetisation axis.
   When B0 is axial (perpendicular to the radial stent magnetisation), B_total
   points mostly along z, so the projection of ∂B_stent/∂x onto B_total is
   suppressed: |∇|B_total|| is *correctly reduced* relative to |∇|B_stent||.
   This is not a numerical error.

4. **Force parameter — the correct capture metric** — The magnetic force on a
   superparamagnetic particle (SPION-labelled cell) with effective volume
   susceptibility χ_eff is:

       F = (V_p · χ_eff / μ₀) · |B_total| · ∇|B_total|
         = (V_p · χ_eff / (2μ₀)) · ∇(|B_total|²)

   The relevant scalar metric for capture distance is therefore the *force
   parameter*:

       FP(r) = |B_total(r)| · |∇|B_total(r)||     [units: T²/m]

   When B0 is applied axially, |B_total| increases from ~30 mT to ~500 mT
   (×17) while |∇|B_total|| decreases by only ~5×, giving a net ~3× gain in
   FP at 200 µm and a ~55× gain at 500 µm.  **Do not use** |∇|B_total||
   alone as a proxy for capture force when B0 ≠ 0.

   Stage 2 will formalise this into a full ``MagneticForce`` class with
   SPION parameters (χ_eff, V_p) and convert FP into a physical force in
   Newtons.  For Stage 1, callers can compute the force parameter directly:

       fp = tf.B_magnitude(x, y, z) * tf.grad_B(x, y, z)
"""

from __future__ import annotations

import numpy as np

from stent_capture.core.field_model import StentRing
from stent_capture.core.gradient import compute_gradient_magnitude


class UniformExternalField:
    """
    A spatially uniform external magnetic field.

    Parameters
    ----------
    B0_vector : array-like, shape (3,)
        Field vector [Bx, By, Bz] in Tesla.
        Example: ``np.array([0, 0, 0.5])`` for 0.5 T along z (axial).

    Examples
    --------
    >>> import numpy as np
    >>> ext = UniformExternalField([0, 0, 0.5])
    >>> pts = np.zeros((4, 3))
    >>> ext.field_at(pts)
    array([[0. , 0. , 0.5],
           [0. , 0. , 0.5],
           [0. , 0. , 0.5],
           [0. , 0. , 0.5]])
    """

    def __init__(self, B0_vector: np.ndarray) -> None:
        self.B0 = np.asarray(B0_vector, dtype=float).ravel()
        if self.B0.shape != (3,):
            raise ValueError("B0_vector must have exactly 3 components.")

    @property
    def magnitude(self) -> float:
        """Scalar magnitude |B0| in Tesla."""
        return float(np.linalg.norm(self.B0))

    def field_at(self, points: np.ndarray) -> np.ndarray:
        """
        Return B0 broadcast to match the observation array.

        Parameters
        ----------
        points : ndarray, shape (N, 3)

        Returns
        -------
        B : ndarray, shape (N, 3)
            Each row is B0 (identical for all points).
        """
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[np.newaxis, :]
        N = pts.shape[0]
        return np.broadcast_to(self.B0, (N, 3)).copy()

    def __repr__(self) -> str:
        return f"UniformExternalField(B0={self.B0} T)"


class TotalField:
    """
    Superposition of a stent ring field and an optional uniform external field.

    B_total(r) = B_stent(r) + B0

    This is the primary field object to pass to force and trajectory
    calculations in later stages.

    Parameters
    ----------
    stent_ring : StentRing
        Magnetised stent geometry.
    external_field : UniformExternalField or None
        If None, TotalField is equivalent to the stent field alone (B0 = 0).

    Notes
    -----
    When ``stent_ring.assume_saturation`` is True the ring's magnetisation M
    is used as-is (the saturation value supplied at construction).  TotalField
    does not modify M; that decision is documented on
    :class:`~stent_capture.core.field_model.StentRing`.
    """

    def __init__(
        self,
        stent_ring: StentRing,
        external_field: UniformExternalField | None = None,
    ) -> None:
        self.stent = stent_ring
        self.external = external_field

    def field_at(self, points: np.ndarray) -> np.ndarray:
        """
        Total B-field at observation points.

        Parameters
        ----------
        points : ndarray, shape (N, 3)
            Observation coordinates [x, y, z] in metres.

        Returns
        -------
        B : ndarray, shape (N, 3)
            B_total = B_stent + B0 in Tesla.
        """
        B = self.stent.field_at(points)
        if self.external is not None:
            B = B + self.external.field_at(points)
        return B

    def B_field(
        self,
        obs_x: np.ndarray,
        obs_y: np.ndarray,
        obs_z: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (Bx, By, Bz) of B_total at observation points."""
        shape = np.broadcast(obs_x, obs_y, obs_z).shape
        pts = np.column_stack([
            np.broadcast_to(obs_x, shape).ravel(),
            np.broadcast_to(obs_y, shape).ravel(),
            np.broadcast_to(obs_z, shape).ravel(),
        ])
        B = self.field_at(pts)
        return (
            B[:, 0].reshape(shape),
            B[:, 1].reshape(shape),
            B[:, 2].reshape(shape),
        )

    def B_magnitude(
        self,
        obs_x: np.ndarray,
        obs_y: np.ndarray,
        obs_z: np.ndarray,
    ) -> np.ndarray:
        """Return |B_total| at observation points."""
        Bx, By, Bz = self.B_field(obs_x, obs_y, obs_z)
        return np.sqrt(Bx**2 + By**2 + Bz**2)

    def grad_B(
        self,
        obs_x: np.ndarray,
        obs_y: np.ndarray,
        obs_z: np.ndarray,
        dx: float = 5e-7,
    ) -> np.ndarray:
        """
        Return |∇|B_total|| via central finite differences.

        Uses ``compute_gradient_magnitude`` with ``self.field_at`` so that
        B0 is correctly included in the gradient computation.
        """
        shape = np.broadcast(obs_x, obs_y, obs_z).shape
        pts = np.column_stack([
            np.broadcast_to(obs_x, shape).ravel(),
            np.broadcast_to(obs_y, shape).ravel(),
            np.broadcast_to(obs_z, shape).ravel(),
        ])
        result = compute_gradient_magnitude(self.field_at, pts, dx=dx)
        return result.reshape(shape)

    def __repr__(self) -> str:
        ext = repr(self.external) if self.external is not None else "None"
        return f"TotalField(stent={self.stent!r}, external={ext})"
