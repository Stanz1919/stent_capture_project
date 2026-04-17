"""
core.field_model
================
3-D magnetic field model for uniformly magnetised rectangular prisms
using the Akoun & Yonnet (1984) / Furlani (2001) closed-form analytical
expressions — the same formulation implemented in magpylib.

Public API
----------
StentRing
    Ring of n identical rectangular struts arranged on a circle of radius R
    in the z = 0 plane.  Each strut is modelled as an independent uniformly
    magnetised rectangular prism; their fields are summed in the global frame.

The low-level kernel `_akoun_yonnet_local` is exposed for testing but is
not part of the stable public API.
"""

from __future__ import annotations

import itertools

import numpy as np
from numpy import pi

MU_0 = 4 * pi * 1e-7  # T·m/A


# ---------------------------------------------------------------------------
# Low-level kernel — Akoun & Yonnet analytical expressions
# ---------------------------------------------------------------------------

def _akoun_yonnet_local(
    ox: np.ndarray,
    oy: np.ndarray,
    oz: np.ndarray,
    a: float,
    b: float,
    c: float,
    Mx: float,
    My: float,
    Mz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Exact 3-D B-field from a uniformly magnetised rectangular prism.

    The prism is centred at the origin with half-dimensions (a, b, c) along
    its local (x, y, z) axes.  Magnetisation components (Mx, My, Mz) in A/m.

    Derivation outline
    ------------------
    Surface magnetic charge density σ_m = M·n̂ on each face creates a scalar
    potential φ_m.  The double integral over each rectangular face evaluates
    analytically to:

        F(u, v, w) = u·ln(v+R) + v·ln(u+R) − w·arctan(uv / (wR))

    where R = sqrt(u²+v²+w²).  Evaluating at all 8 corners with alternating
    signs and taking H = −∇φ_m gives closed-form sums of arctan and ln terms.
    Outside the prism B = μ₀H; the interior correction (+ μ₀M along the
    magnetisation axis) is not applied here — callers should mask interior
    points if needed.

    Parameters
    ----------
    ox, oy, oz : array-like
        Observation coordinates in the prism's local frame (m).
    a, b, c    : float
        Half-dimensions along local x, y, z (m).
    Mx, My, Mz : float
        Magnetisation components (A/m).

    Returns
    -------
    Bx, By, Bz : ndarray
        B-field in Tesla, same shape as broadcast(ox, oy, oz).
    """
    shape = np.broadcast(ox, oy, oz).shape
    ox = np.broadcast_to(ox, shape).ravel().astype(float).copy()
    oy = np.broadcast_to(oy, shape).ravel().astype(float).copy()
    oz = np.broadcast_to(oz, shape).ravel().astype(float).copy()

    Bx = np.zeros(ox.size)
    By = np.zeros(ox.size)
    Bz = np.zeros(ox.size)

    pre = MU_0 / (4 * pi)

    for i, j, k in itertools.product(range(2), repeat=3):
        eps = (-1) ** (i + j + k)
        # Corner displacement: index 0 → +half-dim, index 1 → −half-dim
        U = ox + (1 - 2 * i) * a
        V = oy + (1 - 2 * j) * b
        W = oz + (1 - 2 * k) * c
        R = np.sqrt(U**2 + V**2 + W**2)
        R = np.maximum(R, 1e-30)

        atan_UV = np.arctan2(U * V, W * R)
        atan_VW = np.arctan2(V * W, U * R)
        atan_UW = np.arctan2(U * W, V * R)

        ln_UpR = np.log(np.maximum(U + R, 1e-30))
        ln_VpR = np.log(np.maximum(V + R, 1e-30))
        ln_WpR = np.log(np.maximum(W + R, 1e-30))

        # Mz contribution (charges on ±z faces)
        Bx += eps * Mz * (-ln_VpR)
        By += eps * Mz * (-ln_UpR)
        Bz += eps * Mz * atan_UV

        # Mx contribution (charges on ±x faces)
        Bx += eps * Mx * atan_VW
        By += eps * Mx * (-ln_WpR)
        Bz += eps * Mx * (-ln_VpR)

        # My contribution (charges on ±y faces)
        Bx += eps * My * (-ln_WpR)
        By += eps * My * atan_UW
        Bz += eps * My * (-ln_UpR)

    # The raw corner-sum is the NEGATIVE of the physical field (the derivation
    # of the Furlani / Akoun & Yonnet kernel accumulates a global sign flip
    # relative to the H = -∇φ_m convention when the corner signs are
    # distributed as (-1)^(i+j+k)).  Negating here gives the correct B = μ₀H
    # in the exterior: verified against the far-field dipole limit where
    # B_x(along +M axis) → +μ₀m/(2π r³) > 0.
    return (
        (-pre * Bx).reshape(shape),
        (-pre * By).reshape(shape),
        (-pre * Bz).reshape(shape),
    )


# ---------------------------------------------------------------------------
# StentRing — public high-level class
# ---------------------------------------------------------------------------

class StentRing:
    """
    Ring of uniformly magnetised rectangular struts — exact 3-D field via the
    Akoun & Yonnet (1984) analytical expressions.

    Geometry
    --------
    Global frame: stent axis along z; the ring lies in the z = 0 plane.
    Each strut is a rectangular prism with:
      - Radial half-dimension  a = t / 2  (local x → radial direction)
      - Circumferential half-dim b = w / 2  (local y → circumferential)
      - Axial half-dim          c = L / 2  (local z → axial)
    Strut centres: (R cos θ_i, R sin θ_i, 0) for θ_i = 2πi/n_struts.

    Magnetisation
    -------------
    mag_mode='radial'        : M points outward along the radial direction of
                               each strut (local +x).
    mag_mode='circumferential': M points in the circumferential direction (local +y).
    mag_mode='axial'          : M points along the stent axis (local +z).

    Parameters
    ----------
    n_struts : int
    R        : float  — ring radius (m)
    w        : float  — strut circumferential width (m)
    t        : float  — strut radial thickness (m)
    L        : float  — strut axial length (m)
    M        : float  — magnetisation magnitude (A/m)
    mag_mode : str    — 'radial' | 'circumferential' | 'axial'
    assume_saturation : bool
        If True, the magnetisation M passed at construction is treated as the
        saturation value and used directly, regardless of any applied external
        field magnitude.  This is appropriate when the ferromagnetic material
        is expected to saturate under the applied field (e.g. 304 SS or 430 SS
        at B0 ≳ 0.1 T).  The caller is responsible for supplying the correct
        saturation magnetisation for their material:

            304 SS : M_sat ≈ 1.0 MA/m
            430 SS : M_sat ≈ 1.2 MA/m
            Pure Fe: M_sat ≈ 1.7 MA/m

        This flag does NOT implement a full M(H) hysteresis curve — it simply
        signals the modelling intent and can be queried by wrapper classes.
    """

    def __init__(
        self,
        n_struts: int,
        R: float,
        w: float,
        t: float,
        L: float,
        M: float,
        mag_mode: str = "radial",
        assume_saturation: bool = False,
    ) -> None:
        self.n_struts = n_struts
        self.R = R
        self.w = w
        self.t = t
        self.L = L
        self.M = M
        self.mag_mode = mag_mode
        self.assume_saturation = assume_saturation

        self.angles = np.linspace(0, 2 * pi, n_struts, endpoint=False)
        self.cx = R * np.cos(self.angles)
        self.cy = R * np.sin(self.angles)
        self.cz = np.zeros(n_struts)

        # Rotation matrices: v_global = rot[i] @ v_local
        # Columns are local basis vectors expressed in the global frame.
        self.rot: list[np.ndarray] = []
        for th in self.angles:
            self.rot.append(
                np.array(
                    [
                        [np.cos(th), -np.sin(th), 0.0],
                        [np.sin(th),  np.cos(th), 0.0],
                        [0.0,         0.0,        1.0],
                    ]
                )
            )

        if mag_mode == "radial":
            self.M_local = np.array([M, 0.0, 0.0])
        elif mag_mode == "circumferential":
            self.M_local = np.array([0.0, M, 0.0])
        elif mag_mode == "axial":
            self.M_local = np.array([0.0, 0.0, M])
        else:
            raise ValueError(f"Unknown mag_mode: {mag_mode!r}. "
                             "Choose 'radial', 'circumferential', or 'axial'.")

    # ------------------------------------------------------------------
    # Primary interface: field_at(points)
    # ------------------------------------------------------------------

    def field_at(self, points: np.ndarray) -> np.ndarray:
        """
        B-field at an array of observation points.

        Parameters
        ----------
        points : ndarray, shape (N, 3)
            Observation coordinates [x, y, z] in metres.

        Returns
        -------
        B : ndarray, shape (N, 3)
            Magnetic flux density [Bx, By, Bz] in Tesla.
        """
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[np.newaxis, :]
        N = pts.shape[0]

        obs = pts.T  # (3, N)

        a = self.t / 2
        b = self.w / 2
        c = self.L / 2
        Mlx, Mly, Mlz = self.M_local

        B = np.zeros((3, N))

        for i in range(self.n_struts):
            Rm = self.rot[i]
            center = np.array([self.cx[i], self.cy[i], self.cz[i]])

            obs_loc = Rm.T @ (obs - center[:, np.newaxis])  # (3, N)

            bxl, byl, bzl = _akoun_yonnet_local(
                obs_loc[0], obs_loc[1], obs_loc[2],
                a, b, c, Mlx, Mly, Mlz,
            )

            B_loc = np.stack([bxl.ravel(), byl.ravel(), bzl.ravel()])  # (3, N)
            B += Rm @ B_loc

        return B.T  # (N, 3)

    # ------------------------------------------------------------------
    # Convenience wrappers that match the original script's call style
    # ------------------------------------------------------------------

    def B_field(
        self,
        obs_x: np.ndarray,
        obs_y: np.ndarray,
        obs_z: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (Bx, By, Bz) at observation points (any broadcastable shape)."""
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
        """Return |B| at observation points."""
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
        Return |∇|B|| via 3-D central finite differences.

        Uses ``compute_gradient_magnitude`` internally so that the gradient
        calculation is always consistent with the ``field_at`` interface used
        by :class:`~stent_capture.physics.external_field.TotalField`.
        """
        from stent_capture.core.gradient import compute_gradient_magnitude

        shape = np.broadcast(obs_x, obs_y, obs_z).shape
        pts = np.column_stack([
            np.broadcast_to(obs_x, shape).ravel(),
            np.broadcast_to(obs_y, shape).ravel(),
            np.broadcast_to(obs_z, shape).ravel(),
        ])
        result = compute_gradient_magnitude(self.field_at, pts, dx=dx)
        return result.reshape(shape)
