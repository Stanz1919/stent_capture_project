"""
core field model
================
3-D magnetic field model for uniformly magnetised rectangular prisms
using the Akoun & Yonnet (1984) closed-form analytical
expressions — the same formulation implemented in magpylib.
"""

from __future__ import annotations

import itertools

import numpy as np
from numpy import pi

MU_0 = 4 * pi * 1e-7

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

       
        Bx += eps * Mz * (-ln_VpR)
        By += eps * Mz * (-ln_UpR)
        Bz += eps * Mz * atan_UV

       
        Bx += eps * Mx * atan_VW
        By += eps * Mx * (-ln_WpR)
        Bz += eps * Mx * (-ln_VpR)

   
        Bx += eps * My * (-ln_WpR)
        By += eps * My * atan_UW
        Bz += eps * My * (-ln_UpR)


    return (
        (-pre * Bx).reshape(shape),
        (-pre * By).reshape(shape),
        (-pre * Bz).reshape(shape),
    )


class StentRing:


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

    def field_at(self, points: np.ndarray) -> np.ndarray:

        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[np.newaxis, :]
        N = pts.shape[0]

        obs = pts.T

        a = self.t / 2
        b = self.w / 2
        c = self.L / 2
        Mlx, Mly, Mlz = self.M_local

        B = np.zeros((3, N))

        for i in range(self.n_struts):
            Rm = self.rot[i]
            center = np.array([self.cx[i], self.cy[i], self.cz[i]])

            obs_loc = Rm.T @ (obs - center[:, np.newaxis])

            bxl, byl, bzl = _akoun_yonnet_local(
                obs_loc[0], obs_loc[1], obs_loc[2],
                a, b, c, Mlx, Mly, Mlz,
            )

            B_loc = np.stack([bxl.ravel(), byl.ravel(), bzl.ravel()])
            B += Rm @ B_loc

        return B.T

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

        from stent_capture.core.gradient import compute_gradient_magnitude

        shape = np.broadcast(obs_x, obs_y, obs_z).shape
        pts = np.column_stack([
            np.broadcast_to(obs_x, shape).ravel(),
            np.broadcast_to(obs_y, shape).ravel(),
            np.broadcast_to(obs_z, shape).ravel(),
        ])
        result = compute_gradient_magnitude(self.field_at, pts, dx=dx)
        return result.reshape(shape)
