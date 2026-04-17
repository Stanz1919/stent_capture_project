"""
physics.hydrodynamics
=====================
Hydrodynamic drag on a cell in pulsatile / Poiseuille blood flow.

Models the viscous drag on a spherical cell (Stokes regime) moving through
a Newtonian fluid with a Poiseuille velocity profile.  This is the dominant
competing force to magnetic capture at the stent surface.

Default parameters represent a cerebral arterial vessel (M1-segment MCA-
representative):
- Vessel radius 1.54 mm (matching stent outer surface; ~3 mm diameter)
- Mean velocity 0.2 m/s (representative of distal/diseased-vessel flow;
  Aaslid et al. 1982 report healthy MCA mean ~0.62 m/s)
  Womersley number ~ 2, quasi-steady Poiseuille is a reasonable
  approximation for time-averaged analysis.
- Blood viscosity 4 mPa·s (whole blood, hematocrit ~40%)

References
----------
Aaslid, R. et al. (1982). Noninvasive transcranial Doppler ultrasound
    recording of flow velocity in basal cerebral arteries.
    Journal of Neurosurgery, 57(6), 769-774.
Womersley, J.R. (1955). Method for the calculation of velocity, rate of flow
    and viscous drag in arteries when the pressure gradient is known.
    Journal of Physiology, 127, 553-563.
Pedley, T.J. (1980). The Fluid Mechanics of Large Blood Vessels.
    Cambridge University Press.
Furlani, E.P. & Ng, K.C. (2006). Analytical model of magnetic nanoparticle
    transport and capture in the microvasculature. Physical Review E, 73, 061919.
Tefft, B.J. et al. (2014). Magnetizable stent-grafts enable endothelial cell
    capture. IEEE Transactions on Magnetics, 50(11), 1-4.
"""

from __future__ import annotations

from numpy import pi
import numpy as np


class BloodFlow:
    """
    Laminar Poiseuille flow in a cylindrical vessel, flowing in the +z direction.

    The velocity profile is parabolic:

        v_z(r) = 2 * v_mean * (1 - (r / R_vessel)^2)   for r <= R_vessel
        v_z(r) = 0                                       for r >  R_vessel

    where r = sqrt(x^2 + y^2).

    Geometry note
    -------------
    The stent ring (R = 1.5 mm, radial thickness t = 80 µm) sits embedded in
    the vessel wall with its outer surface at r = R + t/2 = 1.54 mm.  Setting
    ``vessel_radius = 1.54e-3`` m aligns the Poiseuille boundary with the stent
    outer surface.  Blood fills the lumen r < R - t/2 ≈ 1.46 mm; in the thin
    annular region 1.46 < r < 1.54 mm (strut metal) the velocity is zero in
    reality but the Poiseuille formula gives a small non-zero value — this is
    acceptable for the near-wall drag estimate used in the capture criterion.

    Parameters
    ----------
    vessel_radius : float
        Vessel radius (m).  Default 1.54e-3 m (matches stent outer surface).
    mean_velocity : float
        Cross-sectional mean velocity (m/s).  Default 0.2 m/s (distal/diseased
        vessel representative; healthy MCA mean is ~0.62 m/s per Aaslid et al. 1982).
    viscosity : float
        Dynamic viscosity of blood (Pa·s).  Default 4e-3 (whole blood).
    density : float
        Blood density (kg/m³).  Default 1060.
    """

    def __init__(
        self,
        vessel_radius: float = 1.54e-3,
        mean_velocity: float = 0.2,
        viscosity: float = 4e-3,
        density: float = 1060,
    ) -> None:
        self.vessel_radius = vessel_radius
        self.mean_velocity = mean_velocity
        self.viscosity = viscosity
        self.density = density
        self.v_max = 2.0 * mean_velocity   # Poiseuille peak velocity at centre

    # ------------------------------------------------------------------
    # Velocity and shear rate profiles
    # ------------------------------------------------------------------

    def velocity_at(self, points: np.ndarray) -> np.ndarray:
        """
        Blood velocity vector at each observation point.

        Parameters
        ----------
        points : ndarray, shape (N, 3)
            Observation coordinates [x, y, z] in metres.

        Returns
        -------
        v : ndarray, shape (N, 3)
            Velocity [vx=0, vy=0, vz] in m/s.  Non-zero only inside the
            vessel (r < R_vessel).
        """
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[np.newaxis, :]

        r = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
        v = np.zeros_like(pts)
        inside = r < self.vessel_radius
        v[inside, 2] = self.v_max * (1.0 - (r[inside] / self.vessel_radius) ** 2)
        return v   # (N, 3)

    def shear_rate_at(self, points: np.ndarray) -> np.ndarray:
        """
        Magnitude of the radial shear rate |dv_z/dr| at each point (1/s).

        For Poiseuille flow: |dv_z/dr| = 2 * v_max * r / R_vessel^2.

        Parameters
        ----------
        points : ndarray, shape (N, 3)

        Returns
        -------
        gamma : ndarray, shape (N,)  —  shear rate in s⁻¹.
        """
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[np.newaxis, :]

        r = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
        gamma = np.zeros(len(pts))
        inside = r < self.vessel_radius
        gamma[inside] = (2.0 * self.v_max * r[inside]) / self.vessel_radius ** 2
        return gamma   # (N,)

    @property
    def wall_shear_stress(self) -> float:
        """
        Wall shear stress at r = R_vessel (Pa).

        tau_wall = viscosity * |dv_z/dr|_{r=R} = viscosity * 4 * v_mean / R_vessel.
        """
        return self.viscosity * 4.0 * self.mean_velocity / self.vessel_radius

    def __repr__(self) -> str:
        return (
            f"BloodFlow(R={self.vessel_radius*1e3:.2f}mm, "
            f"v_mean={self.mean_velocity:.3f}m/s, "
            f"eta={self.viscosity*1e3:.1f}mPas)"
        )


# ---------------------------------------------------------------------------
# Stokes drag
# ---------------------------------------------------------------------------

def stokes_drag(
    cell,
    blood_flow: BloodFlow,
    points: np.ndarray,
    cell_velocities: np.ndarray | None = None,
) -> np.ndarray:
    """
    Stokes drag force on a spherical cell in Poiseuille blood flow.

        F_drag = 6 * pi * eta * R_cell * (v_blood - v_cell)

    For stationary cells (``cell_velocities=None``), this gives the maximum
    drag force the flow can exert on a cell at each location.

    Parameters
    ----------
    cell : SPIONLabelledCell
        Cell model; only ``cell.radius`` is used.
    blood_flow : BloodFlow
        Flow model providing local velocity.
    points : ndarray, shape (N, 3)
        Observation coordinates in metres.
    cell_velocities : ndarray, shape (N, 3), optional
        Velocity of each cell (m/s).  If None, cells are assumed stationary
        (v_cell = 0), giving the worst-case drag for capture.

    Returns
    -------
    F_drag : ndarray, shape (N, 3)
        Drag force vector in Newtons.  The direction is along
        (v_blood - v_cell), i.e. along +z for stationary cells in Poiseuille
        flow — opposing capture which is primarily radial.
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        pts = pts[np.newaxis, :]

    v_blood = blood_flow.velocity_at(pts)   # (N, 3)
    if cell_velocities is None:
        v_cell = np.zeros_like(v_blood)
    else:
        v_cell = np.asarray(cell_velocities, dtype=float)
        if v_cell.ndim == 1:
            v_cell = np.broadcast_to(v_cell, pts.shape).copy()

    coeff = 6.0 * pi * blood_flow.viscosity * cell.radius   # N·s/m
    return coeff * (v_blood - v_cell)   # (N, 3), Newtons
