"""
physics.magnetic_force
======================
Magnetic force on a SPION-labelled MSC in an external field.

The force on a superparamagnetic particle is derived from the energy of a
magnetic dipole in a non-uniform field (Furlani & Ng 2006, Eq. 2):

    F = mu_0 * V_p * (M_p · nabla) H_total

In the linear susceptibility approximation (M_p = chi_p * H_total, valid for
chi_p << 1 or explicitly modelled for near-saturation via an effective chi):

    F = (V_p * chi_p / mu_0) * (B_total · nabla) B_total

For a spherical particle in a curl-free field (nabla x B = 0 in free space),
this simplifies to the gradient of the magnetic energy density:

    F = (V_p * chi_p / mu_0) * nabla(|B_total|^2 / 2)
      = (V_p * chi_p / mu_0) * |B_total| * nabla|B_total|

The magnitude is:

    |F| = (V_spion * chi_spion / mu_0) * |B_total| * |nabla||B_total||

and the direction is along nabla|B_total| (toward the region of higher |B|,
i.e. toward the strut surface).

This is the standard scalar-force approximation used throughout the magnetic
drug/cell targeting literature (Furlani & Ng 2006, Avilés et al. 2007,
Polyak et al. 2008, Tefft et al. 2014).

References
----------
Furlani, E.P. & Ng, K.C. (2006). Analytical model of magnetic nanoparticle
    transport and capture in the microvasculature. Physical Review E, 73, 061919.
Avilés, M.O. et al. (2007). Theoretical analysis of a transdermal
    ferromagnetic implant for retention of magnetic drug carrier particles.
    Journal of Magnetism and Magnetic Materials, 310, 428-443.
Polyak, B. et al. (2008). High field gradient targeting of magnetic
    nanoparticle-loaded cells to the surfaces of steel stents.
    PNAS, 105(2), 698-703.
Chorny, M. et al. (2007). Targeting stents with locally delivered paclitaxel-
    loaded magnetic nanoparticles prevents neointima formation in porcine
    arteries. FASEB J, 21(3), 2510-2519.
"""

from __future__ import annotations

from numpy import pi
import numpy as np

from stent_capture.core.gradient import compute_gradient_vector

MU_0 = 4 * pi * 1e-7


class SPIONLabelledCell:
    """
    Model of a human bone marrow MSC loaded with superparamagnetic iron oxide
    nanoparticles (SPIONs).

    Default values are representative of human BM-MSCs loaded with
    biodegradable polymeric SPIONs at the dose used by Polyak et al. (2008):

    - Cell radius 10 µm (upper end of 10–20 µm suspended human BM-MSC
      diameter range; Dominici et al. 2006 Cytotherapy — conservative for drag)
    - SPION load 200 pg iron oxide per cell (Polyak et al. 2008 stent study)
    - chi_spion = 2.0 (SI volume susceptibility of magnetite SPIONs in the
      field regime 0.1–1 T; appropriate for partially saturated nanoparticles
      with M ~ chi * B/mu_0 — see Furlani & Ng 2006 Table 1)
    - spion_density 5170 kg/m³ (magnetite, Fe₃O₄)

    Note on chi_spion
    -----------------
    chi_spion = 2.0 is the susceptibility of the raw SPION material (not an
    effective whole-cell value).  The force is computed using V_spion (SPION
    volume only) × chi_spion, which gives the correct scaling.  No volume-
    fraction reduction to a whole-cell chi is applied (Furlani & Ng 2006,
    Eq. 2).

    Parameters
    ----------
    radius : float
        Cell radius (m).  Default 10 µm.
    spion_mass_per_cell : float
        Total mass of SPION material per cell (kg).  Default 200e-15 kg (200 pg).
    spion_susceptibility : float
        SI volume susceptibility of the SPION material (dimensionless).
        Default 2.0.
    spion_density : float
        Density of SPION material (kg/m³).  Default 5170 (magnetite).
    spion_sat_magnetization : float or None
        Saturation magnetization for Langevin saturation model (A/m).
        Default 446e3 (bulk magnetite, Furlani & Ng 2006 Table 1).
        If None, uses constant susceptibility (backward-compatible mode).
    """

    def __init__(
        self,
        radius: float = 10e-6,
        spion_mass_per_cell: float = 200e-15,
        spion_susceptibility: float = 2.0,
        spion_density: float = 5170,
        spion_sat_magnetization: float | None = 446e3,
    ) -> None:
        self.radius = radius
        self.spion_mass = spion_mass_per_cell
        self.spion_susceptibility = spion_susceptibility
        self.spion_density = spion_density
        self.spion_sat_magnetization = spion_sat_magnetization

    @property
    def volume(self) -> float:
        """Cell volume (m³), assumed spherical."""
        return (4.0 / 3.0) * pi * self.radius ** 3

    @property
    def spion_volume(self) -> float:
        """Volume of SPION material within the cell (m³)."""
        return self.spion_mass / self.spion_density

    def __repr__(self) -> str:
        sat_mode = "Langevin" if self.spion_sat_magnetization is not None else "linear"
        return (
            f"SPIONLabelledCell(radius={self.radius*1e6:.1f}um, "
            f"spion_mass={self.spion_mass*1e15:.1f}pg, "
            f"chi={self.spion_susceptibility}, sat={sat_mode})"
        )


def _chi_effective(chi_0: float, M_sat: float | None, B_mag: np.ndarray) -> np.ndarray:
    """
    Compute field-dependent effective susceptibility using the Langevin model.

    For superparamagnetic particles, the magnetization follows the Langevin saturation
    function M(H) = M_sat * L(alpha * H), where L(x) = coth(x) - 1/x is the Langevin
    function. At low fields (x→0), this gives the linear limit M ≈ chi_0 * H.

    Parameters
    ----------
    chi_0 : float
        Initial (low-field) susceptibility at B → 0 (dimensionless).
    M_sat : float or None
        Saturation magnetization of SPION material (A/m).
        If None, returns constant chi_0 (backward-compatible linear mode).
        446e3 A/m for bulk magnetite (Furlani & Ng 2006 Table 1).
    B_mag : ndarray, shape (N,)
        Magnitude of total magnetic field at each evaluation point (T).

    Returns
    -------
    chi_eff : ndarray, shape (N,)
        Effective susceptibility at each point: chi_eff = chi_0 * f_sat(B),
        where f_sat is the Langevin saturation factor (1 at low B, <1 at high B).

    Notes
    -----
    The Langevin model is the classical approach for paramagnetic/superparamagnetic
    particles in ferrofluids (Furlani & Ng 2006; Rosensweig, Ferrofluids, 1985).
    The saturation factor is computed using the Langevin function with numerical
    stability at small ξ via Taylor expansion (ξ < 1e-3).
    """
    if M_sat is None:
        return np.full_like(B_mag, chi_0)

    alpha = 3.0 * chi_0 / M_sat
    xi = alpha * B_mag / MU_0

    small = xi < 1e-3
    f_sat = np.empty_like(xi)

    f_sat[small] = 1.0 - xi[small]**2 / 15.0

    xi_large = xi[~small]
    coth_xi = 1.0 / np.tanh(xi_large)
    L_xi = coth_xi - 1.0 / xi_large
    f_sat[~small] = 3.0 * L_xi / xi_large

    return chi_0 * f_sat


def magnetic_force(
    cell: SPIONLabelledCell,
    total_field,
    points: np.ndarray,
    dx: float = 5e-7,
) -> np.ndarray:
    """
    Magnetic force vector on a SPION-labelled cell at each observation point.

    Uses the scalar-gradient approximation (Furlani & Ng 2006):

        F_vec = (V_spion * chi_spion / mu_0) * |B_total| * nabla|B_total|

    where nabla|B_total| is the vector gradient (pointing toward the strut).

    Parameters
    ----------
    cell : SPIONLabelledCell
        Cell model providing V_spion and chi_spion.
    total_field : TotalField or StentRing
        Any object with a ``field_at(points)`` method returning ``(N, 3)``
        B-vectors.
    points : ndarray, shape (N, 3)
        Observation coordinates [x, y, z] in metres.
    dx : float, optional
        Finite-difference step for gradient (m).  Default 500 nm.

    Returns
    -------
    F : ndarray, shape (N, 3)
        Force vector [Fx, Fy, Fz] in Newtons at each point.
        Positive components indicate force in the +x/+y/+z direction
        (toward the strut for points radially outside it).

    Notes
    -----
    The force magnitude |F| equals
    ``(V_spion * chi / mu_0) * |B_total| * |nabla||B_total||``.

    The force direction is along nabla|B_total|, pointing from low to high |B|
    (i.e. attracted toward the strut surface).  For cells outside the stent
    ring (r > R + t/2), the radial force component will be negative (inward).
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim == 1:
        pts = pts[np.newaxis, :]

    B_vecs  = total_field.field_at(pts)
    B_mag   = np.linalg.norm(B_vecs, axis=1)

    chi_eff = _chi_effective(cell.spion_susceptibility,
                             cell.spion_sat_magnetization, B_mag)

    grad_v  = compute_gradient_vector(total_field.field_at, pts, dx=dx)

    return (cell.spion_volume / MU_0) * chi_eff[:, np.newaxis] * B_mag[:, np.newaxis] * grad_v
