"""
paracrine.secretion
====================
Map captured-cell positions to a spatial VEGF source field S(x, z) for the
:class:`~stent_capture.paracrine.transport.ParacrineField` solver.

Each captured cell is modelled as a Gaussian source:

    S_i(x, z) = q_vol · exp( −[(x−x_i)² + (z−z_i)²] / (2σ²) )

where *q_vol* is the peak volumetric secretion rate in ng mL⁻¹ s⁻¹ and σ
is the cell radius.  The total VEGF mass secreted per cell is q_cell (g s⁻¹)
distributed over the Gaussian footprint and an effective tissue thickness *h*.

Literature values
-----------------
q_cell  = 0.068 molecules / cell / s  (VEGF₁₆₅)
        = 5.08 × 10⁻²¹ g / s         (MW_VEGF = 45 kDa)
    Stefanini MO et al. (2008) PLoS ONE 3(11):e3565.
    (Calibrated against two-compartment mouse model; originally per
    myonuclear domain, here applied per single endothelial cell.)

Independent measurement:
    0.001 pg / cell / day  (= 1.16 × 10⁻²⁰ g/s)  — retinal endothelial cells.
    Li J et al. (2006) Curr Eye Res 31(4):353–61.

σ  = 10 µm  (endothelial cell radius).
h  = 20 µm  (effective tissue-slab thickness, ≈ 2 cell layers).
"""

from __future__ import annotations

import numpy as np

MW_VEGF  = 45_000.0     # g / mol
N_A      = 6.022e23      # molecules / mol

Q_CELL_MOL_PER_S = 0.068                      # molecules / cell / s  (Stefanini 2008)
Q_CELL_G_PER_S   = Q_CELL_MOL_PER_S * MW_VEGF / N_A   # ≈ 5.08e-21 g/s

SIGMA_DEFAULT     = 10e-6     # m  — cell radius
H_TISSUE_DEFAULT  = 20e-6     # m  — slab thickness


class VEGFSource:
    """
    Build the spatial source field for one or more captured cells.

    Parameters
    ----------
    q_cell : float
        VEGF secretion rate per cell (g s⁻¹).
        Default from Stefanini et al. (2008).
    sigma : float
        Gaussian half-width (m).  Default 10 µm (cell radius).
    h_tissue : float
        Effective tissue thickness (m).  Default 20 µm.
    """

    def __init__(
        self,
        q_cell:   float = Q_CELL_G_PER_S,
        sigma:    float = SIGMA_DEFAULT,
        h_tissue: float = H_TISSUE_DEFAULT,
    ) -> None:
        self.q_cell   = q_cell
        self.sigma    = sigma
        self.h_tissue = h_tissue

    def source_field(
        self,
        X:              np.ndarray,
        Z:              np.ndarray,
        cell_positions: np.ndarray,
    ) -> np.ndarray:
        """
        Compute the volumetric source S(x, z) in ng mL⁻¹ s⁻¹.

        Parameters
        ----------
        X, Z : ndarray, shape (Nx, Nz)
            Mesh-grid arrays of spatial coordinates (m).
        cell_positions : ndarray, shape (N_cells, 2)
            Each row is (x_i, z_i) of a captured cell (m).

        Returns
        -------
        S : ndarray, shape (Nx, Nz)
            Source field in ng mL⁻¹ s⁻¹.
        """
        s  = self.sigma
        h  = self.h_tissue
        # Peak volumetric rate: q_cell / (2π σ² h) in g / (m³ · s)
        # Convert to ng / mL:  1 g/m³ = 10⁹ ng / 10⁶ mL = 10³ ng/mL
        q_vol_peak = (self.q_cell / (2.0 * np.pi * s * s * h)) * 1e3  # ng/(mL·s)

        S = np.zeros_like(X)
        positions = np.atleast_2d(cell_positions)
        for xi, zi in positions:
            r2 = (X - xi) ** 2 + (Z - zi) ** 2
            S += q_vol_peak * np.exp(-r2 / (2.0 * s * s))
        return S

    @staticmethod
    def from_molecules_per_s(rate: float) -> float:
        """Convert a secretion rate in molecules/cell/s to g/s."""
        return rate * MW_VEGF / N_A

    @staticmethod
    def cell_positions_on_ring(
        stent_ring,
        z_positions: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Generate captured-cell positions on the unrolled stent surface.

        Maps the cylindrical strut centres (cx, cy) to circumferential
        coordinate θ → x = R·θ on the flat 2-D domain.

        Parameters
        ----------
        stent_ring : StentRing
            Provides cx, cy, R, n_struts.
        z_positions : ndarray or None
            Axial positions for each strut ring.
            If None, all cells placed at z = 0.

        Returns
        -------
        positions : ndarray, shape (n_struts, 2)
            Columns are (x_circ, z).
        """
        theta = np.arctan2(stent_ring.cy, stent_ring.cx)
        theta = theta % (2.0 * np.pi)
        x_circ = stent_ring.R * theta
        if z_positions is None:
            z_positions = np.zeros_like(x_circ)
        return np.column_stack([x_circ, z_positions])
