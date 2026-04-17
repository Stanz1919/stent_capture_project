"""
paracrine.therapeutic
======================
Therapeutic-zone analysis for VEGF concentration fields produced by
:class:`~stent_capture.paracrine.transport.ParacrineField`.

Key metric: the *therapeutic zone radius* — the maximum distance from the
nearest captured cell at which the steady-state VEGF concentration still
exceeds a biologically meaningful threshold.

Literature thresholds
---------------------
C_therapeutic ∈ [5, 25] ng mL⁻¹  — Engineering choice; not directly from literature.
    Ozawa et al. (2004) reports cell VEGF production rates (5-70 ng/10⁶ cells/day
    for normal angiogenesis), not tissue microenvironmental concentrations (ng/mL).
    The 5-25 ng/mL used here is an order-of-magnitude approximation for the
    paracrine simulation domain.

C_aberrant    > 100  ng mL⁻¹      — Engineering choice; no direct literature source.
    Related to Ozawa's finding that ≥100 ng/10⁶ cells/day production triggers
    aberrant angiogenesis, but this is cell production rate, not tissue concentration.

    Ozawa CR et al. (2004) J Clin Invest 113(4):516–27.
    doi: 10.1172/JCI18420

Note on residence-time shortcuts
---------------------------------
A formula ``t_res = distance / v_slip`` has been proposed as a quick proxy for
therapeutic exposure.  While dimensionally consistent ([s] = [m] / [m s⁻¹]),
it is **physically invalid** in this context: VEGF transport is
diffusion-dominated (Péclet number Pe = v·L/D ≪ 1 for L < 1 mm and v ~
µm/s slip velocities).  The residence time of a passively advected parcel is
irrelevant when the concentration field is set by the diffusion–reaction
balance.  This module solves the PDE instead.
"""

from __future__ import annotations

import numpy as np


C_THERAPEUTIC_LOW  = 5.0    # ng/mL  — Engineering approximation (not directly from literature)
C_THERAPEUTIC_HIGH = 25.0   # ng/mL  — Order-of-magnitude choice for paracrine domain
C_ABERRANT         = 100.0  # ng/mL  — Related to Ozawa threshold but different units (ng/10⁶ cells/day)


def therapeutic_zone_radius(
    C:              np.ndarray,
    X:              np.ndarray,
    Z:              np.ndarray,
    cell_positions: np.ndarray,
    threshold:      float = C_THERAPEUTIC_LOW,
) -> float:
    """
    Maximum distance from the nearest captured cell where C ≥ threshold.

    Parameters
    ----------
    C : ndarray, shape (Nx, Nz)
        Steady-state VEGF concentration (ng mL⁻¹).
    X, Z : ndarray, shape (Nx, Nz)
        Mesh-grid coordinate arrays (m).
    cell_positions : ndarray, shape (N_cells, 2)
        Captured-cell positions (x, z) in m.
    threshold : float
        Concentration threshold (ng mL⁻¹).  Default 5 ng/mL.

    Returns
    -------
    r_max : float
        Therapeutic zone radius (m).  Zero if threshold never exceeded.
    """
    above = C >= threshold
    if not np.any(above):
        return 0.0

    pts = np.column_stack([X[above], Z[above]])
    cells = np.atleast_2d(cell_positions)

    min_dists = np.full(len(pts), np.inf)
    for cx, cz in cells:
        d = np.sqrt((pts[:, 0] - cx) ** 2 + (pts[:, 1] - cz) ** 2)
        np.minimum(min_dists, d, out=min_dists)

    return float(np.max(min_dists))


def concentration_vs_distance(
    C:              np.ndarray,
    X:              np.ndarray,
    Z:              np.ndarray,
    cell_positions: np.ndarray,
    n_bins:         int = 80,
    r_max:          float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Azimuthally-averaged C(r) profile, where r = distance from nearest cell.

    Parameters
    ----------
    C, X, Z : ndarray, shape (Nx, Nz)
    cell_positions : ndarray, shape (N_cells, 2)
    n_bins : int
    r_max : float or None
        Maximum radial distance.  Defaults to 3 × L_diffusion.

    Returns
    -------
    r_centres : ndarray, shape (n_bins,)
        Bin centres (m).
    C_mean : ndarray, shape (n_bins,)
        Mean concentration in each radial bin (ng mL⁻¹).
    """
    from stent_capture.paracrine.transport import L_DIFFUSION

    cells = np.atleast_2d(cell_positions)
    flat_x = X.ravel()
    flat_z = Z.ravel()
    flat_c = C.ravel()

    min_dist = np.full(len(flat_x), np.inf)
    for cx, cz in cells:
        d = np.sqrt((flat_x - cx) ** 2 + (flat_z - cz) ** 2)
        np.minimum(min_dist, d, out=min_dist)

    if r_max is None:
        r_max = 3.0 * L_DIFFUSION

    edges = np.linspace(0.0, r_max, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    C_mean  = np.zeros(n_bins)

    for b in range(n_bins):
        mask = (min_dist >= edges[b]) & (min_dist < edges[b + 1])
        if np.any(mask):
            C_mean[b] = np.mean(flat_c[mask])

    return centres, C_mean


def time_to_threshold(
    C_snaps:        np.ndarray,
    t_snaps:        np.ndarray,
    X:              np.ndarray,
    Z:              np.ndarray,
    probe_positions: np.ndarray,
    threshold:      float = C_THERAPEUTIC_LOW,
) -> np.ndarray:
    """
    First time at which C ≥ threshold at specified probe points.

    Parameters
    ----------
    C_snaps : ndarray, shape (n_snap, Nx, Nz)
        Transient concentration snapshots.
    t_snaps : ndarray, shape (n_snap,)
        Times of snapshots (s).
    X, Z : ndarray, shape (Nx, Nz)
        Mesh-grid coordinate arrays (m).
    probe_positions : ndarray, shape (n_probes, 2)
        Probe locations (x, z) in m.
    threshold : float
        Concentration threshold (ng mL⁻¹).

    Returns
    -------
    t_thresh : ndarray, shape (n_probes,)
        Time to reach threshold (s).  ``np.inf`` if never reached.
    """
    probes = np.atleast_2d(probe_positions)
    n_probes = len(probes)
    t_thresh = np.full(n_probes, np.inf)

    Nx, Nz = X.shape
    dx = X[1, 0] - X[0, 0] if Nx > 1 else 1.0
    dz = Z[0, 1] - Z[0, 0] if Nz > 1 else 1.0

    for p, (px, pz) in enumerate(probes):
        ix = int(np.clip(np.round((px - X[0, 0]) / dx), 0, Nx - 1))
        iz = int(np.clip(np.round((pz - Z[0, 0]) / dz), 0, Nz - 1))
        for s, t in enumerate(t_snaps):
            if C_snaps[s, ix, iz] >= threshold:
                t_thresh[p] = t
                break

    return t_thresh
