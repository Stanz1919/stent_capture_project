"""
paracrine.transport
====================
Two-dimensional reaction–diffusion solver for VEGF transport in the
vessel-wall tissue plane surrounding a magnetised stent.

Model
-----
The tissue is represented as a thin 2-D slab (thickness *h*) wrapping the
inner vessel surface.  Coordinates:

    x — circumferential (periodic, 0 … 2πR)
    z — axial            (zero-flux at boundaries)

VEGF secreted by captured cells diffuses and degrades:

    ∂C/∂t = D (∂²C/∂x² + ∂²C/∂z²)  +  S(x, z)  −  k · C

where *C* is the VEGF concentration (ng mL⁻¹), *D* is the effective
tissue diffusion coefficient, *k* is the first-order degradation rate,
and *S* is the volumetric source term from captured cells.

Steady state is obtained by solving the sparse linear system directly
(``scipy.sparse.linalg.spsolve``).  Transient solutions use explicit
Euler with automatic CFL-limited time stepping.

Literature parameters
---------------------
D = 1.04 × 10⁻¹¹ m²/s  (VEGF₁₆₅ in tissue ECM)
    Mac Gabhann F & Popel AS (2006) PLoS Comput Biol 2(10):e127.

k = 1.93 × 10⁻⁴ s⁻¹    (half-life ≈ 60 min in tissue)
    Stefanini MO et al. (2008) PLoS ONE 3(11):e3565.

Characteristic diffusion length: L_D = √(D/k) ≈ 232 µm.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


# ── Default physical constants ──────────────────────────────────────
D_VEGF_TISSUE  = 1.04e-11   # m²/s  — Mac Gabhann & Popel 2006
K_DEG_TISSUE   = 1.93e-4    # s⁻¹   — Stefanini et al. 2008
L_DIFFUSION    = np.sqrt(D_VEGF_TISSUE / K_DEG_TISSUE)  # ≈ 232 µm


class ParacrineField:
    """
    2-D finite-difference solver for VEGF reaction–diffusion.

    Parameters
    ----------
    Lx, Lz : float
        Domain size (m).  Lx = circumferential, Lz = axial.
    Nx, Nz : int
        Grid points in each direction.
    D : float
        Diffusion coefficient (m²/s).
    k_deg : float
        First-order degradation rate (s⁻¹).
    periodic_x : bool
        If True, x-direction has periodic BCs (default True).
    """

    def __init__(
        self,
        Lx:         float,
        Lz:         float,
        Nx:         int   = 200,
        Nz:         int   = 200,
        D:          float = D_VEGF_TISSUE,
        k_deg:      float = K_DEG_TISSUE,
        periodic_x: bool  = True,
    ) -> None:
        self.Lx, self.Lz = Lx, Lz
        self.Nx, self.Nz = Nx, Nz
        self.D    = D
        self.k_deg = k_deg
        self.periodic_x = periodic_x

        self.dx = Lx / Nx
        self.dz = Lz / Nz
        self.x  = np.linspace(0.5 * self.dx, Lx - 0.5 * self.dx, Nx)
        self.z  = np.linspace(0.5 * self.dz, Lz - 0.5 * self.dz, Nz)
        self.X, self.Z = np.meshgrid(self.x, self.z, indexing='ij')

        self._source = np.zeros((Nx, Nz))

    # ── Source term ─────────────────────────────────────────────────

    @property
    def source(self) -> np.ndarray:
        """Current source field S(x, z) in ng mL⁻¹ s⁻¹."""
        return self._source

    @source.setter
    def source(self, S: np.ndarray) -> None:
        if S.shape != (self.Nx, self.Nz):
            raise ValueError(f"Source shape {S.shape} != grid ({self.Nx}, {self.Nz})")
        self._source = S.copy()

    # ── Steady-state solve ──────────────────────────────────────────

    def solve_steady_state(self) -> np.ndarray:
        """
        Solve  D ∇²C − k C = −S  for the steady-state concentration.

        Returns C(x, z) in ng mL⁻¹, shape (Nx, Nz).
        """
        N  = self.Nx * self.Nz
        dx2 = self.dx ** 2
        dz2 = self.dz ** 2
        D   = self.D
        k   = self.k_deg

        diag   = np.full(N, -2.0 * D / dx2 - 2.0 * D / dz2 - k)
        off_x  = np.full(N, D / dx2)
        off_z  = np.full(N, D / dz2)

        diags  = [diag, off_x[:N-1], off_x[:N-1], off_z[:N-self.Nz], off_z[:N-self.Nz]]
        offsets = [0, 1, -1, self.Nz, -self.Nz]

        A = sparse.diags(diags, offsets, shape=(N, N), format='lil')

        Nz = self.Nz
        Nx = self.Nx

        for i in range(Nx):
            # z-direction zero-flux: mirror at j=0 and j=Nz-1
            idx0 = i * Nz
            A[idx0, idx0] += D / dz2          # ghost = interior → add back
            idx_end = i * Nz + Nz - 1
            A[idx_end, idx_end] += D / dz2

        if self.periodic_x:
            for j in range(Nz):
                i_first = j                     # i=0
                i_last  = (Nx - 1) * Nz + j     # i=Nx-1
                A[i_first, i_last] += D / dx2
                A[i_last, i_first] += D / dx2
        else:
            for j in range(Nz):
                A[j, j] += D / dx2
                A[(Nx-1)*Nz + j, (Nx-1)*Nz + j] += D / dx2

        A = A.tocsc()
        rhs = -self._source.ravel()
        C = spsolve(A, rhs)
        return C.reshape(self.Nx, self.Nz)

    # ── Transient solve ─────────────────────────────────────────────

    def solve_transient(
        self,
        t_final:    float,
        n_snapshots: int   = 50,
        dt:         float | None = None,
        C0:         np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Explicit-Euler time integration from t=0 to t_final.

        Parameters
        ----------
        t_final : float
            End time (s).
        n_snapshots : int
            Number of evenly-spaced snapshots to store.
        dt : float or None
            Time step.  If None, set to 0.4 × CFL limit.
        C0 : ndarray or None
            Initial concentration.  Defaults to zeros.

        Returns
        -------
        C_snaps : ndarray, shape (n_snapshots, Nx, Nz)
        t_snaps : ndarray, shape (n_snapshots,)
        """
        dx2, dz2 = self.dx ** 2, self.dz ** 2
        dt_cfl = 0.5 / (self.D * (1.0 / dx2 + 1.0 / dz2) + 0.5 * self.k_deg)
        if dt is None:
            dt = 0.4 * dt_cfl
        if dt > dt_cfl:
            raise ValueError(f"dt={dt:.3e} exceeds CFL limit {dt_cfl:.3e}")

        n_steps   = int(np.ceil(t_final / dt))
        dt        = t_final / n_steps
        save_every = max(1, n_steps // n_snapshots)

        C = np.zeros((self.Nx, self.Nz)) if C0 is None else C0.copy()
        snaps, times = [], []
        D, k, S = self.D, self.k_deg, self._source
        per_x = self.periodic_x

        for step in range(n_steps + 1):
            if step % save_every == 0 or step == n_steps:
                snaps.append(C.copy())
                times.append(step * dt)

            # Laplacian — x direction
            Cx = np.empty_like(C)
            Cx[1:-1, :] = (C[2:, :] - 2.0 * C[1:-1, :] + C[:-2, :]) / dx2
            if per_x:
                Cx[0, :]  = (C[1, :] - 2.0 * C[0, :] + C[-1, :]) / dx2
                Cx[-1, :] = (C[0, :] - 2.0 * C[-1, :] + C[-2, :]) / dx2
            else:
                Cx[0, :]  = (C[1, :] - C[0, :]) / dx2
                Cx[-1, :] = (C[-2, :] - C[-1, :]) / dx2

            # Laplacian — z direction (zero-flux)
            Cz = np.empty_like(C)
            Cz[:, 1:-1] = (C[:, 2:] - 2.0 * C[:, 1:-1] + C[:, :-2]) / dz2
            Cz[:, 0]    = (C[:, 1] - C[:, 0]) / dz2
            Cz[:, -1]   = (C[:, -2] - C[:, -1]) / dz2

            C += dt * (D * (Cx + Cz) + S - k * C)
            np.maximum(C, 0.0, out=C)

        return np.array(snaps), np.array(times)
