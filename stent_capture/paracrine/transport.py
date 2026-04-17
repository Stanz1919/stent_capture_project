"""
paracrine.transport
====================
Two-dimensional reaction–diffusion-advection solver for VEGF transport in the
vessel-wall tissue plane surrounding a magnetised stent.

Model
-----
The tissue is represented as a thin 2-D slab (thickness *h*) wrapping the
inner vessel surface.  Coordinates:

    x — circumferential (periodic, 0 … 2πR)
    z — axial            (zero-flux at boundaries)

VEGF secreted by captured cells diffuses, advects in tissue interstitial flow,
and degrades:

    ∂C/∂t = D (∂²C/∂x² + ∂²C/∂z²)  −  ∇·(u·C)  +  S(x, z)  −  k · C

where *C* is the VEGF concentration (ng mL⁻¹), *D* is the effective
tissue diffusion coefficient, *u* is the tissue interstitial velocity,
*k* is the first-order degradation rate, and *S* is the volumetric source
term from captured cells.

Advection is included to account for tissue interstitial flow. Typical values:
~1 µm/s in normal tissue, 10–100 µm/s in inflamed/pathological tissue.
The Péclet number Pe = |u| L_D / D, where L_D ≈ 730 µm (VEGF diffusion length),
ranges Pe ≈ 7 (normal) to Pe ≈ 350 (inflamed), so advection can range from
secondary to dominant depending on tissue condition. Advection affects both
peak concentration and spatial distribution.

Steady state is obtained by solving the sparse linear system directly
(``scipy.sparse.linalg.spsolve``).  Transient solutions use explicit
Euler with automatic CFL-limited time stepping.

Literature parameters
---------------------
D = 1.04 × 10⁻¹⁰ m²/s  (VEGF₁₆₅ in tissue ECM, in vivo skeletal muscle)
    Mac Gabhann F & Popel AS (2006) PLoS Comput Biol 2(10):e127.

k = 1.93 × 10⁻⁴ s⁻¹    (half-life ≈ 60 min in tissue)
    In vitro (Kleinheinz et al. 2010, Chen et al.; via PLOS ONE 2011 mouse model).

Characteristic diffusion length: L_D = √(D/k) ≈ 734 µm.

Tissue interstitial velocity: u_z ~ 50–100 µm/s (axial, order-of-magnitude).
    Piwnica-Worms et al. (1999); tissue fluid velocity is typically 0.1–10 µm/s
    in normal tissue, can reach 50–100 µm/s in inflamed/pathological tissue.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


# ── Default physical constants ──────────────────────────────────────
D_VEGF_TISSUE  = 1.04e-10   # m²/s  — Mac Gabhann & Popel 2006 (VEGF164 in vivo)
K_DEG_TISSUE   = 1.93e-4    # s⁻¹   — In vitro (Kleinheinz et al. 2010, Chen et al.; via PLOS ONE 2011 mouse model)
L_DIFFUSION    = np.sqrt(D_VEGF_TISSUE / K_DEG_TISSUE)  # ≈ 734 µm


class ParacrineField:
    """
    2-D finite-difference solver for VEGF reaction–diffusion-advection.

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
    u_axial : float or None
        Tissue interstitial velocity in axial direction (m/s).
        Default None = no advection. Typical ~1e-4 m/s (100 µm/s).
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
        u_axial:    float | None = None,
    ) -> None:
        self.Lx, self.Lz = Lx, Lz
        self.Nx, self.Nz = Nx, Nz
        self.D    = D
        self.k_deg = k_deg
        self.periodic_x = periodic_x
        self.u_axial = u_axial

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
        Solve  D ∇²C − u·∇C − k C = −S  for the steady-state concentration.

        Advection term is included if u_axial is set (axial z-direction).
        x is circumferential (periodic if ``periodic_x``), z is axial (zero-flux).

        Returns C(x, z) in ng mL⁻¹, shape (Nx, Nz).
        """
        Nx, Nz = self.Nx, self.Nz
        N   = Nx * Nz
        dx2 = self.dx ** 2
        dz2 = self.dz ** 2
        D   = self.D
        k   = self.k_deg
        u_z = self.u_axial if self.u_axial is not None else 0.0

        # ── 1-D Laplacian in x (circumferential) ──
        if self.periodic_x:
            main_x = np.full(Nx, -2.0)
            off_x  = np.ones(Nx - 1)
            Lx_1d = sparse.diags([off_x, main_x, off_x], [-1, 0, 1],
                                 shape=(Nx, Nx), format='lil')
            Lx_1d[0, Nx - 1] = 1.0
            Lx_1d[Nx - 1, 0] = 1.0
        else:
            main_x = np.full(Nx, -2.0)
            main_x[0] = -1.0       # zero-flux ghost cell
            main_x[-1] = -1.0
            off_x  = np.ones(Nx - 1)
            Lx_1d = sparse.diags([off_x, main_x, off_x], [-1, 0, 1],
                                 shape=(Nx, Nx), format='lil')
        Lx_1d = Lx_1d.tocsr() / dx2

        # ── 1-D Laplacian in z (zero-flux) ──
        main_z = np.full(Nz, -2.0)
        main_z[0]  = -1.0          # zero-flux ghost cell: -2C+C = -C
        main_z[-1] = -1.0
        off_z  = np.ones(Nz - 1)
        Lz_1d = sparse.diags([off_z, main_z, off_z], [-1, 0, 1],
                             shape=(Nz, Nz), format='csr') / dz2

        # ── 1-D advection in z (central difference, zero-flux at ends) ──
        # ∂C/∂z stencil: upper = +1/(2dz), lower = -1/(2dz); boundaries use one-sided
        if u_z != 0.0:
            main_a = np.zeros(Nz)
            upper  = np.full(Nz - 1,  1.0 / (2.0 * self.dz))
            lower  = np.full(Nz - 1, -1.0 / (2.0 * self.dz))
            # Forward diff at j=0: (C[1]-C[0])/dz; backward diff at j=Nz-1
            main_a[0]  = -1.0 / self.dz
            upper[0]   =  1.0 / self.dz
            main_a[-1] =  1.0 / self.dz
            lower[-1]  = -1.0 / self.dz
            Adv_1d = sparse.diags([lower, main_a, upper], [-1, 0, 1],
                                  shape=(Nz, Nz), format='csr')
        else:
            Adv_1d = sparse.csr_matrix((Nz, Nz))

        # ── 2-D operator: L = D*(Lx ⊗ I + I ⊗ Lz) − u_z*(I ⊗ Adv) − k*I ──
        I_x = sparse.eye(Nx, format='csr')
        I_z = sparse.eye(Nz, format='csr')

        A = (D * (sparse.kron(Lx_1d, I_z) + sparse.kron(I_x, Lz_1d))
             - u_z * sparse.kron(I_x, Adv_1d)
             - k * sparse.eye(N, format='csr'))

        A = A.tocsc()
        rhs = -self._source.ravel()
        C = spsolve(A, rhs)
        return C.reshape(Nx, Nz)

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
            Time step.  If None, set to 0.4 × CFL limit (accounts for both
            diffusion and advection).
        C0 : ndarray or None
            Initial concentration.  Defaults to zeros.

        Returns
        -------
        C_snaps : ndarray, shape (n_snapshots, Nx, Nz)
        t_snaps : ndarray, shape (n_snapshots,)
        """
        dx2, dz2 = self.dx ** 2, self.dz ** 2
        u_z = self.u_axial if self.u_axial is not None else 0.0
        # CFL bound: max(|u_z|/dz, D/dz², D/dx²) scaled by dt; includes advection if present
        dt_cfl = 0.5 / (self.D * (1.0 / dx2 + 1.0 / dz2) + abs(u_z) / self.dz)
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

            # Advection term: -u_z * ∂C/∂z (central difference, all x, interior z)
            Adv = np.zeros_like(C)
            if u_z != 0.0:
                Adv[:, 1:-1] = -u_z * (C[:, 2:] - C[:, :-2]) / (2.0 * self.dz)

            C += dt * (D * (Cx + Cz) + Adv + S - k * C)
            np.maximum(C, 0.0, out=C)

        return np.array(snaps), np.array(times)
