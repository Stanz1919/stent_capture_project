"""
Magnetic Stent Field Gradient — Complete Analysis & Figures
============================================================
Self-contained script for computing and visualising B-field gradients
around magnetised stent struts for cell capture applications.

Physics (2D): Magnetic surface charge model for near-field,
              dipole approximation for far-field. 2D cross-section.

Physics (3D): Akoun & Yonnet (1984) / Furlani (2001) closed-form
              analytical expressions for uniformly magnetised rectangular
              prisms. Exact solution for B at any 3D observation point.
              Same formulation used in magpylib.
"""

import numpy as np
from numpy import pi, sqrt
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from pathlib import Path
import sys

MU_0 = 4 * pi * 1e-7  # T·m/A

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)

# Dissertation-quality plot style
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 11,
    'axes.labelsize': 13, 'axes.titlesize': 14,
    'legend.fontsize': 9, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.grid': True, 'grid.alpha': 0.3, 'lines.linewidth': 1.8,
})

# Default stent parameters
DEFAULTS = dict(
    R=1.5e-3, w=100e-6, t=80e-6, L=500e-6,
    M=1.0e6, n_struts=8, mag_mode='radial'
)

THRESHOLDS = {'40 T/m': 40, '100 T/m': 100, '300 T/m': 300}
TH_COLORS = ['#2ecc71', '#e74c3c', '#3498db']


# ===========================================================
# Core physics: B-field from uniformly magnetised rectangle
# ===========================================================

def B_rect_charge_2D(obs_x, obs_y, cx, cy, w, h, Mx, My, n_q=200):
    """
    2D B-field from a uniformly magnetised rectangle using surface charges.
    σ_m = M·n̂ on each face, discretised into line charges.
    
    obs_x, obs_y: observation point arrays (any shape, will be flattened)
    cx, cy: rectangle centre
    w, h: width (x), height (y)
    Mx, My: magnetisation components
    n_q: charges per face
    """
    shape = obs_x.shape
    x = obs_x.ravel()
    y = obs_y.ravel()
    Bx = np.zeros_like(x)
    By = np.zeros_like(y)
    
    # Four faces: left, right, bottom, top
    # Outward normals: (-1,0), (+1,0), (0,-1), (0,+1)
    # Surface charge σ = M·n̂
    
    faces = [
        # (start_x, start_y, end_x, end_y, normal_x, normal_y)
        (cx - w/2, cy - h/2, cx - w/2, cy + h/2, -1,  0),  # left
        (cx + w/2, cy - h/2, cx + w/2, cy + h/2, +1,  0),  # right
        (cx - w/2, cy - h/2, cx + w/2, cy - h/2,  0, -1),  # bottom
        (cx - w/2, cy + h/2, cx + w/2, cy + h/2,  0, +1),  # top
    ]
    
    for sx, sy, ex, ey, nx, ny in faces:
        sigma = Mx * nx + My * ny  # surface charge density
        if abs(sigma) < 1e-20:
            continue
        
        # Discretise face into point charges
        ts = np.linspace(0, 1, n_q + 1)
        ts = (ts[:-1] + ts[1:]) / 2  # midpoints
        
        qx = sx + (ex - sx) * ts
        qy = sy + (ey - sy) * ts
        
        # Face length
        face_len = sqrt((ex - sx)**2 + (ey - sy)**2)
        dl = face_len / n_q
        
        # Charge strength per element (2D: line charge, units A/m * m = A)
        q = sigma * dl  # equivalent charge strength
        
        for i in range(n_q):
            dx = x - qx[i]
            dy = y - qy[i]
            r2 = dx**2 + dy**2
            r2 = np.maximum(r2, 1e-30)
            
            # 2D: H from line charge = σ*dl/(2π) * r̂/r  (scalar potential gradient)
            # B = μ₀*H outside magnet
            # For surface charge model: H = σ/(2π) * r̂/|r| per unit length
            factor = MU_0 * q / (2 * pi * r2)
            Bx += factor * dx
            By += factor * dy
    
    return Bx.reshape(shape), By.reshape(shape)


def B_dipole_2D(obs_x, obs_y, cx, cy, mx, my):
    """2D magnetic dipole field. m = M * volume * direction."""
    dx = obs_x - cx
    dy = obs_y - cy
    r2 = dx**2 + dy**2
    r2 = np.maximum(r2, 1e-30)
    r = np.sqrt(r2)
    
    m_dot_r = mx * dx + my * dy
    
    # 2D dipole: B = μ₀/(2π) * (2(m·r̂)r̂ - m) / r²
    prefactor = MU_0 / (2 * pi * r2)
    Bx = prefactor * (2 * m_dot_r * dx / r2 - mx)  # wait, let me use 3D
    By = prefactor * (2 * m_dot_r * dy / r2 - my)
    
    # Actually use 3D dipole projected to 2D (more physical for finite-length strut)
    # B = μ₀/(4π) * (3(m·r̂)r̂ - m) / r³
    r3 = r2 * r
    prefactor = MU_0 / (4 * pi * r3)
    m_dot_rhat = m_dot_r / r
    Bx = prefactor * (3 * m_dot_rhat * dx / r - mx)
    By = prefactor * (3 * m_dot_rhat * dy / r - my)
    
    return Bx, By


class StentRing2D:
    """Ring of magnetised rectangular struts in the xy-plane."""
    
    def __init__(self, n_struts, R, w, t, M, mag_mode='radial', use_dipole=False):
        self.n_struts = n_struts
        self.R = R
        self.w = w  # strut width (circumferential)
        self.t = t  # strut thickness (radial)
        self.M = M
        self.mag_mode = mag_mode
        self.use_dipole = use_dipole
        
        self.angles = np.linspace(0, 2*pi, n_struts, endpoint=False)
        self.strut_cx = R * np.cos(self.angles)
        self.strut_cy = R * np.sin(self.angles)
        
        # Magnetisation direction per strut
        if mag_mode == 'radial':
            self.Mx = M * np.cos(self.angles)
            self.My = M * np.sin(self.angles)
        elif mag_mode == 'circumferential':
            self.Mx = M * (-np.sin(self.angles))
            self.My = M * np.cos(self.angles)
        else:  # axial — no in-plane component, skip
            self.Mx = np.zeros_like(self.angles)
            self.My = np.zeros_like(self.angles)
    
    def B_field(self, obs_x, obs_y, n_q=150):
        """Compute total B-field from all struts."""
        Bx_tot = np.zeros_like(obs_x)
        By_tot = np.zeros_like(obs_y)
        
        for i in range(self.n_struts):
            if self.use_dipole:
                # Dipole moment: m = M * area (2D: per unit length)
                area = self.w * self.t
                mx = self.Mx[i] * area
                my = self.My[i] * area
                bx, by = B_dipole_2D(obs_x, obs_y, 
                                      self.strut_cx[i], self.strut_cy[i],
                                      mx, my)
            else:
                bx, by = B_rect_charge_2D(
                    obs_x, obs_y,
                    self.strut_cx[i], self.strut_cy[i],
                    self.w, self.t,
                    self.Mx[i], self.My[i],
                    n_q=n_q
                )
            Bx_tot += bx
            By_tot += by
        
        return Bx_tot, By_tot
    
    def B_magnitude(self, obs_x, obs_y, **kw):
        Bx, By = self.B_field(obs_x, obs_y, **kw)
        return np.sqrt(Bx**2 + By**2)
    
    def grad_B(self, obs_x, obs_y, dx=5e-7, **kw):
        """Compute |∇|B|| using central differences."""
        Bp_x = self.B_magnitude(obs_x + dx, obs_y, **kw)
        Bm_x = self.B_magnitude(obs_x - dx, obs_y, **kw)
        Bp_y = self.B_magnitude(obs_x, obs_y + dx, **kw)
        Bm_y = self.B_magnitude(obs_x, obs_y - dx, **kw)
        
        dBdx = (Bp_x - Bm_x) / (2 * dx)
        dBdy = (Bp_y - Bm_y) / (2 * dx)
        
        return np.sqrt(dBdx**2 + dBdy**2)


def make_ring(**overrides):
    """Create a StentRing2D with defaults + overrides."""
    p = {**DEFAULTS, **overrides}
    return StentRing2D(p['n_struts'], p['R'], p['w'], p['t'], p['M'],
                       p.get('mag_mode', 'radial'),
                       p.get('use_dipole', False))


# ===========================================================
# 3D physics: Akoun & Yonnet / Furlani closed-form expressions
# ===========================================================

def _B_cuboid_local(ox, oy, oz, a, b, c, Mx, My, Mz):
    """
    Exact 3D B-field from a uniformly magnetised rectangular prism using the
    Akoun & Yonnet (1984) / Furlani (2001) analytical expressions — the same
    formulation implemented in magpylib.

    The prism is centred at the origin with half-dimensions (a, b, c) along
    the local (x, y, z) axes.  Magnetisation components (Mx, My, Mz) in A/m.

    Derivation outline
    ------------------
    Surface magnetic charge density σ_m = M·n̂ on each face creates a scalar
    potential φ_m.  The double integral over each rectangular face evaluates
    analytically to:

        F(u, v, w) = u·ln(v+R) + v·ln(u+R) − w·arctan(uv / (wR))

    where R = sqrt(u²+v²+w²).  Evaluating at the 8 corners with alternating
    signs and taking H = −∇φ_m gives closed-form sums of arctan and ln terms.
    Outside the prism B = μ₀H; the inside correction (+ μ₀M along mag. axis)
    is not applied here — mask interior points in callers if needed.

    Parameters
    ----------
    ox, oy, oz : array-like — observation coords in local frame (m)
    a, b, c    : float      — half-dimensions along local x, y, z (m)
    Mx, My, Mz : float      — magnetisation components (A/m)

    Returns
    -------
    Bx, By, Bz : ndarray — field in Tesla, same shape as inputs
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
        # Corner displacements: i=0 → +half-dim, i=1 → −half-dim
        U = ox + (1 - 2 * i) * a
        V = oy + (1 - 2 * j) * b
        W = oz + (1 - 2 * k) * c
        R = np.sqrt(U ** 2 + V ** 2 + W ** 2)
        R = np.maximum(R, 1e-30)

        # arctan terms — np.arctan2 is well-defined when denominator is zero
        atan_UV = np.arctan2(U * V, W * R)   # parallel: Mz → Bz
        atan_VW = np.arctan2(V * W, U * R)   # parallel: Mx → Bx
        atan_UW = np.arctan2(U * W, V * R)   # parallel: My → By

        # ln terms — argument (coord + R) ≥ 0 since R ≥ |coord|;
        # clamp to avoid log(0) on degenerate edges
        ln_UpR = np.log(np.maximum(U + R, 1e-30))
        ln_VpR = np.log(np.maximum(V + R, 1e-30))
        ln_WpR = np.log(np.maximum(W + R, 1e-30))

        # Mz contribution: charges on ±z faces
        Bx += eps * Mz * (-ln_VpR)
        By += eps * Mz * (-ln_UpR)
        Bz += eps * Mz * atan_UV

        # Mx contribution: charges on ±x faces
        Bx += eps * Mx * atan_VW
        By += eps * Mx * (-ln_WpR)
        Bz += eps * Mx * (-ln_VpR)

        # My contribution: charges on ±y faces
        Bx += eps * My * (-ln_WpR)
        By += eps * My * atan_UW
        Bz += eps * My * (-ln_UpR)

    return (pre * Bx).reshape(shape), (pre * By).reshape(shape), (pre * Bz).reshape(shape)


class StentRing3D:
    """
    Ring of magnetised rectangular struts — exact 3D field via
    Akoun & Yonnet analytical expressions.

    Coordinate convention
    ---------------------
    Global frame: stent axis along z, ring in the z = 0 plane.
    Each strut has a local frame where:
      local x → radial direction   (cos θ, sin θ, 0)
      local y → circumferential    (−sin θ, cos θ, 0)
      local z → axial              (0, 0, 1)
    Half-dimensions in local frame: (t/2, w/2, L/2).
    """

    def __init__(self, n_struts, R, w, t, L, M, mag_mode='radial'):
        self.n_struts = n_struts
        self.R = R
        self.w = w    # circumferential width
        self.t = t    # radial thickness
        self.L = L    # axial length
        self.M = M
        self.mag_mode = mag_mode

        self.angles = np.linspace(0, 2 * pi, n_struts, endpoint=False)
        self.cx = R * np.cos(self.angles)
        self.cy = R * np.sin(self.angles)
        self.cz = np.zeros(n_struts)

        # Rotation matrices R such that v_global = R @ v_local
        # Columns are local basis vectors expressed in global frame
        self.rot = []
        for th in self.angles:
            self.rot.append(np.array([
                [ np.cos(th), -np.sin(th), 0.],
                [ np.sin(th),  np.cos(th), 0.],
                [ 0.,          0.,         1.],
            ]))

        # Magnetisation in local frame
        if mag_mode == 'radial':
            self.M_local = np.array([M, 0., 0.])
        elif mag_mode == 'circumferential':
            self.M_local = np.array([0., M, 0.])
        elif mag_mode == 'axial':
            self.M_local = np.array([0., 0., M])
        else:
            raise ValueError(f"Unknown mag_mode: {mag_mode!r}")

    def B_field(self, obs_x, obs_y, obs_z):
        """Return (Bx, By, Bz) at observation points (any broadcastable shape)."""
        shape = np.broadcast(obs_x, obs_y, obs_z).shape
        obs = np.stack([
            np.broadcast_to(obs_x, shape).ravel().astype(float),
            np.broadcast_to(obs_y, shape).ravel().astype(float),
            np.broadcast_to(obs_z, shape).ravel().astype(float),
        ])  # (3, N)

        Bx = np.zeros(obs.shape[1])
        By = np.zeros(obs.shape[1])
        Bz = np.zeros(obs.shape[1])

        a = self.t / 2   # radial half-dim  (local x)
        b = self.w / 2   # circumf. half-dim (local y)
        c = self.L / 2   # axial half-dim    (local z)
        Mlx, Mly, Mlz = self.M_local

        for i in range(self.n_struts):
            Rm = self.rot[i]
            center = np.array([self.cx[i], self.cy[i], self.cz[i]])

            # Transform observations to strut's local frame
            obs_loc = Rm.T @ (obs - center[:, None])  # (3, N)

            bxl, byl, bzl = _B_cuboid_local(
                obs_loc[0], obs_loc[1], obs_loc[2],
                a, b, c, Mlx, Mly, Mlz,
            )

            # Rotate field back to global frame
            B_loc = np.stack([bxl, byl, bzl])   # (3, N)
            B_glob = Rm @ B_loc
            Bx += B_glob[0]
            By += B_glob[1]
            Bz += B_glob[2]

        return Bx.reshape(shape), By.reshape(shape), Bz.reshape(shape)

    def B_magnitude(self, obs_x, obs_y, obs_z):
        Bx, By, Bz = self.B_field(obs_x, obs_y, obs_z)
        return np.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)

    def grad_B(self, obs_x, obs_y, obs_z, dx=5e-7):
        """Full 3D |∇|B|| via central differences."""
        Bp_x = self.B_magnitude(obs_x + dx, obs_y, obs_z)
        Bm_x = self.B_magnitude(obs_x - dx, obs_y, obs_z)
        Bp_y = self.B_magnitude(obs_x, obs_y + dx, obs_z)
        Bm_y = self.B_magnitude(obs_x, obs_y - dx, obs_z)
        Bp_z = self.B_magnitude(obs_x, obs_y, obs_z + dx)
        Bm_z = self.B_magnitude(obs_x, obs_y, obs_z - dx)
        dBdx = (Bp_x - Bm_x) / (2 * dx)
        dBdy = (Bp_y - Bm_y) / (2 * dx)
        dBdz = (Bp_z - Bm_z) / (2 * dx)
        return np.sqrt(dBdx ** 2 + dBdy ** 2 + dBdz ** 2)


def make_ring_3d(**overrides):
    """Create a StentRing3D with defaults + overrides."""
    p = {**DEFAULTS, **overrides}
    return StentRing3D(p['n_struts'], p['R'], p['w'], p['t'], p['L'],
                       p['M'], p.get('mag_mode', 'radial'))


# ===========================================================
# Figure generation
# ===========================================================

def fig1_single_strut():
    """Single strut: B and gradient vs distance, charge vs dipole."""
    print("  Fig 1: Single strut radial profile...")
    
    p = DEFAULTS
    # Single strut at origin, magnetised in +x
    ring1 = StentRing2D(1, 0, p['w'], p['t'], p['M'], 'radial')
    # Override: place strut at origin with M in +x
    ring1.strut_cx = np.array([0.0])
    ring1.strut_cy = np.array([0.0])
    ring1.Mx = np.array([p['M']])
    ring1.My = np.array([0.0])
    ring1.angles = np.array([0.0])
    
    # Also dipole version
    ring1d = StentRing2D(1, 0, p['w'], p['t'], p['M'], 'radial', use_dipole=True)
    ring1d.strut_cx = np.array([0.0])
    ring1d.strut_cy = np.array([0.0])
    ring1d.Mx = np.array([p['M']])
    ring1d.My = np.array([0.0])
    ring1d.angles = np.array([0.0])
    
    d = np.linspace(10e-6, 2e-3, 300)
    surface = p['t'] / 2
    obs_x = d + surface
    obs_y = np.zeros_like(obs_x)
    
    B_charge = ring1.B_magnitude(obs_x, obs_y)
    G_charge = ring1.grad_B(obs_x, obs_y)
    B_dipole = ring1d.B_magnitude(obs_x, obs_y)
    G_dipole = ring1d.grad_B(obs_x, obs_y)
    
    d_um = d * 1e6
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    ax = axes[0]
    ax.semilogy(d_um, B_charge*1000, 'b-', lw=2, label='Surface charge model')
    ax.semilogy(d_um, B_dipole*1000, 'r--', lw=1.5, alpha=0.8, label='Dipole approx.')
    ax.set_xlabel('Distance from strut surface (μm)')
    ax.set_ylabel('|B| (mT)')
    ax.set_title('(a) Magnetic flux density')
    ax.legend()
    ax.set_xlim(0, 2000)
    
    ax = axes[1]
    ax.semilogy(d_um, G_charge, 'b-', lw=2, label='Surface charge model')
    ax.semilogy(d_um, G_dipole, 'r--', lw=1.5, alpha=0.8, label='Dipole approx.')
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=':', lw=1.5, alpha=0.7, label=lbl)
    ax.set_xlabel('Distance from strut surface (μm)')
    ax.set_ylabel('|∇|B|| (T/m)')
    ax.set_title('(b) Field gradient')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 2000); ax.set_ylim(0.1, None)
    
    fig.suptitle(f'Single Magnetised Strut (w={p["w"]*1e6:.0f} μm, '
                 f't={p["t"]*1e6:.0f} μm, M={p["M"]/1e6:.1f} MA/m)',
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig1_single_strut.png')
    fig.savefig(OUT / 'fig1_single_strut.pdf')
    plt.close()


def fig2_ring_heatmaps():
    """Stent ring cross-section: B and gradient heatmaps."""
    print("  Fig 2: Ring cross-section heatmaps...")
    
    ring = make_ring()
    ext = 3.5e-3
    n = 200
    x = np.linspace(-ext, ext, n)
    y = np.linspace(-ext, ext, n)
    X, Y = np.meshgrid(x, y)
    
    Bmag = ring.B_magnitude(X, Y)
    Gmag = ring.grad_B(X, Y)
    
    # Mask inside struts
    for i in range(ring.n_struts):
        dist = np.sqrt((X - ring.strut_cx[i])**2 + (Y - ring.strut_cy[i])**2)
        mask = dist < max(ring.w, ring.t) * 1.2
        Bmag[mask] = np.nan
        Gmag[mask] = np.nan
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax = axes[0]
    im = ax.pcolormesh(X*1e3, Y*1e3, Bmag*1000,
                        norm=LogNorm(vmin=0.01, vmax=100),
                        cmap='inferno', shading='auto')
    fig.colorbar(im, ax=ax, label='|B| (mT)', shrink=0.85)
    circ = Circle((0,0), DEFAULTS['R']*1e3, fill=False, color='w', ls='--', lw=1)
    ax.add_patch(circ)
    for i in range(ring.n_struts):
        ax.plot(ring.strut_cx[i]*1e3, ring.strut_cy[i]*1e3, 'ws', ms=4)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.set_title('(a) Magnetic flux density |B|')
    ax.set_aspect('equal')
    
    ax = axes[1]
    im = ax.pcolormesh(X*1e3, Y*1e3, Gmag,
                        norm=LogNorm(vmin=1, vmax=5000),
                        cmap='inferno', shading='auto')
    fig.colorbar(im, ax=ax, label='|∇|B|| (T/m)', shrink=0.85)
    circ2 = Circle((0,0), DEFAULTS['R']*1e3, fill=False, color='w', ls='--', lw=1)
    ax.add_patch(circ2)
    for i in range(ring.n_struts):
        ax.plot(ring.strut_cx[i]*1e3, ring.strut_cy[i]*1e3, 'ws', ms=4)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.set_title('(b) Field gradient |∇|B||')
    ax.set_aspect('equal')
    
    fig.suptitle(f'Stent Ring Cross-Section ({DEFAULTS["n_struts"]} struts, '
                 f'R={DEFAULTS["R"]*1e3:.1f} mm, M={DEFAULTS["M"]/1e6:.1f} MA/m)',
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig2_ring_heatmaps.png')
    fig.savefig(OUT / 'fig2_ring_heatmaps.pdf')
    plt.close()


def fig3_gradient_vs_distance():
    """Gradient vs distance with thresholds — through strut and between struts."""
    print("  Fig 3: Gradient vs distance...")
    
    ring = make_ring()
    R = DEFAULTS['R']; t = DEFAULTS['t']
    r_outer = R + t/2
    
    d = np.linspace(5e-6, 2e-3, 300)
    
    # Through a strut (angle=0)
    obs_x_through = d + r_outer
    obs_y_through = np.zeros_like(d)
    G_through = ring.grad_B(obs_x_through, obs_y_through)
    
    # Between struts
    angle_bet = pi / DEFAULTS['n_struts']
    obs_x_bet = (d + r_outer) * np.cos(angle_bet)
    obs_y_bet = (d + r_outer) * np.sin(angle_bet)
    G_between = ring.grad_B(obs_x_bet, obs_y_bet)
    
    d_um = d * 1e6
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    for ax, yscale in zip(axes, ['linear', 'log']):
        plot_fn = ax.plot if yscale == 'linear' else ax.semilogy
        plot_fn(d_um, G_through, 'b-', lw=2, label='Through strut')
        plot_fn(d_um, G_between, 'b--', lw=1.5, alpha=0.7, label='Between struts')
        for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
            ax.axhline(val, color=c, ls=':', lw=1.5, alpha=0.7, label=lbl)
        ax.set_xlabel('Distance from stent surface (μm)')
        ax.set_ylabel('|∇|B|| (T/m)')
        ax.legend(fontsize=8)
        ax.set_xlim(0, 2000)
    
    axes[0].set_title('(a) Linear scale')
    axes[0].set_ylim(0, min(3000, np.nanmax(G_through)*1.1))
    axes[1].set_title('(b) Log scale')
    axes[1].set_ylim(0.1, 1e4)
    
    fig.suptitle(f'Gradient vs Distance from Stent Surface '
                 f'({DEFAULTS["n_struts"]} struts, M={DEFAULTS["M"]/1e6:.1f} MA/m)',
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig3_gradient_vs_distance.png')
    fig.savefig(OUT / 'fig3_gradient_vs_distance.pdf')
    plt.close()


def fig4_magnetisation_sweep():
    """Effect of magnetisation on gradient and capture distance."""
    print("  Fig 4: Magnetisation sweep...")
    
    M_vals = [0.2e6, 0.5e6, 0.8e6, 1.0e6, 1.2e6]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(M_vals)))
    
    R = DEFAULTS['R']; t = DEFAULTS['t']
    r_outer = R + t/2
    d = np.linspace(5e-6, 1.5e-3, 250)
    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    d_um = d * 1e6
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    cap_dists = {lbl: [] for lbl in THRESHOLDS}
    
    for M, col in zip(M_vals, colors):
        ring = make_ring(M=M)
        G = ring.grad_B(obs_x, obs_y)
        axes[0].semilogy(d_um, G, color=col, lw=1.8, label=f'M = {M/1e6:.1f} MA/m')
        
        for lbl, thr in THRESHOLDS.items():
            above = np.where(G >= thr)[0]
            cap_dists[lbl].append(d_um[above[-1]] if len(above) > 0 else 0)
    
    ax = axes[0]
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=':', lw=1.2, alpha=0.5)
    ax.set_xlabel('Distance from stent surface (μm)')
    ax.set_ylabel('|∇|B|| (T/m)')
    ax.set_title('(a) Gradient profiles')
    ax.legend(fontsize=8); ax.set_xlim(0, 1500); ax.set_ylim(0.1, 1e4)
    
    ax = axes[1]
    M_arr = np.array(M_vals) / 1e6
    for (lbl, dists), c in zip(cap_dists.items(), TH_COLORS):
        ax.plot(M_arr, dists, 'o-', color=c, lw=2, ms=8, label=lbl)
    ax.set_xlabel('Magnetisation M (MA/m)')
    ax.set_ylabel('Capture distance (μm)')
    ax.set_title('(b) Capture distance vs magnetisation')
    ax.legend(fontsize=8)
    
    fig.suptitle('Effect of Magnetisation on Cell Capture Range', fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig4_magnetisation_sweep.png')
    fig.savefig(OUT / 'fig4_magnetisation_sweep.pdf')
    plt.close()


def fig5_strut_dimensions():
    """Effect of strut thickness and width on gradient."""
    print("  Fig 5: Strut dimensions sweep...")
    
    R = DEFAULTS['R']
    
    thicknesses = [40e-6, 60e-6, 80e-6, 100e-6, 120e-6]
    widths = [50e-6, 80e-6, 100e-6, 150e-6, 200e-6]
    colors_t = plt.cm.plasma(np.linspace(0.15, 0.85, len(thicknesses)))
    colors_w = plt.cm.cividis(np.linspace(0.15, 0.85, len(widths)))
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    d = np.linspace(5e-6, 1e-3, 200)
    
    ax = axes[0]
    for thick, col in zip(thicknesses, colors_t):
        r_outer = R + thick/2
        obs_x = d + r_outer
        obs_y = np.zeros_like(d)
        ring = make_ring(t=thick)
        G = ring.grad_B(obs_x, obs_y)
        ax.semilogy(d*1e6, G, color=col, lw=1.8, label=f't = {thick*1e6:.0f} μm')
    
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=':', lw=1.2, alpha=0.5)
    ax.set_xlabel('Distance from strut surface (μm)')
    ax.set_ylabel('|∇|B|| (T/m)')
    ax.set_title(f'(a) Varying thickness (w={DEFAULTS["w"]*1e6:.0f} μm)')
    ax.legend(fontsize=8); ax.set_xlim(0, 1000); ax.set_ylim(0.1, 1e4)
    
    ax = axes[1]
    r_outer = R + DEFAULTS['t']/2
    for w, col in zip(widths, colors_w):
        obs_x = d + r_outer
        obs_y = np.zeros_like(d)
        ring = make_ring(w=w)
        G = ring.grad_B(obs_x, obs_y)
        ax.semilogy(d*1e6, G, color=col, lw=1.8, label=f'w = {w*1e6:.0f} μm')
    
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=':', lw=1.2, alpha=0.5)
    ax.set_xlabel('Distance from strut surface (μm)')
    ax.set_ylabel('|∇|B|| (T/m)')
    ax.set_title(f'(b) Varying width (t={DEFAULTS["t"]*1e6:.0f} μm)')
    ax.legend(fontsize=8); ax.set_xlim(0, 1000); ax.set_ylim(0.1, 1e4)
    
    fig.suptitle('Effect of Strut Dimensions on Field Gradient', fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig5_strut_dimensions.png')
    fig.savefig(OUT / 'fig5_strut_dimensions.pdf')
    plt.close()


def fig6_n_struts():
    """Effect of number of struts on gradient and angular uniformity."""
    print("  Fig 6: Number of struts study...")
    
    n_vals = [4, 6, 8, 10, 12]
    colors = plt.cm.tab10(np.linspace(0, 0.5, len(n_vals)))
    
    R = DEFAULTS['R']; t = DEFAULTS['t']
    r_outer = R + t/2
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    # (a) Gradient through strut
    d = np.linspace(5e-6, 1e-3, 200)
    ax = axes[0]
    for ns, col in zip(n_vals, colors):
        ring = make_ring(n_struts=ns)
        obs_x = d + r_outer
        obs_y = np.zeros_like(d)
        G = ring.grad_B(obs_x, obs_y)
        ax.semilogy(d*1e6, G, color=col, lw=1.8, label=f'{ns} struts')
    
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=':', lw=1.2, alpha=0.5)
    ax.set_xlabel('Distance from stent surface (μm)')
    ax.set_ylabel('|∇|B|| (T/m)')
    ax.set_title('(a) Gradient through strut')
    ax.legend(fontsize=8); ax.set_xlim(0, 1000); ax.set_ylim(0.1, 1e4)
    
    # (b) Angular variation at fixed distance
    fixed_d = 200e-6
    r_eval = r_outer + fixed_d
    n_angles = 200
    angles = np.linspace(0, 2*pi, n_angles, endpoint=False)
    
    ax = axes[1]
    for ns, col in zip(n_vals, colors):
        ring = make_ring(n_struts=ns)
        obs_x = r_eval * np.cos(angles)
        obs_y = r_eval * np.sin(angles)
        G = ring.grad_B(obs_x, obs_y)
        ax.plot(np.degrees(angles), G, color=col, lw=1.5, label=f'{ns} struts')
    
    ax.set_xlabel('Angle (degrees)')
    ax.set_ylabel('|∇|B|| (T/m)')
    ax.set_title(f'(b) Angular variation at {fixed_d*1e6:.0f} μm from surface')
    ax.legend(fontsize=8); ax.set_xlim(0, 360)
    
    fig.suptitle('Effect of Number of Struts on Gradient Distribution', fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig6_n_struts.png')
    fig.savefig(OUT / 'fig6_n_struts.pdf')
    plt.close()


def fig7_gradient_contours():
    """Cross-section with gradient contours at capture thresholds."""
    print("  Fig 7: Gradient heatmap with contours...")
    
    ring = make_ring()
    ext = 3.0e-3
    n = 200
    x = np.linspace(-ext, ext, n)
    y = np.linspace(-ext, ext, n)
    X, Y = np.meshgrid(x, y)
    
    Gmag = ring.grad_B(X, Y)
    
    for i in range(ring.n_struts):
        dist = np.sqrt((X - ring.strut_cx[i])**2 + (Y - ring.strut_cy[i])**2)
        Gmag[dist < max(ring.w, ring.t) * 1.2] = np.nan
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    
    im = ax.pcolormesh(X*1e3, Y*1e3, Gmag,
                        norm=LogNorm(vmin=5, vmax=5000),
                        cmap='inferno', shading='auto')
    
    # Contours at thresholds
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        # Need non-NaN data for contours
        Gmag_filled = np.nan_to_num(Gmag, nan=0)
        cs = ax.contour(X*1e3, Y*1e3, Gmag_filled, levels=[val],
                        colors=[c], linewidths=2)
        if cs.allsegs[0]:
            ax.clabel(cs, fmt=f'{val} T/m', fontsize=9, colors=[c])
    
    circ = Circle((0,0), DEFAULTS['R']*1e3, fill=False, color='w', ls='--', lw=1.5)
    ax.add_patch(circ)
    for i in range(ring.n_struts):
        ax.plot(ring.strut_cx[i]*1e3, ring.strut_cy[i]*1e3, 'ws', ms=5)
    
    fig.colorbar(im, ax=ax, label='|∇|B|| (T/m)', shrink=0.85)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.set_title(f'Gradient Cross-Section with Capture Thresholds\n'
                 f'({DEFAULTS["n_struts"]} struts, R={DEFAULTS["R"]*1e3:.1f} mm, '
                 f'M={DEFAULTS["M"]/1e6:.1f} MA/m)')
    ax.set_aspect('equal')
    plt.tight_layout()
    fig.savefig(OUT / 'fig7_gradient_contours.png')
    fig.savefig(OUT / 'fig7_gradient_contours.pdf')
    plt.close()


def fig8_force_parameter():
    """B·∇|B| force parameter — proportional to magnetic capture force."""
    print("  Fig 8: Force parameter...")
    
    ring = make_ring()
    R = DEFAULTS['R']; t = DEFAULTS['t']
    r_outer = R + t/2
    
    d = np.linspace(5e-6, 1.5e-3, 250)
    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    
    Bmag = ring.B_magnitude(obs_x, obs_y)
    Gmag = ring.grad_B(obs_x, obs_y)
    F_param = Bmag * Gmag  # T²/m
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    ax = axes[0]
    ax.semilogy(d*1e6, F_param, 'b-', lw=2)
    ax.set_xlabel('Distance from stent surface (μm)')
    ax.set_ylabel('B·∇|B| (T²/m)')
    ax.set_title('(a) Force parameter vs distance')
    
    ax = axes[1]
    l1, = ax.semilogy(d*1e6, Bmag*1000, 'b-', lw=2, label='|B| (mT)')
    ax2 = ax.twinx()
    l2, = ax2.semilogy(d*1e6, Gmag, 'r-', lw=2, label='|∇|B|| (T/m)')
    ax.set_xlabel('Distance from stent surface (μm)')
    ax.set_ylabel('|B| (mT)', color='b')
    ax2.set_ylabel('|∇|B|| (T/m)', color='r')
    ax.set_title('(b) B and gradient decomposition')
    ax.legend([l1, l2], [l1.get_label(), l2.get_label()], fontsize=9)
    
    fig.suptitle('Magnetic Force Parameter for Cell Capture', fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig8_force_parameter.png')
    fig.savefig(OUT / 'fig8_force_parameter.pdf')
    plt.close()


# ===========================================================
# 3D figure generation
# ===========================================================

def fig9_axial_profile():
    """
    Axial (z) variation of gradient at the stent wall for several strut lengths.
    Only possible with the 3D model — shows finite-length end effects.
    """
    print("  Fig 9: Axial gradient profile (3D)...")

    R = DEFAULTS['R']
    t = DEFAULTS['t']
    r_obs = R + t / 2 + 200e-6   # 200 µm outside stent wall, through a strut (angle=0)

    L_vals = [200e-6, 500e-6, 1e-3, 2e-3, 5e-3]
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(L_vals)))

    z = np.linspace(-4e-3, 4e-3, 300)
    obs_x = np.full_like(z, r_obs)
    obs_y = np.zeros_like(z)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for L, col in zip(L_vals, colors):
        ring = make_ring_3d(L=L)
        G = ring.grad_B(obs_x, obs_y, z)
        lbl = f'L = {L*1e3:.1f} mm'
        axes[0].semilogy(z * 1e3, G, color=col, lw=1.8, label=lbl)
        axes[1].plot(z * 1e3, G, color=col, lw=1.8, label=lbl)

    for ax in axes:
        for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
            ax.axhline(val, color=c, ls=':', lw=1.2, alpha=0.6, label=lbl)
        ax.axvline(0, color='gray', ls='--', lw=1, alpha=0.5)
        ax.set_xlabel('Axial position z (mm)')
        ax.set_ylabel('|∇|B|| (T/m)')
        ax.set_xlim(-4, 4)
        ax.legend(fontsize=8)

    axes[0].set_title('(a) Log scale')
    axes[0].set_ylim(0.1, 1e4)
    axes[1].set_title('(b) Linear scale')
    axes[1].set_ylim(0, 1200)

    fig.suptitle(f'Axial Gradient Profile — 3D Akoun & Yonnet Model\n'
                 f'(obs. point: r = R + t/2 + 200 µm, through strut)',
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig9_axial_profile.png')
    fig.savefig(OUT / 'fig9_axial_profile.pdf')
    plt.close()


def fig10_2d_vs_3d():
    """
    Gradient vs radial distance: 2D surface-charge model vs 3D Akoun & Yonnet
    at the midplane (z=0) for several strut lengths.
    Long struts → 3D converges to 2D; short struts → 3D is weaker.
    """
    print("  Fig 10: 2D vs 3D comparison...")

    R = DEFAULTS['R']
    t = DEFAULTS['t']
    d = np.linspace(5e-6, 1.5e-3, 250)
    r_outer = R + t / 2
    obs_x = d + r_outer
    obs_y = np.zeros_like(d)
    d_um = d * 1e6

    # 2D reference
    ring_2d = make_ring()
    G_2d = ring_2d.grad_B(obs_x, obs_y)

    L_vals = [200e-6, 500e-6, 1e-3, 5e-3]
    colors_3d = plt.cm.viridis(np.linspace(0.1, 0.9, len(L_vals)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    for ax, yscale in zip(axes, ['log', 'linear']):
        fn = ax.semilogy if yscale == 'log' else ax.plot
        fn(d_um, G_2d, 'k-', lw=2.5, label='2D surface-charge (∞ length)', zorder=5)
        for L, col in zip(L_vals, colors_3d):
            ring3 = make_ring_3d(L=L)
            G_3d = ring3.grad_B(obs_x, obs_y, np.zeros_like(obs_x))
            fn(d_um, G_3d, '--', color=col, lw=1.8,
               label=f'3D Akoun-Yonnet  L={L*1e3:.1f} mm')
        for (_, val), c in zip(THRESHOLDS.items(), TH_COLORS):
            ax.axhline(val, color=c, ls=':', lw=1.2, alpha=0.5)
        ax.set_xlabel('Distance from stent surface (µm)')
        ax.set_ylabel('|∇|B|| (T/m)')
        ax.set_xlim(0, 1500)
        ax.legend(fontsize=8)

    axes[0].set_title('(a) Log scale')
    axes[0].set_ylim(0.1, 1e4)
    axes[1].set_title('(b) Linear scale')
    axes[1].set_ylim(0, 2500)

    fig.suptitle('2D Surface-Charge vs 3D Akoun & Yonnet — Midplane (z = 0)',
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / 'fig10_2d_vs_3d.png')
    fig.savefig(OUT / 'fig10_2d_vs_3d.pdf')
    plt.close()


def fig11_rz_heatmap():
    """
    Gradient in the r-z (axial cross-section) plane — only possible with 3D model.
    Reveals the axial extent of the capture zone and end-effects at strut tips.
    """
    print("  Fig 11: r-z heatmap (3D)...")

    ring = make_ring_3d()
    R = DEFAULTS['R']
    t = DEFAULTS['t']
    L = DEFAULTS['L']

    # r from just inside wall to 3 mm outside; z from -3 mm to +3 mm
    r_vals = np.linspace(R - t, R + 3e-3, 120)
    z_vals = np.linspace(-3e-3, 3e-3, 200)
    R_grid, Z_grid = np.meshgrid(r_vals, z_vals)

    # Evaluate along the θ=0 direction (through a strut)
    obs_x = R_grid
    obs_y = np.zeros_like(R_grid)
    obs_z = Z_grid

    Gmag = ring.grad_B(obs_x, obs_y, obs_z)

    # Mask strut interior (radial band, axial extent)
    inside = (np.abs(R_grid - R) < t / 2) & (np.abs(Z_grid) < L / 2)
    Gmag[inside] = np.nan

    fig, ax = plt.subplots(figsize=(10, 7))

    im = ax.pcolormesh((R_grid - R) * 1e3, Z_grid * 1e3, Gmag,
                       norm=LogNorm(vmin=5, vmax=5000),
                       cmap='inferno', shading='auto')
    fig.colorbar(im, ax=ax, label='|∇|B|| (T/m)', shrink=0.85)

    # Threshold contours
    Gmag_filled = np.nan_to_num(Gmag, nan=0)
    for (_, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        cs = ax.contour((R_grid - R) * 1e3, Z_grid * 1e3, Gmag_filled,
                        levels=[val], colors=[c], linewidths=1.8)
        if cs.allsegs[0]:
            ax.clabel(cs, fmt=f'{val} T/m', fontsize=9, colors=[c])

    # Mark strut outline
    from matplotlib.patches import Rectangle
    strut = Rectangle((-t / 2 * 1e3, -L / 2 * 1e3),
                      t * 1e3, L * 1e3,
                      linewidth=1.2, edgecolor='white', facecolor='none',
                      linestyle='--')
    ax.add_patch(strut)

    ax.set_xlabel('Radial distance from stent wall (mm)')
    ax.set_ylabel('Axial position z (mm)')
    ax.set_title(f'Gradient in r-z Plane — 3D Akoun & Yonnet\n'
                 f'({DEFAULTS["n_struts"]} struts, L={L*1e3:.1f} mm, '
                 f'M={DEFAULTS["M"]/1e6:.1f} MA/m, through-strut slice)')
    plt.tight_layout()
    fig.savefig(OUT / 'fig11_rz_heatmap.png')
    fig.savefig(OUT / 'fig11_rz_heatmap.pdf')
    plt.close()


# ===========================================================
# Main
# ===========================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Stent Magnetic Field Gradient — Figure Generation")
    print("=" * 60)
    
    figs_2d = [
        fig1_single_strut,
        fig2_ring_heatmaps,
        fig3_gradient_vs_distance,
        fig4_magnetisation_sweep,
        fig5_strut_dimensions,
        fig6_n_struts,
        fig7_gradient_contours,
        fig8_force_parameter,
    ]

    figs_3d = [
        fig9_axial_profile,
        fig10_2d_vs_3d,
        fig11_rz_heatmap,
    ]

    print("\n--- 2D surface-charge model ---")
    for i, fn in enumerate(figs_2d, 1):
        fn()
        print(f"  ✓ Figure {i} complete")

    print("\n--- 3D Akoun & Yonnet model ---")
    for i, fn in enumerate(figs_3d, len(figs_2d) + 1):
        fn()
        print(f"  ✓ Figure {i} complete")

    print(f"\nAll figures saved to: {OUT}")
    print("=" * 60)
