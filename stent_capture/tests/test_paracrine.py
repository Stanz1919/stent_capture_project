"""
test_paracrine — unit tests for the paracrine signalling module.

Four tests:
1. Zero source → concentration stays zero everywhere.
2. Point source, no decay → 1/r-like fall-off at steady state.
3. Pure decay (D=0 analogue) → exponential decay matches exp(−k t).
4. Therapeutic zone falls within biologically plausible range (100–500 µm)
   for 8 captured cells with literature parameters.
"""

from __future__ import annotations

import numpy as np
import pytest

from stent_capture.paracrine.transport import ParacrineField, D_VEGF_TISSUE, K_DEG_TISSUE
from stent_capture.paracrine.secretion import VEGFSource, Q_CELL_G_PER_S
from stent_capture.paracrine.therapeutic import (
    therapeutic_zone_radius,
    concentration_vs_distance,
    C_THERAPEUTIC_LOW,
)


# ── 1. Zero source → zero concentration ────────────────────────────

def test_zero_source_gives_zero():
    """With no source, steady-state C must be identically zero."""
    pf = ParacrineField(Lx=1e-3, Lz=1e-3, Nx=30, Nz=30)
    C = pf.solve_steady_state()
    assert np.allclose(C, 0.0, atol=1e-30)


# ── 2. Point source, no decay → 1/r decay profile ─────────────────

def test_point_source_no_decay_radial_falloff():
    """
    Single Gaussian source with k=0 (no degradation).  Steady-state
    concentration must decrease monotonically with distance from source.
    """
    Lx = 2e-3
    Lz = 2e-3
    N  = 80
    pf = ParacrineField(Lx=Lx, Lz=Lz, Nx=N, Nz=N,
                        D=D_VEGF_TISSUE, k_deg=1e-20, periodic_x=False)

    src = VEGFSource(sigma=20e-6)
    centre = np.array([[Lx / 2, Lz / 2]])
    pf.source = src.source_field(pf.X, pf.Z, centre)
    C = pf.solve_steady_state()

    r_centres, C_mean = concentration_vs_distance(
        C, pf.X, pf.Z, centre, n_bins=30, r_max=Lx / 3,
    )

    valid = C_mean > 0
    assert np.sum(valid) > 5, "Not enough bins with nonzero concentration"

    diffs = np.diff(C_mean[valid])
    assert np.all(diffs <= 0), (
        "Concentration should decrease monotonically with distance; "
        f"got positive diffs at indices {np.where(diffs > 0)[0]}"
    )


# ── 3. Pure decay → exponential in time ───────────────────────────

def test_pure_decay_exponential():
    """
    Start with uniform C₀ and D ≈ 0.  Only degradation acts:
    C(t) = C₀ exp(−k t).  Check at t = 1/k (≈ 1 e-fold).
    """
    k = K_DEG_TISSUE
    C0_val = 10.0  # ng/mL
    N = 20
    pf = ParacrineField(Lx=1e-3, Lz=1e-3, Nx=N, Nz=N,
                        D=1e-20, k_deg=k, periodic_x=True)

    C0 = np.full((N, N), C0_val)
    t_final = 1.0 / k  # one e-folding time

    dt = 0.1 / k  # keep dt·k = 0.1 for Euler accuracy
    C_snaps, t_snaps = pf.solve_transient(t_final=t_final, n_snapshots=5, dt=dt, C0=C0)

    C_final = C_snaps[-1]
    # Compare against Euler-exact: C = C0 * (1 - k*dt_actual)^n_steps
    n_steps = int(np.round(t_final / dt))
    dt_actual = t_final / n_steps
    expected_euler = C0_val * (1.0 - k * dt_actual) ** n_steps

    np.testing.assert_allclose(
        C_final.mean(), expected_euler,
        rtol=1e-6,
        err_msg=f"Pure decay: expected {expected_euler:.6f}, got {C_final.mean():.6f}",
    )


# ── 4. Therapeutic zone biologically plausible ─────────────────────

def test_therapeutic_zone_plausible():
    """
    Eight captured cells (one per strut) with literature parameters must
    produce a therapeutic zone (C ≥ 5 ng/mL) in the range 0–500 µm.
    Zero is acceptable (source too weak); but if nonzero, must not
    exceed ~1 mm (that would indicate a unit error).
    """
    from stent_capture.figures.common import make_ring, DEFAULTS

    ring = make_ring()
    R = DEFAULTS["R"]
    Lx = 2.0 * np.pi * R
    Lz = 3.0 * DEFAULTS["L"]

    pf = ParacrineField(Lx=Lx, Lz=Lz, Nx=120, Nz=120)

    src = VEGFSource()
    z_c = Lz / 2.0
    cell_pos = src.cell_positions_on_ring(
        ring, z_positions=np.full(ring.n_struts, z_c),
    )
    pf.source = src.source_field(pf.X, pf.Z, cell_pos)
    C = pf.solve_steady_state()

    r_tz = therapeutic_zone_radius(C, pf.X, pf.Z, cell_pos,
                                   threshold=C_THERAPEUTIC_LOW)

    assert r_tz < 1e-3, (
        f"Therapeutic zone radius {r_tz*1e6:.0f} µm exceeds 1 mm — "
        "likely a unit error"
    )
