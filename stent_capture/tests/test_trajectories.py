"""
tests.test_trajectories
=======================
Five tests for simulation.trajectories.integrate_trajectory.

Tests cover:
1. Pure advection in zero magnetic field → straight-line trajectory.
2. Axial velocity matches Poiseuille prediction at intermediate points.
3. Capture near the stent inner surface with 200 pg SPION loading.
4. Escape at high flow velocity with 10 pg loading (centreline injection).
5. Timeout safety: max_time terminates the integration cleanly.

Run with::

    python -m pytest stent_capture/tests/test_trajectories.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from stent_capture.core.field_model import StentRing
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.simulation.trajectories import integrate_trajectory

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

_STENT = dict(n_struts=8, R=1.5e-3, w=100e-6, t=80e-6, L=500e-6)
_R_VES = 1.54e-3


def _make_ring(M: float = 1.0e6) -> StentRing:
    return StentRing(
        _STENT["n_struts"], _STENT["R"], _STENT["w"],
        _STENT["t"], _STENT["L"], M,
        mag_mode="radial", assume_saturation=True,
    )


def _make_flow(v_mean: float = 0.2) -> BloodFlow:
    return BloodFlow(vessel_radius=_R_VES, mean_velocity=v_mean)


# ===========================================================================
# Test 1 — Straight-line advection in zero magnetic field
# ===========================================================================

class TestStraightLineFlow:
    """
    With M = 0 and no external field, F_mag = 0 everywhere.
    The cell velocity equals the local Poiseuille blood velocity (axial only).
    A cell injected at a fixed radial position should remain at that position
    and escape at z = z_end.
    """

    def test_straight_line_no_lateral_drift(self):
        ring  = _make_ring(M=0.0)
        tf    = TotalField(ring, None)   # no external field
        flow  = _make_flow(0.2)
        cell  = SPIONLabelledCell()

        r_inject = np.array([0.5e-3, 0.0, -2e-3])
        traj = integrate_trajectory(
            cell, tf, flow, ring, r_inject,
            z_end=2e-3, max_time=0.5,
        )

        assert traj.status == 'escaped', (
            f"Expected 'escaped', got {traj.status!r}"
        )

        final = traj.positions[-1]

        # x must not drift (tolerance: 100 nm)
        assert abs(final[0] - 0.5e-3) < 1e-7, (
            f"x drift = {abs(final[0] - 0.5e-3)*1e9:.0f} nm (limit 100 nm)"
        )
        # y must not drift (tolerance: 100 nm)
        assert abs(final[1]) < 1e-7, (
            f"y drift = {abs(final[1])*1e9:.0f} nm (limit 100 nm)"
        )
        # z must stop at z_end (tolerance: 1 µm)
        assert abs(final[2] - 2e-3) < 1e-6, (
            f"|z_final - z_end| = {abs(final[2] - 2e-3)*1e6:.2f} µm (limit 1 µm)"
        )


# ===========================================================================
# Test 2 — Axial velocity matches Poiseuille at intermediate points
# ===========================================================================

class TestPoiseuillVelocity:
    """
    In zero magnetic field the ODE reduces to dr/dt = v_blood(r).
    The cell's axial velocity at any stored position must equal the
    Poiseuille prediction for that radial coordinate.
    """

    def test_velocity_matches_poiseuille(self):
        ring  = _make_ring(M=0.0)
        tf    = TotalField(ring, None)
        flow  = _make_flow(0.2)
        cell  = SPIONLabelledCell()

        r_inject = np.array([0.5e-3, 0.0, -2e-3])
        # max_step=1e-3 forces ≥ 10 output points over the ~11 ms transit;
        # without it RK45 takes 2-3 large steps on this smooth ODE.
        traj = integrate_trajectory(
            cell, tf, flow, ring, r_inject,
            z_end=2e-3, max_time=0.5, max_step=1e-3,
        )

        n = len(traj.positions)
        assert n >= 10, f"Too few trajectory points ({n}); need ≥ 10"

        # 10 evenly-spaced interior points (skip first and last)
        idx = np.linspace(1, n - 2, 10, dtype=int)
        pos_sample = traj.positions[idx]

        v_computed  = traj.velocities[idx, 2]          # axial component (m/s)
        v_poiseuille = flow.velocity_at(pos_sample)[:, 2]

        np.testing.assert_allclose(
            v_computed, v_poiseuille, rtol=1e-3,
            err_msg=(
                "Axial velocity deviates from Poiseuille by > 0.1% at "
                "one or more intermediate trajectory points."
            ),
        )


# ===========================================================================
# Test 3 — Capture near the stent inner surface
# ===========================================================================

class TestCaptureNearStrut:
    """
    A cell injected 5 µm inside the lumen from the stent inner surface
    (r = R − t/2 − 5 µm) with 200 pg SPION loading and B0 = 0.5 T should
    be captured quickly: the magnetic drift velocity at this distance far
    exceeds the flow drag, and the cell travels only 5 µm radially to reach
    the wall event boundary at r = R − t/2.
    """

    def test_captured_with_200pg(self):
        ring  = _make_ring()
        tf    = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        flow  = _make_flow(0.05)
        cell  = SPIONLabelledCell(spion_mass_per_cell=200e-15)

        r_inner  = _STENT["R"] - _STENT["t"] / 2   # 1.46 mm
        r_inject = np.array([r_inner - 5e-6, 0.0, 0.0])

        traj = integrate_trajectory(
            cell, tf, flow, ring, r_inject,
            z_end=2e-3, max_time=0.1,
        )

        assert traj.status == 'captured', (
            f"Expected 'captured', got {traj.status!r}"
        )

        # Capture position must be close to the stent inner surface (< 20 µm)
        cap = traj.capture_position
        r_cap = np.sqrt(cap[0]**2 + cap[1]**2)
        assert abs(r_cap - r_inner) < 20e-6, (
            f"Capture radius {r_cap*1e6:.1f} µm deviates from "
            f"inner surface {r_inner*1e6:.1f} µm by more than 20 µm"
        )


# ===========================================================================
# Test 4 — Escape at high flow velocity
# ===========================================================================

class TestEscapeHighFlow:
    """
    A cell at the vessel centreline with 10 pg loading (Polyak default) under
    v_mean = 0.5 m/s should escape.  At the centreline:
      - axial velocity = 2 × v_mean = 1.0 m/s
      - net radial force ≈ 0 (8-strut ring is symmetric about the axis)
    The cell travels from z = −2 mm to z = 2 mm in ~4 ms.
    """

    def test_centreline_escape(self):
        ring  = _make_ring()
        tf    = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        flow  = _make_flow(0.5)
        cell  = SPIONLabelledCell()   # 10 pg

        r_inject = np.array([0.0, 0.0, -2e-3])
        traj = integrate_trajectory(
            cell, tf, flow, ring, r_inject,
            z_end=2e-3, max_time=0.1,
        )

        assert traj.status == 'escaped', (
            f"Expected 'escaped', got {traj.status!r}"
        )
        # Final z should be within 1 µm of z_end
        assert abs(traj.positions[-1, 2] - 2e-3) < 1e-6, (
            f"|z_final − z_end| = {abs(traj.positions[-1,2] - 2e-3)*1e6:.2f} µm"
        )


# ===========================================================================
# Test 5 — Timeout safety (no infinite loop)
# ===========================================================================

class TestTimeoutSafety:
    """
    With very slow flow (v_mean = 0.001 m/s) the transit time exceeds 1 s.
    Setting max_time = 0.01 s must terminate cleanly with status = 'error'
    and t_final ≤ max_time + one adaptive step.  This guards against the
    integration hanging on degenerate flow conditions.
    """

    def test_max_time_terminates(self):
        ring  = _make_ring()
        tf    = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        flow  = _make_flow(0.001)
        cell  = SPIONLabelledCell()

        # Transit time ≈ 2e-3 / (2 × 0.001) = 1 s >> max_time
        r_inject = np.array([0.0, 0.0, -1e-3])
        traj = integrate_trajectory(
            cell, tf, flow, ring, r_inject,
            z_end=2e-3, max_time=0.01,
        )

        assert traj.status == 'error', (
            f"Expected 'error' (timeout), got {traj.status!r}"
        )
        # Allow one extra RK45 step beyond max_time (solver may overshoot slightly)
        assert traj.times[-1] <= 0.012, (
            f"Integration overran max_time: t_final = {traj.times[-1]*1e3:.1f} ms"
        )
