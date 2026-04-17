"""
tests.test_hydrodynamics
========================
Tests for BloodFlow and stokes_drag().

Run with::

    python -m pytest stent_capture/tests/test_hydrodynamics.py -v
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy import pi

from stent_capture.physics.hydrodynamics import BloodFlow, stokes_drag
from stent_capture.physics.magnetic_force import SPIONLabelledCell

# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------

_R = 1.54e-3    # vessel radius (m)
_V = 0.2        # mean velocity (m/s)
_ETA = 4e-3     # viscosity (Pa·s)


def _make_flow(**kw) -> BloodFlow:
    p = dict(vessel_radius=_R, mean_velocity=_V, viscosity=_ETA)
    p.update(kw)
    return BloodFlow(**p)


def _make_cell(**kw) -> SPIONLabelledCell:
    return SPIONLabelledCell(**kw)


# ===========================================================================
# Test 1: Velocity profile correctness
# ===========================================================================

class TestPoiseuille:
    """Verify the analytical Poiseuille profile at key locations."""

    def test_centreline_velocity(self):
        """At r=0 the velocity should equal v_max = 2 * v_mean."""
        flow = _make_flow()
        pt = np.array([[0.0, 0.0, 0.0]])
        v = flow.velocity_at(pt)
        assert abs(v[0, 2] - 2 * _V) < 1e-12 * _V, (
            f"Centreline v_z = {v[0,2]:.6f} m/s, expected {2*_V:.6f} m/s"
        )

    def test_wall_velocity_zero(self):
        """At r = R_vessel the velocity must be zero (no-slip)."""
        flow = _make_flow()
        # Slightly inside the wall to stay within the domain
        r = _R * (1 - 1e-9)
        pts = np.array([[r, 0.0, 0.0]])
        v = flow.velocity_at(pts)
        assert abs(v[0, 2]) < 1e-6, (
            f"Near-wall v_z = {v[0,2]:.3e} m/s, expected ~0"
        )

    def test_outside_vessel_zero(self):
        """Points beyond R_vessel should have zero velocity."""
        flow = _make_flow()
        pts = np.array([[_R * 1.5, 0.0, 0.0], [0.0, _R * 2.0, 0.0]])
        v = flow.velocity_at(pts)
        np.testing.assert_array_equal(v, 0.0)

    def test_velocity_only_axial(self):
        """All velocity is along +z; x and y components must be zero."""
        flow = _make_flow()
        rng = np.random.default_rng(42)
        r_vals = rng.uniform(0, _R * 0.99, size=20)
        angles = rng.uniform(0, 2 * pi, size=20)
        pts = np.column_stack([r_vals * np.cos(angles),
                                r_vals * np.sin(angles),
                                np.zeros(20)])
        v = flow.velocity_at(pts)
        np.testing.assert_array_equal(v[:, 0], 0.0)
        np.testing.assert_array_equal(v[:, 1], 0.0)

    def test_mean_velocity_integral(self):
        """
        Cross-sectional average of v_z over many rings should equal v_mean.
        Uses the analytical result: mean = v_max / 2 = v_mean.
        """
        flow = _make_flow()
        # Sample on a grid and area-weight
        n = 200
        r_vals = np.linspace(0, _R, n + 1)[:-1] + _R / (2 * n)
        pts = np.column_stack([r_vals, np.zeros(n), np.zeros(n)])
        v = flow.velocity_at(pts)
        # Area-weighted average: integral(v_z * 2pi r dr) / (pi R²)
        dr = _R / n
        area_weighted = np.sum(v[:, 2] * 2 * pi * r_vals * dr) / (pi * _R ** 2)
        assert abs(area_weighted - _V) / _V < 5e-4, (
            f"Area-weighted mean velocity = {area_weighted:.5f} m/s, "
            f"expected {_V:.5f} m/s"
        )


# ===========================================================================
# Test 2: Shear rate profile
# ===========================================================================

class TestShearRate:
    def test_shear_rate_at_centreline_zero(self):
        """Shear rate at r=0 is zero (velocity maximum)."""
        flow = _make_flow()
        pt = np.array([[0.0, 0.0, 0.0]])
        assert abs(flow.shear_rate_at(pt)[0]) < 1e-10

    def test_shear_rate_increases_with_r(self):
        """Shear rate is monotonically increasing from centre to wall."""
        flow = _make_flow()
        r_vals = np.linspace(0.0, _R * 0.99, 20)
        pts = np.column_stack([r_vals, np.zeros(20), np.zeros(20)])
        gamma = flow.shear_rate_at(pts)
        assert np.all(np.diff(gamma) >= 0), "Shear rate should increase with r"

    def test_wall_shear_stress_formula(self):
        """wall_shear_stress = 4 * eta * v_mean / R_vessel."""
        flow = _make_flow()
        expected = 4 * _ETA * _V / _R
        assert abs(flow.wall_shear_stress - expected) / expected < 1e-12


# ===========================================================================
# Test 3: Stokes drag
# ===========================================================================

class TestStokesDrag:
    def test_drag_zero_at_vessel_wall(self):
        """At r = R_vessel, blood velocity is zero, so drag on stationary cell = 0."""
        flow = _make_flow()
        cell = _make_cell()
        r = _R * (1 - 1e-9)
        pt = np.array([[r, 0.0, 0.0]])
        F = stokes_drag(cell, flow, pt)
        assert np.linalg.norm(F) < 1e-15, (
            f"Drag near wall should be ~0, got {np.linalg.norm(F):.3e} N"
        )

    def test_drag_direction_axial(self):
        """Drag on a stationary cell is along +z (axial flow direction)."""
        flow = _make_flow()
        cell = _make_cell()
        # Point at centreline
        pt = np.array([[0.0, 0.0, 0.0]])
        F = stokes_drag(cell, flow, pt)
        assert abs(F[0, 0]) < 1e-20 and abs(F[0, 1]) < 1e-20, (
            "Drag should be purely axial (Fx=Fy=0 for axial Poiseuille)"
        )
        assert F[0, 2] > 0, "Drag on stationary cell should be along +z"

    def test_drag_scales_linearly_with_radius(self):
        """Stokes formula: F_drag proportional to cell radius."""
        flow = _make_flow()
        pt = np.array([[0.0, 0.0, 0.0]])   # centreline
        cell1 = _make_cell(radius=5e-6)
        cell2 = _make_cell(radius=10e-6)
        F1 = np.linalg.norm(stokes_drag(cell1, flow, pt))
        F2 = np.linalg.norm(stokes_drag(cell2, flow, pt))
        np.testing.assert_allclose(F2 / F1, 2.0, rtol=1e-10)

    def test_drag_scales_linearly_with_velocity(self):
        """Drag proportional to blood velocity (and hence mean_velocity)."""
        pt = np.array([[0.0, 0.0, 0.0]])   # centreline
        cell = _make_cell()
        flow1 = _make_flow(mean_velocity=0.1)
        flow2 = _make_flow(mean_velocity=0.2)
        F1 = np.linalg.norm(stokes_drag(cell, flow1, pt))
        F2 = np.linalg.norm(stokes_drag(cell, flow2, pt))
        np.testing.assert_allclose(F2 / F1, 2.0, rtol=1e-10)

    def test_drag_moving_cell_reduces_force(self):
        """A cell already moving with the flow experiences less drag."""
        flow = _make_flow()
        cell = _make_cell()
        pt = np.array([[0.0, 0.0, 0.0]])   # centreline
        F_stationary = stokes_drag(cell, flow, pt)
        v_with_flow = np.array([[0.0, 0.0, flow.v_max]])   # same as flow
        F_comoving = stokes_drag(cell, flow, pt, cell_velocities=v_with_flow)
        # Co-moving cell should have zero drag
        assert np.linalg.norm(F_comoving) < 1e-20

    def test_stokes_drag_formula_centreline(self):
        """Verify F_drag = 6 pi eta R v_z at centreline."""
        flow = _make_flow()
        cell = _make_cell(radius=10e-6)
        pt = np.array([[0.0, 0.0, 0.0]])
        F = stokes_drag(cell, flow, pt)
        expected = 6 * pi * _ETA * 10e-6 * flow.v_max
        np.testing.assert_allclose(F[0, 2], expected, rtol=1e-12)
