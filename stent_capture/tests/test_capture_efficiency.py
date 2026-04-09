"""
tests.test_capture_efficiency
==============================
Five tests for simulation.capture_efficiency.

Tests
-----
1. Zero-field sweep: all cells escape (efficiency = 0).
2. Summary integrity: keys present, counts sum to n_total.
3. Non-zero capture: 200 pg cells near the lumen wall are captured under slow flow.
4. Velocity monotonicity: efficiency at low v_mean ≥ efficiency at high v_mean.
5. Loading monotonicity: efficiency at high loading ≥ efficiency at low loading.

All tests use ``n_workers=1`` (serial) to avoid subprocess-spawn overhead in CI
and to keep wall-clock time under ~30 s.

Run with::

    python -m pytest stent_capture/tests/test_capture_efficiency.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from stent_capture.core.field_model import StentRing
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.simulation.capture_efficiency import (
    sweep_injection_line,
    capture_efficiency_vs_velocity,
    capture_efficiency_vs_loading,
)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

_R_VES   = 1.54e-3
_R_INNER = 1.5e-3 - 80e-6 / 2   # R − t/2 = 1.46 mm

def _make_ring(M: float = 1.0e6) -> StentRing:
    return StentRing(8, 1.5e-3, 100e-6, 80e-6, 500e-6, M,
                     mag_mode="radial", assume_saturation=True)


def _make_tf(M: float = 1.0e6, B0_z: float = 0.5):
    ring = _make_ring(M)
    return ring, TotalField(ring, UniformExternalField([0.0, 0.0, B0_z]))


def _make_flow(v: float) -> BloodFlow:
    return BloodFlow(vessel_radius=_R_VES, mean_velocity=v)


# Injection line just inside the lumen wall: x ∈ [1.20 mm, 1.45 mm], y = z = 0
_LINE_START = np.array([1.20e-3, 0.0, -2e-3])
_LINE_END   = np.array([1.45e-3, 0.0, -2e-3])

_KW_FAST = dict(z_end=2e-3, max_time=0.30)   # fast kwargs for tests


# ===========================================================================
# Test 1 — Zero magnetic field: all cells must escape
# ===========================================================================

class TestZeroFieldAllEscape:
    """
    With M = 0 (no stent magnetisation) and B0 = 0 there is no magnetic force.
    Every cell follows the Poiseuille streamline and escapes at z = z_end.

    Injection line kept away from the slow near-wall region (r ≤ 1.2 mm) so
    that every cell transits z = −2 mm → +2 mm in < 150 ms < max_time = 0.30 s.
    """

    # Injection line from r = 0.3 mm to r = 1.2 mm — safely away from the
    # slow near-wall zone (v_axial at 1.2 mm ≈ 39 mm/s → transit ~102 ms).
    _START = np.array([0.3e-3, 0.0, -2e-3])
    _END   = np.array([1.2e-3, 0.0, -2e-3])

    def test_efficiency_is_zero(self):
        ring = _make_ring(M=0.0)
        tf   = TotalField(ring, None)    # no external field
        flow = _make_flow(0.05)
        cell = SPIONLabelledCell(spion_mass_per_cell=200e-15)

        trajs, summary = sweep_injection_line(
            cell, tf, flow, ring,
            self._START, self._END,
            n_points=5,
            n_workers=1,
            **_KW_FAST,
        )

        assert summary['efficiency'] == 0.0, (
            f"Expected zero capture in zero field; "
            f"got {summary['n_captured']}/{summary['n_total']}"
        )
        assert all(t.status == 'escaped' for t in trajs), (
            f"All trajectories should be 'escaped' in zero field; "
            f"got statuses {[t.status for t in trajs]}"
        )


# ===========================================================================
# Test 2 — Summary dictionary integrity
# ===========================================================================

class TestSummaryIntegrity:
    """
    The summary dict returned by sweep_injection_line must contain all required
    keys, and the counts must be consistent: n_captured + n_escaped + n_error
    == n_total.
    """

    def test_summary_keys_and_counts(self):
        ring = _make_ring(M=0.0)
        tf   = TotalField(ring, None)
        flow = _make_flow(0.2)
        cell = SPIONLabelledCell()

        trajs, summary = sweep_injection_line(
            cell, tf, flow, ring,
            _LINE_START, _LINE_END,
            n_points=6,
            n_workers=1,
            **_KW_FAST,
        )

        required_keys = {
            'n_total', 'n_captured', 'n_escaped', 'n_error',
            'efficiency', 'injection_points',
        }
        assert required_keys.issubset(summary.keys()), (
            f"Missing summary keys: {required_keys - summary.keys()}"
        )

        # Count consistency
        assert (summary['n_captured'] + summary['n_escaped'] + summary['n_error']
                == summary['n_total']), "Counts do not sum to n_total."

        # efficiency is the captured fraction
        if summary['n_total'] > 0:
            expected_eff = summary['n_captured'] / summary['n_total']
            assert abs(summary['efficiency'] - expected_eff) < 1e-12

        # injection_points shape
        assert summary['injection_points'].shape == (6, 3), (
            f"Expected injection_points shape (6, 3), "
            f"got {summary['injection_points'].shape}"
        )

        # Trajectory list length matches n_total
        assert len(trajs) == summary['n_total']


# ===========================================================================
# Test 3 — Non-zero capture near the stent wall
# ===========================================================================

class TestNonZeroCapture:
    """
    With 200 pg SPION loading, B0 = 0.5 T, and slow flow (v_mean = 0.05 m/s),
    cells injected close to the lumen wall (r ≈ 1.40–1.45 mm) must be captured.
    Capture efficiency must be > 0.
    """

    def test_some_cells_captured(self):
        ring, tf = _make_tf()
        flow = _make_flow(0.05)
        cell = SPIONLabelledCell(spion_mass_per_cell=200e-15)

        trajs, summary = sweep_injection_line(
            cell, tf, flow, ring,
            _LINE_START, _LINE_END,
            n_points=5,
            n_workers=1,
            **_KW_FAST,
        )

        assert summary['n_captured'] > 0, (
            f"Expected at least one capture near the lumen wall with 200 pg "
            f"SPION loading and v_mean = 0.05 m/s; got 0/{summary['n_total']}."
        )


# ===========================================================================
# Test 4 — Velocity monotonicity
# ===========================================================================

class TestVelocityMonotonicity:
    """
    At fixed SPION loading (200 pg) capture efficiency must decrease (or stay
    equal) as mean blood velocity increases.  Specifically, efficiency at
    v_mean = 0.02 m/s must be ≥ efficiency at v_mean = 0.50 m/s.
    """

    def test_efficiency_decreases_with_velocity(self):
        ring, tf = _make_tf()
        cell = SPIONLabelledCell(spion_mass_per_cell=200e-15)

        results = capture_efficiency_vs_velocity(
            cell, tf, ring,
            velocities=[0.02, 0.50],
            line_start=_LINE_START,
            line_end=_LINE_END,
            n_cells_per_velocity=5,
            vessel_radius=_R_VES,
            n_workers=1,
            **_KW_FAST,
        )

        eff_low, eff_high = results['efficiencies']
        assert eff_low >= eff_high, (
            f"Expected efficiency at v=0.02 ({eff_low:.3f}) ≥ "
            f"efficiency at v=0.50 ({eff_high:.3f})."
        )


# ===========================================================================
# Test 5 — Loading monotonicity
# ===========================================================================

class TestLoadingMonotonicity:
    """
    At fixed flow velocity (0.05 m/s) capture efficiency must increase (or stay
    equal) with SPION loading.  Specifically, efficiency at 200 pg must be ≥
    efficiency at 10 pg.
    """

    def test_efficiency_increases_with_loading(self):
        ring, tf = _make_tf()
        flow = _make_flow(0.05)

        def cell_factory(m_kg: float) -> SPIONLabelledCell:
            return SPIONLabelledCell(spion_mass_per_cell=m_kg)

        results = capture_efficiency_vs_loading(
            cell_factory, tf, flow, ring,
            loadings_kg=[10e-15, 200e-15],
            line_start=_LINE_START,
            line_end=_LINE_END,
            n_cells_per_loading=5,
            n_workers=1,
            **_KW_FAST,
        )

        eff_low, eff_high = results['efficiencies']
        assert eff_high >= eff_low, (
            f"Expected efficiency at 200 pg ({eff_high:.3f}) ≥ "
            f"efficiency at 10 pg ({eff_low:.3f})."
        )
