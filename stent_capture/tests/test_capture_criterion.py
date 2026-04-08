"""
tests.test_capture_criterion
============================
Tests for capture_map() and capture_distance().

Run with::

    python -m pytest stent_capture/tests/test_capture_criterion.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from stent_capture.core.field_model import StentRing
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.physics.capture_criterion import capture_map, capture_distance

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

_STENT = dict(n_struts=8, R=1.5e-3, w=100e-6, t=80e-6, L=500e-6, M=1.0e6)


def _make_ring(**kw) -> StentRing:
    p = {**_STENT, **kw}
    return StentRing(p["n_struts"], p["R"], p["w"], p["t"], p["L"],
                     p["M"], mag_mode="radial", assume_saturation=True)


def _default_cell() -> SPIONLabelledCell:
    return SPIONLabelledCell()


def _default_flow(**kw) -> BloodFlow:
    p = dict(vessel_radius=1.54e-3, mean_velocity=0.2)
    p.update(kw)
    return BloodFlow(**p)


# ===========================================================================
# Test 1: capture_map dict structure
# ===========================================================================

class TestCapatureMapStructure:
    def test_keys_present(self):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()
        flow = _default_flow()
        pt = np.array([[_STENT["R"] + _STENT["t"] / 2 + 50e-6, 0.0, 0.0]])
        result = capture_map(cell, tf, flow, pt)
        for key in ["F_mag_vec", "F_drag_vec", "F_mag", "F_drag", "margin", "captured"]:
            assert key in result, f"Key '{key}' missing from capture_map result"

    def test_shapes_consistent(self):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()
        flow = _default_flow()
        N = 10
        pts = np.zeros((N, 3))
        pts[:, 0] = _STENT["R"] + _STENT["t"] / 2 + np.linspace(10e-6, 300e-6, N)
        result = capture_map(cell, tf, flow, pts)
        assert result["F_mag_vec"].shape  == (N, 3)
        assert result["F_drag_vec"].shape == (N, 3)
        assert result["F_mag"].shape      == (N,)
        assert result["F_drag"].shape     == (N,)
        assert result["margin"].shape     == (N,)
        assert result["captured"].shape   == (N,)

    def test_margin_equals_difference(self):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()
        flow = _default_flow()
        pt = np.array([[_STENT["R"] + _STENT["t"] / 2 + 100e-6, 0.0, 0.0]])
        result = capture_map(cell, tf, flow, pt)
        np.testing.assert_allclose(
            result["margin"], result["F_mag"] - result["F_drag"], rtol=1e-12
        )


# ===========================================================================
# Test 2: Capture at stent surface under default B0
# ===========================================================================

class TestCaptureNearSurface:
    """
    Very close to the stent surface (< 10 µm), the gradient is enormous
    (>> 1000 T/m) and F_mag must greatly exceed typical Stokes drag.
    """

    def test_captured_at_5um(self):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()
        flow = _default_flow()

        R = _STENT["R"]
        t = _STENT["t"]
        pt = np.array([[R + t / 2 + 5e-6, 0.0, 0.0]])
        result = capture_map(cell, tf, flow, pt)
        assert result["captured"][0], (
            f"Cell should be captured at 5 µm from strut surface. "
            f"|F_mag|={result['F_mag'][0]*1e12:.1f} pN, "
            f"|F_drag|={result['F_drag'][0]*1e12:.1f} pN"
        )


# ===========================================================================
# Test 3: No capture far from stent (regardless of B0)
# ===========================================================================

class TestNoCaptureFarFromStent:
    """
    Far from the stent (r >> R), |∇B| -> 0 so F_mag -> 0.
    Even modest blood flow will then dominate.
    """

    @pytest.mark.parametrize("B0_z", [0.0, 0.5, 1.0])
    def test_not_captured_at_1mm(self, B0_z):
        ring = _make_ring()
        ext = UniformExternalField([0.0, 0.0, B0_z]) if B0_z > 0 else None
        tf = TotalField(ring, ext)
        cell = _default_cell()
        flow = _default_flow()

        # 1 mm from strut surface, still inside vessel (lumen radius ~1.46mm)
        R = _STENT["R"]
        t = _STENT["t"]
        # Use a point well inside the lumen far from strut
        pt = np.array([[0.5e-3, 0.0, 0.0]])   # centreline, 1 mm from axis
        result = capture_map(cell, tf, flow, pt)
        # |F_drag| at centreline with v_mean=0.2 m/s is ~15 pN (Stokes)
        # |F_mag| at 1 mm from axis (stent is at R=1.5mm) should be tiny
        assert not result["captured"][0], (
            f"Cell at centreline should NOT be captured "
            f"(|F_mag|={result['F_mag'][0]*1e12:.3f} pN, "
            f"|F_drag|={result['F_drag'][0]*1e12:.3f} pN)"
        )


# ===========================================================================
# Test 4: Capture range decreases with increasing flow velocity
# ===========================================================================

class TestCaptureRangeVsVelocity:
    def test_range_decreases_with_velocity(self):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()

        v_vals = [0.05, 0.2, 0.5]
        distances = []
        for v in v_vals:
            flow = _default_flow(mean_velocity=v)
            d = capture_distance(cell, tf, flow)
            distances.append(d)

        for i in range(len(distances) - 1):
            assert distances[i] >= distances[i + 1], (
                f"Capture distance should decrease with velocity: "
                f"v={v_vals[i]} -> d={distances[i]*1e6:.1f}µm, "
                f"v={v_vals[i+1]} -> d={distances[i+1]*1e6:.1f}µm"
            )


# ===========================================================================
# Test 5: Capture range increases with B0 (force parameter enhancement)
# ===========================================================================

class TestCaptureRangeVsB0:
    def test_range_increases_with_B0(self):
        """
        Applying axial B0 increases the force parameter |B|*|nabla B|, which
        should extend the capture distance (or at minimum not reduce it).
        Note: this tests the FORCE PARAMETER, not just |nabla B| — the latter
        decreases (see fig12 physics note) but the former increases (fig13).
        """
        ring = _make_ring()
        cell = _default_cell()
        flow = _default_flow(mean_velocity=0.2)

        d_no_B0 = capture_distance(
            cell, TotalField(ring), flow
        )
        d_B0_05 = capture_distance(
            cell, TotalField(ring, UniformExternalField([0.0, 0.0, 0.5])), flow
        )
        d_B0_10 = capture_distance(
            cell, TotalField(ring, UniformExternalField([0.0, 0.0, 1.0])), flow
        )

        assert d_B0_05 >= d_no_B0, (
            f"B0=0.5T should extend capture distance vs B0=0 "
            f"(got {d_B0_05*1e6:.1f} vs {d_no_B0*1e6:.1f} µm)"
        )
        assert d_B0_10 >= d_B0_05, (
            f"B0=1.0T should extend capture distance vs B0=0.5T "
            f"(got {d_B0_10*1e6:.1f} vs {d_B0_05*1e6:.1f} µm)"
        )
