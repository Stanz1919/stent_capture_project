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


# ===========================================================================
# Test 6: Maximum force ratio in lumen < 1.0 at physiological flow (Fig 16)
# ===========================================================================

class TestNoStaticCapturePhysiological:
    """
    Documents the key finding visualised in Fig 16: at cerebral arterial
    flow (v_mean = 0.05 m/s, the slowest case), the static Furlani & Ng force
    balance criterion is NEVER satisfied anywhere in the lumen for a 10 pg
    SPION-loaded cell with B0 = 0.5 T axial.

    The maximum force ratio |F_mag|/|F_drag| is between 0.01 and 1.0, meaning
    drag always exceeds magnetic force by at least one order of magnitude in
    the lumen interior (away from the stent surface).

    This test prevents silent regression if parameters (SPION mass, flow,
    B0) are changed later, and motivates Stage 3 trajectory integration.
    """

    def test_max_force_ratio_below_unity_at_slow_flow(self):
        """Max |F_mag|/|F_drag| in lumen < 1.0 at v_mean=0.05 m/s, B0=0.5T."""
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()
        flow = _default_flow(mean_velocity=0.05)

        # Sample a radial line inside the lumen (from near stent surface to centre)
        R = _STENT["R"]
        t = _STENT["t"]
        r_lumen = R - t / 2        # stent inner surface
        r_vals = np.linspace(10e-6, r_lumen - 10e-6, 50)
        pts = np.column_stack([r_vals, np.zeros_like(r_vals), np.zeros_like(r_vals)])

        result = capture_map(cell, tf, flow, pts)
        ratio = result["F_mag"] / np.where(result["F_drag"] > 0, result["F_drag"], np.inf)
        max_ratio = float(ratio.max())

        assert max_ratio < 1.0, (
            f"Static force balance should predict NO capture in lumen at "
            f"v_mean=0.05 m/s; expected max ratio < 1.0, got {max_ratio:.4f}. "
            "If this fails, parameters changed significantly — update Fig 16 caption."
        )
        assert max_ratio > 0.01, (
            f"Max ratio unexpectedly low ({max_ratio:.6f}); check SPION parameters."
        )


# ===========================================================================
# Test 7: Capture distance is monotonically non-decreasing with SPION loading
# ===========================================================================

class TestLoadingSweep:
    """
    Documents the SPION loading sweep analysed in Fig 17.
    Uses direction='inward' so the sweep covers the actual lumen region
    (r < R - t/2) where cells flow.
    """

    def test_capture_distance_monotonic_in_loading(self):
        """Capture distance increases monotonically with SPION loading."""
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        loadings = [10e-15, 30e-15, 100e-15, 300e-15]   # kg
        flow = _default_flow(mean_velocity=0.05)

        distances = []
        for m in loadings:
            cell = SPIONLabelledCell(spion_mass_per_cell=m)
            d = capture_distance(cell, tf, flow, direction="inward")
            distances.append(d)

        for i in range(len(distances) - 1):
            assert distances[i + 1] >= distances[i], (
                f"Capture distance should be non-decreasing with loading: "
                f"{loadings[i]*1e15:.0f} pg → {distances[i]*1e6:.1f} µm, "
                f"{loadings[i+1]*1e15:.0f} pg → {distances[i+1]*1e6:.1f} µm"
            )

    def test_capture_distance_positive_at_100pg(self):
        """
        At 100 pg SPION loading with v_mean=0.05 m/s and B0=0.5T, there should
        be meaningful capture (> 30 µm from stent inner surface). This confirms
        that static capture is physically predicted in the upper experimental
        loading regime and that the code is not accidentally returning zero.

        Threshold lowered from 50 µm → 30 µm when Langevin saturation became
        the default (commit 53f7ff9): at B ≈ 0.5 T the Langevin factor
        3L(ξ)/ξ ≈ 0.70 attenuates χ_eff and shortens the capture zone.
        """
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = SPIONLabelledCell(spion_mass_per_cell=100e-15)
        flow = _default_flow(mean_velocity=0.05)

        d = capture_distance(cell, tf, flow, direction="inward")
        assert d > 30e-6, (
            f"Expected capture distance > 30 µm at 100 pg, v=0.05 m/s, B0=0.5T; "
            f"got {d*1e6:.1f} µm. Check SPION parameters and field model."
        )
