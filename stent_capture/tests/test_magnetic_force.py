"""
tests.test_magnetic_force
=========================
Tests for SPIONLabelledCell and magnetic_force().

Run with::

    python -m pytest stent_capture/tests/test_magnetic_force.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from stent_capture.core.field_model import StentRing
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force, MU_0
from numpy import pi

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_DEFAULT_STENT = dict(n_struts=8, R=1.5e-3, w=100e-6, t=80e-6, L=500e-6, M=1.0e6)


def _make_ring(**kw) -> StentRing:
    p = {**_DEFAULT_STENT, **kw}
    return StentRing(p["n_struts"], p["R"], p["w"], p["t"], p["L"],
                     p["M"], mag_mode="radial", assume_saturation=True)


def _default_cell() -> SPIONLabelledCell:
    return SPIONLabelledCell()


# ===========================================================================
# Test 1: Force is zero when B is spatially uniform (no gradient)
# ===========================================================================

class _UniformField:
    """Minimal mock of a spatially-uniform B field."""
    def __init__(self, B0_vec):
        self.B0 = np.asarray(B0_vec, dtype=float)

    def field_at(self, points):
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[np.newaxis, :]
        return np.broadcast_to(self.B0, (len(pts), 3)).copy()


class TestZeroGradientField:
    """Uniform field has zero gradient, so force must be zero everywhere."""

    @pytest.mark.parametrize("B0_vec", [
        [0.5, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.3, 0.2, 0.1],
    ])
    def test_force_zero_in_uniform_field(self, B0_vec):
        cell = _default_cell()
        field = _UniformField(B0_vec)
        pts = np.array([[1e-3, 0.0, 0.0], [0.0, 1e-3, 0.0], [0.0, 0.0, 1e-3]])
        F = magnetic_force(cell, field, pts)
        # Force should be < 1 fN (FD noise only)
        assert np.all(np.abs(F) < 1e-15), (
            f"Force in uniform field should be ~0, got max |F| = {np.abs(F).max():.3e} N"
        )


# ===========================================================================
# Test 2: Force points toward strut (radially inward) outside the stent
# ===========================================================================

class TestForceDirection:
    """
    At a point radially outside the stent outer surface (along +x), the
    gradient of |B| points toward the strut (−x direction), so the radial
    component of the magnetic force must be negative.
    """

    def test_radial_force_inward(self):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()

        R = _DEFAULT_STENT["R"]
        t = _DEFAULT_STENT["t"]
        r_outer = R + t / 2
        pt = np.array([[r_outer + 200e-6, 0.0, 0.0]])

        F = magnetic_force(cell, tf, pt)
        Fx = float(F[0, 0])
        assert Fx < 0.0, (
            f"Radial force at r = r_outer + 200um should be negative (inward), "
            f"got Fx = {Fx:.3e} N"
        )

    def test_radial_force_inward_no_B0(self):
        ring = _make_ring()
        tf = TotalField(ring)
        cell = _default_cell()

        R = _DEFAULT_STENT["R"]
        t = _DEFAULT_STENT["t"]
        pt = np.array([[R + t / 2 + 100e-6, 0.0, 0.0]])

        F = magnetic_force(cell, tf, pt)
        assert float(F[0, 0]) < 0.0, "Radial force without B0 should also be inward"


# ===========================================================================
# Test 3: Force scales linearly with chi_spion
# ===========================================================================

class TestChiScaling:
    @pytest.mark.parametrize("chi_ratio", [0.5, 2.0, 4.0])
    def test_force_scales_linearly_with_chi(self, chi_ratio):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))

        cell_ref = SPIONLabelledCell(spion_susceptibility=1.0)
        cell_scaled = SPIONLabelledCell(spion_susceptibility=1.0 * chi_ratio)

        pt = np.array([[_DEFAULT_STENT["R"] + _DEFAULT_STENT["t"] / 2 + 200e-6,
                        0.0, 0.0]])

        F_ref    = magnetic_force(cell_ref,    tf, pt)
        F_scaled = magnetic_force(cell_scaled, tf, pt)

        np.testing.assert_allclose(
            np.linalg.norm(F_scaled),
            np.linalg.norm(F_ref) * chi_ratio,
            rtol=1e-10,
            err_msg=f"Force should scale linearly with chi (ratio {chi_ratio})",
        )


# ===========================================================================
# Test 4: Force scales linearly with V_spion (and hence spion_mass)
# ===========================================================================

class TestVolumeScaling:
    @pytest.mark.parametrize("mass_ratio", [0.1, 2.0, 10.0])
    def test_force_scales_linearly_with_spion_mass(self, mass_ratio):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))

        m_ref = 10e-15   # 10 pg
        cell_ref    = SPIONLabelledCell(spion_mass_per_cell=m_ref)
        cell_scaled = SPIONLabelledCell(spion_mass_per_cell=m_ref * mass_ratio)

        pt = np.array([[_DEFAULT_STENT["R"] + _DEFAULT_STENT["t"] / 2 + 200e-6,
                        0.0, 0.0]])

        F_ref    = magnetic_force(cell_ref,    tf, pt)
        F_scaled = magnetic_force(cell_scaled, tf, pt)

        np.testing.assert_allclose(
            np.linalg.norm(F_scaled),
            np.linalg.norm(F_ref) * mass_ratio,
            rtol=1e-10,
            err_msg=f"Force should scale linearly with spion_mass (ratio {mass_ratio})",
        )


# ===========================================================================
# Test 5: SPIONLabelledCell property consistency
# ===========================================================================

class TestCellProperties:
    def test_volume_spherical(self):
        cell = SPIONLabelledCell(radius=10e-6)
        expected = (4.0 / 3.0) * pi * (10e-6) ** 3
        assert abs(cell.volume - expected) / expected < 1e-12

    def test_spion_volume_from_mass_and_density(self):
        cell = SPIONLabelledCell(spion_mass_per_cell=10e-15, spion_density=5170)
        expected = 10e-15 / 5170
        assert abs(cell.spion_volume - expected) / expected < 1e-12

    def test_spion_volume_fraction_physical(self):
        """SPION volume fraction should be << 1 for a physically realistic cell."""
        cell = _default_cell()
        phi = cell.spion_volume / cell.volume
        assert phi < 0.01, (
            f"SPION volume fraction {phi:.4f} is unrealistically high (>1%); "
            f"check spion_mass_per_cell units (should be kg)"
        )


# ===========================================================================
# Test 6: Order-of-magnitude check vs. literature
# ===========================================================================

class TestOrderOfMagnitude:
    """
    At 100 µm from the stent outer surface with B0 = 0.5 T axial, the force
    on a default SPION-labelled cell should be physically meaningful (between
    1 pN and 10 nN).

    Expected (from force parameter calculations in Stage 1):
        |B| ~ 0.51 T, |nabla|B|| ~ 218 T/m  -->  FP ~ 111 T²/m
        prefactor = V_spion * chi / mu_0
                  = (10e-15/5170) * 2.0 / (4pi * 1e-7)
                  ~ 3.1e-12  [m² A/T]
        |F| ~ 3.1e-12 * 111 ~ 340 pN   (within [1 pN, 10 nN])

    This is consistent with the Furlani & Ng 2006 / Polyak 2008 regime for
    short capture distances (< 200 µm from strut surface) in a strong-gradient
    stent field.  The upper bound of 10 nN rules out arithmetic errors
    (e.g. using ng instead of pg for spion_mass).
    """

    def test_force_in_plausible_range_at_100um(self):
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))
        cell = _default_cell()

        R = _DEFAULT_STENT["R"]
        t = _DEFAULT_STENT["t"]
        pt = np.array([[R + t / 2 + 100e-6, 0.0, 0.0]])

        F = magnetic_force(cell, tf, pt)
        F_mag = float(np.linalg.norm(F))

        assert 1e-12 < F_mag < 1e-8, (
            f"|F| = {F_mag * 1e12:.1f} pN at 100 µm with B0=0.5T. "
            f"Expected 1 pN to 10 nN. "
            f"If outside range, check spion_mass_per_cell is in kg (default 10e-15 kg = 10 pg)."
        )
