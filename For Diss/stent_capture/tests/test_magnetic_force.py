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


_DEFAULT_STENT = dict(n_struts=12, R=1.5e-3, w=100e-6, t=80e-6, L=500e-6, M=1.0e6)


def _make_ring(**kw) -> StentRing:
    p = {**_DEFAULT_STENT, **kw}
    return StentRing(p["n_struts"], p["R"], p["w"], p["t"], p["L"],
                     p["M"], mag_mode="radial", assume_saturation=True)


def _default_cell() -> SPIONLabelledCell:
    return SPIONLabelledCell()



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
        assert np.all(np.abs(F) < 1e-15), (
            f"Force in uniform field should be ~0, got max |F| = {np.abs(F).max():.3e} N"
        )



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



class TestChiScaling:
    @pytest.mark.parametrize("chi_ratio", [0.5, 2.0, 4.0])
    def test_force_scales_linearly_with_chi(self, chi_ratio):
        """Test linear chi scaling in constant-susceptibility (non-saturating) mode."""
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.5]))

        cell_ref = SPIONLabelledCell(spion_susceptibility=1.0, spion_sat_magnetization=None)
        cell_scaled = SPIONLabelledCell(spion_susceptibility=1.0 * chi_ratio, spion_sat_magnetization=None)

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



class TestLangevinSaturation:
    """
    Test the Langevin saturation model for SPION susceptibility.

    At low fields (B→0), chi_eff should match chi_0.
    At high fields (B >> M_sat/chi_0), chi_eff should be much less than chi_0.
    This is the physically realistic model for magnetite SPIONs.
    """

    def test_chi_eff_at_zero_field(self):
        """At B→0, chi_eff should equal chi_0."""
        ring = _make_ring()
        tf = TotalField(ring, UniformExternalField([0.0, 0.0, 0.001]))

        chi_0 = 2.0
        cell = SPIONLabelledCell(spion_susceptibility=chi_0,
                                spion_sat_magnetization=446e3)

        pt = np.array([[_DEFAULT_STENT["R"] + _DEFAULT_STENT["t"] / 2 + 100e-6,
                        0.0, 0.0]])

        from stent_capture.physics.magnetic_force import _chi_effective
        B_mag = 0.001
        chi_eff = _chi_effective(chi_0, 446e3, np.array([B_mag]))[0]
        assert abs(chi_eff - chi_0) / chi_0 < 0.01, (
            f"chi_eff at B=1mT should be ~chi_0 (within 1%), "
            f"got {chi_eff:.4f} (chi_0={chi_0})"
        )

    def test_chi_eff_at_high_field_is_reduced(self):
        """At B=1.5T, chi_eff should be significantly less than chi_0."""
        chi_0 = 2.0
        M_sat = 446e3
        from stent_capture.physics.magnetic_force import _chi_effective

        chi_eff = _chi_effective(chi_0, M_sat, np.array([1.5]))[0]
        assert chi_eff < chi_0, (
            f"chi_eff should be less than chi_0 at 1.5T, "
            f"got chi_eff={chi_eff:.4f}, chi_0={chi_0}"
        )
        expected = 0.35
        assert abs(chi_eff - expected) / expected < 0.10, (
            f"chi_eff at 1.5T should be ~{expected:.2f}, got {chi_eff:.4f}"
        )

    def test_chi_eff_monotonically_decreasing(self):
        """chi_eff should decrease monotonically with increasing B."""
        chi_0 = 2.0
        M_sat = 446e3
        from stent_capture.physics.magnetic_force import _chi_effective

        B_array = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0])
        chi_eff_array = _chi_effective(chi_0, M_sat, B_array)

        for i in range(len(chi_eff_array) - 1):
            assert chi_eff_array[i] > chi_eff_array[i+1], (
                f"chi_eff should decrease with B: {chi_eff_array}"
            )

    def test_constant_chi_mode_unaffected(self):
        """With spion_sat_magnetization=None, chi should remain constant."""
        chi_0 = 2.0
        from stent_capture.physics.magnetic_force import _chi_effective

        B_array = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0])
        chi_eff_array = _chi_effective(chi_0, None, B_array)

        np.testing.assert_allclose(
            chi_eff_array, chi_0,
            rtol=1e-12,
            err_msg="With M_sat=None, chi_eff should equal chi_0 everywhere"
        )



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
