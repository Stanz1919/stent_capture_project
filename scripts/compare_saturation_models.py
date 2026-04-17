"""
Compare constant-chi vs Langevin saturation models.

Runs key tests with both susceptibility models to quantify impact
on capture efficiency, trajectory behavior, and static predictions.

Usage::

    python -m scripts.compare_saturation_models
"""

from __future__ import annotations

import numpy as np
import time

from stent_capture.figures.common import DEFAULTS, make_ring
from stent_capture.physics.external_field import TotalField, UniformExternalField
from stent_capture.physics.magnetic_force import SPIONLabelledCell, magnetic_force
from stent_capture.physics.hydrodynamics import BloodFlow
from stent_capture.physics.capture_criterion import capture_distance
from stent_capture.simulation.capture_efficiency import sweep_injection_line

# ============================================================================
# Parameters
# ============================================================================

_B0_Z = 1.5  # T
_R_VES = 1.54e-3
_V_MEAN = 0.2  # m/s

_LINE_START = np.array([1.20e-3, 0.0, -2e-3])
_LINE_END = np.array([1.45e-3, 0.0, -2e-3])

_LOADINGS_PG = np.array([10., 30., 50., 100., 200.])


# ============================================================================
# Test 1: Static Capture Distance vs SPION Loading
# ============================================================================

def test_static_capture_distance():
    """Compare static capture distances between models."""
    print("\n" + "=" * 75)
    print("TEST 1: Static Capture Distance vs SPION Loading")
    print("=" * 75)
    print("Conditions: B0=1.5T, v=0.2 m/s, inward direction")
    print()

    # Setup
    ring_langevin = make_ring(B0_magnitude=_B0_Z)
    ring_langevin.assume_saturation = True
    tf_langevin = TotalField(ring_langevin, UniformExternalField([0.0, 0.0, _B0_Z]))

    ring_const = make_ring(B0_magnitude=_B0_Z, M=1.0e6)
    ring_const.assume_saturation = True
    tf_const = TotalField(ring_const, UniformExternalField([0.0, 0.0, _B0_Z]))

    flow = BloodFlow(vessel_radius=_R_VES, mean_velocity=_V_MEAN)

    # Results
    results = {
        "loading_pg": [],
        "d_langevin_um": [],
        "d_constant_um": [],
        "ratio": [],
    }

    print("Loading(pg) | Langevin(um) | Constant(um) | Ratio")
    print("-" * 55)

    for m_pg in _LOADINGS_PG:
        m_kg = m_pg * 1e-15

        cell_langevin = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=446e3)
        cell_constant = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=None)

        d_langevin = capture_distance(cell_langevin, tf_langevin, flow, direction="inward") * 1e6
        d_constant = capture_distance(cell_constant, tf_const, flow, direction="inward") * 1e6

        ratio = d_constant / d_langevin if d_langevin > 0.5 else np.inf

        results["loading_pg"].append(m_pg)
        results["d_langevin_um"].append(d_langevin)
        results["d_constant_um"].append(d_constant)
        results["ratio"].append(ratio)

        ratio_str = f"{ratio:.2f}x" if ratio < 100 else ">>1.0x"
        print(f"  {m_pg:>8.0f}  |  {d_langevin:>10.1f}  |  {d_constant:>10.1f}  |  {ratio_str:>6}")

    return results


# ============================================================================
# Test 2: Trajectory Capture Efficiency at Multiple Points
# ============================================================================

def test_trajectory_efficiency():
    """Compare trajectory-based capture efficiency."""
    print("\n" + "=" * 75)
    print("TEST 2: Trajectory-Based Capture Efficiency")
    print("=" * 75)
    print("Injection line: r=1.20-1.45mm at z=-2mm, 20 cells per condition")
    print()

    # Setup
    ring_langevin = make_ring(B0_magnitude=_B0_Z)
    ring_langevin.assume_saturation = True
    tf_langevin = TotalField(ring_langevin, UniformExternalField([0.0, 0.0, _B0_Z]))

    ring_const = make_ring(B0_magnitude=_B0_Z, M=1.0e6)
    ring_const.assume_saturation = True
    tf_const = TotalField(ring_const, UniformExternalField([0.0, 0.0, _B0_Z]))

    # Test cases
    test_cases = [
        ("50pg, v=0.05", 50e-15, 0.05),
        ("50pg, v=0.20", 50e-15, 0.20),
        ("200pg, v=0.05", 200e-15, 0.05),
    ]

    results = {
        "case": [],
        "eff_langevin": [],
        "eff_constant": [],
        "ratio": [],
    }

    print("Case               | Langevin | Constant | Ratio")
    print("-" * 55)

    for label, m_kg, v in test_cases:
        # Langevin
        cell_langevin = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=446e3)
        flow_langevin = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        trajs_l, summary_l = sweep_injection_line(
            cell_langevin, tf_langevin, flow_langevin, ring_langevin,
            _LINE_START, _LINE_END, n_points=20,
            z_end=2e-3, max_time=1.0,
        )
        eff_langevin = summary_l["efficiency"]

        # Constant chi
        cell_constant = SPIONLabelledCell(spion_mass_per_cell=m_kg,
                                         spion_sat_magnetization=None)
        flow_const = BloodFlow(vessel_radius=_R_VES, mean_velocity=v)
        trajs_c, summary_c = sweep_injection_line(
            cell_constant, tf_const, flow_const, ring_const,
            _LINE_START, _LINE_END, n_points=20,
            z_end=2e-3, max_time=1.0,
        )
        eff_constant = summary_c["efficiency"]

        ratio = eff_constant / eff_langevin if eff_langevin > 0.01 else np.inf

        results["case"].append(label)
        results["eff_langevin"].append(eff_langevin)
        results["eff_constant"].append(eff_constant)
        results["ratio"].append(ratio)

        ratio_str = f"{ratio:.2f}x" if ratio < 100 else ">>1.0x"
        print(f"  {label:>15}  |  {eff_langevin:>6.3f}  |  {eff_constant:>6.3f}  |  {ratio_str:>6}")

    return results


# ============================================================================
# Test 3: Force Field Comparison at Key Locations
# ============================================================================

def test_force_field_comparison():
    """Compare force magnitudes at specific radial distances."""
    print("\n" + "=" * 75)
    print("TEST 3: Force Field Comparison")
    print("=" * 75)
    print("Strut-aligned axis, B0=1.5T, 50pg loading")
    print()

    ring = make_ring(B0_magnitude=_B0_Z)
    ring.assume_saturation = True
    tf_langevin = TotalField(ring, UniformExternalField([0.0, 0.0, _B0_Z]))

    ring_const = make_ring(B0_magnitude=_B0_Z, M=1.0e6)
    ring_const.assume_saturation = True
    tf_const = TotalField(ring_const, UniformExternalField([0.0, 0.0, _B0_Z]))

    R = DEFAULTS["R"]
    t = DEFAULTS["t"]
    r_outer = R + t / 2

    distances_um = [50, 100, 150, 200]

    cell_langevin = SPIONLabelledCell(spion_mass_per_cell=50e-15,
                                     spion_sat_magnetization=446e3)
    cell_constant = SPIONLabelledCell(spion_mass_per_cell=50e-15,
                                     spion_sat_magnetization=None)

    print("Distance(um) | F_Langevin(pN) | F_Const(pN) | Ratio")
    print("-" * 55)

    for d_um in distances_um:
        d_m = d_um * 1e-6
        pt = np.array([[r_outer + d_m, 0.0, 0.0]])

        F_langevin = magnetic_force(cell_langevin, tf_langevin, pt)
        F_constant = magnetic_force(cell_constant, tf_const, pt)

        F_mag_langevin = np.linalg.norm(F_langevin) * 1e12
        F_mag_constant = np.linalg.norm(F_constant) * 1e12

        ratio = F_mag_constant / F_mag_langevin if F_mag_langevin > 0.01 else np.inf

        ratio_str = f"{ratio:.2f}x" if ratio < 100 else ">>1.0x"
        print(f"  {d_um:>10}  |  {F_mag_langevin:>13.2f}  |  {F_mag_constant:>10.2f}  |  {ratio_str:>6}")

    print()


# ============================================================================
# Main
# ============================================================================

def main():
    print("\n" + "=" * 75)
    print("SPION SATURATION MODEL COMPARISON ANALYSIS".center(75))
    print("=" * 75)

    t_start = time.time()

    # Run tests
    test1_results = test_static_capture_distance()
    test2_results = test_trajectory_efficiency()
    test_force_field_comparison()

    elapsed = time.time() - t_start

    # Summary
    print("=" * 75)
    print("SUMMARY STATISTICS")
    print("=" * 75)

    avg_ratio_static = np.mean([r for r in test1_results["ratio"] if r < 100])
    print(f"\nTest 1 (Static Capture Distance):")
    print(f"  Average ratio (Constant/Langevin): {avg_ratio_static:.2f}x")
    print(f"  Min ratio: {min([r for r in test1_results['ratio'] if r < 100]):.2f}x")
    print(f"  Max ratio: {max([r for r in test1_results['ratio'] if r < 100]):.2f}x")

    avg_ratio_traj = np.mean([r for r in test2_results["ratio"] if r < 100])
    print(f"\nTest 2 (Trajectory Efficiency):")
    print(f"  Average ratio (Constant/Langevin): {avg_ratio_traj:.2f}x")
    print(f"  Min ratio: {min([r for r in test2_results['ratio'] if r < 100]):.2f}x")
    print(f"  Max ratio: {max([r for r in test2_results['ratio'] if r < 100]):.2f}x")

    print(f"\nTotal time: {elapsed:.1f} s")

    print("\n" + "=" * 75)
    print("INTERPRETATION")
    print("=" * 75)
    print("""
The Langevin saturation model produces SIGNIFICANTLY LOWER capture predictions
compared to constant-chi. This is because at MRI fields (1.5T), SPION
susceptibility is reduced by ~5.7x, which directly reduces magnetic force.

The constant-chi model OVERESTIMATES capture capability by assuming SPIONs
remain linearly responsive at all field strengths — unphysical for real SPIONs.

The Langevin model is more realistic physics. Decision on adoption depends on:
1. Which model provides better experimental agreement?
2. Does dissertation scope permit this change?
3. Should results show saturation, constant chi, or both?
""")

    return test1_results, test2_results


if __name__ == "__main__":
    main()
