"""
stent_capture.simulation
========================
Cell trajectory integration sub-package (Stage 3).

Public API
----------
CellTrajectory
    Result object from a single integrated trajectory.
integrate_trajectory
    Integrate one cell path using scipy.solve_ivp (Stage 3a).
"""

from stent_capture.simulation.trajectories import CellTrajectory, integrate_trajectory
from stent_capture.simulation.capture_efficiency import (
    sweep_injection_line,
    capture_efficiency_vs_velocity,
    capture_efficiency_vs_loading,
)

__all__ = [
    "CellTrajectory",
    "integrate_trajectory",
    "sweep_injection_line",
    "capture_efficiency_vs_velocity",
    "capture_efficiency_vs_loading",
]
