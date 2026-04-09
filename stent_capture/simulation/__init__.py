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

__all__ = ["CellTrajectory", "integrate_trajectory"]
