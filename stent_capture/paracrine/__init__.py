"""
paracrine
=========
VEGF paracrine signalling from captured MSCs on a magnetised
stent.  Solves 2-D reaction–diffusion on the vessel-wall tissue plane and
computes therapeutic-zone metrics.

Modules
-------
transport   — ``ParacrineField`` finite-difference PDE solver.
secretion   — ``VEGFSource`` maps captured-cell positions to a source field.
therapeutic — threshold analysis and time-to-therapeutic helpers.
"""

from stent_capture.paracrine.transport import ParacrineField
from stent_capture.paracrine.secretion import VEGFSource
from stent_capture.paracrine.therapeutic import (
    therapeutic_zone_radius,
    time_to_threshold,
)

__all__ = [
    "ParacrineField",
    "VEGFSource",
    "therapeutic_zone_radius",
    "time_to_threshold",
]
