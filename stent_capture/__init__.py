"""
stent_capture — magnetic stent field-gradient analysis package.

Stage 1: 3D Akoun & Yonnet field model + uniform external field support.
"""

from .core.field_model import StentRing
from .core.gradient import compute_gradient_magnitude
from .physics.external_field import UniformExternalField, TotalField

__version__ = "0.1.0"
__all__ = ["StentRing", "compute_gradient_magnitude", "UniformExternalField", "TotalField"]
