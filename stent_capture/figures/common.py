"""
figures.common
==============
Shared constants, default parameters, and helper factories used by every
figure module.  Import this before any figure-specific code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# Apply dissertation plot style
import stent_capture.figures.style  # noqa: F401 — side-effect import

from stent_capture.core.field_model import StentRing

# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------
# results/ lives at the project root (two levels above this file's package).
OUT = Path(__file__).parents[2] / "results"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Default stent parameters (match original stent_analysis.py exactly)
# ---------------------------------------------------------------------------
DEFAULTS: dict = dict(
    R=1.5e-3,       # ring radius (m)
    w=100e-6,       # strut circumferential width (m)
    t=80e-6,        # strut radial thickness (m)
    L=500e-6,       # strut axial length (m)
    M=1.0e6,        # magnetisation (A/m) — ~304 SS saturation
    n_struts=12,   # V2-2C (12-cell) — headline COMSOL geometry
    mag_mode="radial",
)

# ---------------------------------------------------------------------------
# Capture-threshold styling
# ---------------------------------------------------------------------------
THRESHOLDS: dict[str, float] = {"40 T/m": 40, "100 T/m": 100, "300 T/m": 300}
TH_COLORS: list[str] = ["#2ecc71", "#e74c3c", "#3498db"]

# ---------------------------------------------------------------------------
# COMSOL calibration — V2-2C (12-cell), B0 = 1.5 T MRI-strength
# ---------------------------------------------------------------------------
# The COMSOL model (FEM, mu_r = 2 ferromagnet, B0 = 1.5 T axial) concentrates
# flux axially, so grad_B is dominated by dB_concentrated/dr with no suppression.
# The code's radial-M permanent-magnet model plus separate B0 has the same radial
# gradient character when B0 = 0 (no suppression).  At B0 = 0, calibrating M so
# the code's 100 T/m crossing matches COMSOL's gives M_eff = 0.619 MA/m.
#
# Calibration result (B0 = 0, n = 12, through-strut):
#   300 T/m:  code_cal = 0.149 mm   COMSOL = 0.120 mm  (+24%)
#   100 T/m:  code_cal = 0.240 mm   COMSOL = 0.240 mm  (exact)
#    40 T/m:  code_cal = 0.333 mm   COMSOL = 0.380 mm  (-12%)
M_COMSOL_EFF: float = 0.619e6   # A/m

# COMSOL V2-2C threshold crossing distances from stent surface (m)
COMSOL_CROSSINGS: dict[str, float] = {
    "300 T/m": 0.120e-3,
    "100 T/m": 0.240e-3,
    "40 T/m":  0.380e-3,
}


# ---------------------------------------------------------------------------
# Factory helpers
# ---------------------------------------------------------------------------

def make_ring(**overrides) -> StentRing:
    """Return a StentRing built from DEFAULTS merged with *overrides*."""
    p = {**DEFAULTS, **overrides}
    return StentRing(
        p["n_struts"], p["R"], p["w"], p["t"], p["L"],
        p["M"], p.get("mag_mode", "radial"),
    )


def make_ring_comsol(**overrides) -> StentRing:
    """Return a StentRing with M calibrated to match COMSOL V2-2C crossings (B0=0).

    Uses M_COMSOL_EFF = 0.619 MA/m so the through-strut gradient profile aligns
    with COMSOL's 100 T/m crossing at 0.240 mm from the stent surface.
    """
    p = {**DEFAULTS, "M": M_COMSOL_EFF, **overrides}
    return StentRing(
        p["n_struts"], p["R"], p["w"], p["t"], p["L"],
        p["M"], p.get("mag_mode", "radial"),
    )


def threshold_lines(ax, yscale: str = "log") -> None:
    """Draw horizontal threshold lines on *ax*."""
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)


def save_fig(fig, stem: str) -> None:
    """Save figure to OUT as both PNG and PDF."""
    fig.savefig(OUT / f"{stem}.png")
    fig.savefig(OUT / f"{stem}.pdf")
