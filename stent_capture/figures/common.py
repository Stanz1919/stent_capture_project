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
# COMSOL calibration — V2-2C (12-cell)
# ---------------------------------------------------------------------------
# The code's radial-M permanent-magnet model approximates COMSOL's soft-
# ferromagnet (μᵣ = 2) behaviour at B0 = 1.5 T (MRI-strength).
#
# At B0 = 1.5 T, COMSOL's ferromagnet is strongly magnetised; the FEM solver
# computes the resulting concentrated-flux gradient profile. The effective
# permanent magnetisation M_eff = 2.20 MA/m reproduces COMSOL's gradient
# profile with <1% error at the 100 T/m and 40 T/m threshold crossings
# (see fig25). This is the sole validated calibration point.

M_COMSOL_EFF_B15: float = 2.20e6   # A/m, calibrated at B0 = 1.5 T (MRI)

# COMSOL V2-2C threshold crossing distances from stent surface (m)
COMSOL_CROSSINGS: dict[str, float] = {
    "300 T/m": 0.120e-3,
    "100 T/m": 0.240e-3,
    "40 T/m":  0.380e-3,
}


# ---------------------------------------------------------------------------
# Factory helpers
# ---------------------------------------------------------------------------

def make_ring(B0_magnitude: float | None = None, **overrides) -> StentRing:
    """Return a StentRing built from DEFAULTS merged with *overrides*.

    Parameters
    ----------
    B0_magnitude : float, optional
        External magnetic field strength in Tesla. If provided and ≈ 1.5 T
        (MRI strength), automatically uses M_COMSOL_EFF_B15 to match COMSOL's
        ferromagnetic response. Otherwise uses DEFAULTS["M"] (1.0 MA/m by default).
        Explicit M in *overrides* always takes precedence.

    **overrides
        Key-value pairs to override DEFAULTS. If "M" is in overrides, it takes
        precedence over the B0_magnitude-based logic.
    """
    p = {**DEFAULTS, **overrides}

    # Adaptive magnetisation: at B0 ≈ 1.5 T, use calibrated M_COMSOL_EFF_B15
    # unless M is explicitly overridden
    if B0_magnitude is not None and "M" not in overrides:
        if abs(B0_magnitude - 1.5) < 0.01:  # within 10 mT of 1.5 T
            p["M"] = M_COMSOL_EFF_B15

    return StentRing(
        p["n_struts"], p["R"], p["w"], p["t"], p["L"],
        p["M"], p.get("mag_mode", "radial"),
    )


def make_ring_comsol_b15(**overrides) -> StentRing:
    """Return a StentRing with M calibrated to match COMSOL V2-2C at B0=1.5 T.

    Uses M_COMSOL_EFF_B15 = 2.20 MA/m so the through-strut gradient profile
    matches COMSOL's threshold crossings with <1% error at 100 T/m and 40 T/m.
    """
    p = {**DEFAULTS, "M": M_COMSOL_EFF_B15, **overrides}
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
