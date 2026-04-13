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
# ferromagnet (μᵣ = 2) behaviour. Two calibration points are defined:
#
# (1) B0 = 0 (no external field):
#     COMSOL's ferromagnet is unmagnetised at B0 = 0; the code uses an
#     effective permanent M to match COMSOL's gradient profile shape.
#     Calibration: M_eff = 0.619 MA/m matches COMSOL's 100 T/m crossing.
#
# (2) B0 = 1.5 T (MRI-strength, COMSOL's actual operating point):
#     COMSOL's ferromagnet is strongly magnetised; the FEM solver computes
#     the resulting concentrated-flux gradient profile. The code uses an
#     effective permanent M to match COMSOL's threshold crossings exactly.
#     Calibration: M_eff = 2.20 MA/m reproduces COMSOL's gradient profile
#     with <1% error at 100 T/m and 40 T/m crossings.
#
# The ratio 2.20 / 0.619 ≈ 3.55 quantifies the amplification of the
# gradient due to the B0 = 1.5 T external field in COMSOL's model.

M_COMSOL_EFF: float = 0.619e6   # A/m, calibrated at B0 = 0
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


def make_ring_comsol(**overrides) -> StentRing:
    """Return a StentRing with M calibrated to match COMSOL V2-2C at B0=0.

    Uses M_COMSOL_EFF = 0.619 MA/m so the through-strut gradient profile aligns
    with COMSOL's 100 T/m crossing at 0.240 mm from the stent surface.

    Note: For B0 = 1.5 T (MRI), prefer make_ring(B0_magnitude=1.5) which
    automatically uses M_COMSOL_EFF_B15.
    """
    p = {**DEFAULTS, "M": M_COMSOL_EFF, **overrides}
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
