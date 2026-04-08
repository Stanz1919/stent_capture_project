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
    n_struts=8,
    mag_mode="radial",
)

# ---------------------------------------------------------------------------
# Capture-threshold styling
# ---------------------------------------------------------------------------
THRESHOLDS: dict[str, float] = {"40 T/m": 40, "100 T/m": 100, "300 T/m": 300}
TH_COLORS: list[str] = ["#2ecc71", "#e74c3c", "#3498db"]


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


def threshold_lines(ax, yscale: str = "log") -> None:
    """Draw horizontal threshold lines on *ax*."""
    for (lbl, val), c in zip(THRESHOLDS.items(), TH_COLORS):
        ax.axhline(val, color=c, ls=":", lw=1.5, alpha=0.7, label=lbl)


def save_fig(fig, stem: str) -> None:
    """Save figure to OUT as both PNG and PDF."""
    fig.savefig(OUT / f"{stem}.png")
    fig.savefig(OUT / f"{stem}.pdf")
