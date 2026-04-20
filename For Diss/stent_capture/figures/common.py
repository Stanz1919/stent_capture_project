"""
figures.common
==============
Shared constants, default parameters, and helper factories used by every
figure module.  Import this before any figure-specific code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import stent_capture.figures.style

from stent_capture.core.field_model import StentRing

OUT = Path(__file__).parents[2] / "results"
OUT.mkdir(exist_ok=True)

DEFAULTS: dict = dict(
    R=1.5e-3,
    w=100e-6,
    t=80e-6,
    L=500e-6,
    M=1.0e6,
    n_struts=12,
    mag_mode="radial",
)

THRESHOLDS: dict[str, float] = {"40 T/m": 40, "100 T/m": 100, "300 T/m": 300}
TH_COLORS: list[str] = ["#2ecc71", "#e74c3c", "#3498db"]

M_COMSOL_EFF_B15: float = 2.20e6

COMSOL_CROSSINGS: dict[str, float] = {
    "300 T/m": 0.120e-3,
    "100 T/m": 0.240e-3,
    "40 T/m":  0.380e-3,
}


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

    if B0_magnitude is not None and "M" not in overrides:
        if abs(B0_magnitude - 1.5) < 0.01:
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
