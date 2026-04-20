"""
figures.style
=============
Dissertation-quality matplotlib rcParams and consistent styling utilities.

Usage (at import):
    import stent_capture.figures.style  # applies rcParams immediately

Usage (selective):
    from stent_capture.figures.style import apply_style, COLORS_CONCENTRATION
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_RCPARAMS = {
    "font.family":       "serif",
    "font.size":         11,
    "axes.labelsize":    13,
    "axes.titlesize":    14,
    "legend.fontsize":   9,
    "xtick.labelsize":   11,
    "ytick.labelsize":   11,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "axes.grid":         True,
    "grid.alpha":        0.3,
    "lines.linewidth":   1.8,
}

COLORS_CONCENTRATION = "viridis"

COLORS_THRESHOLD = "#2ecc71"
COLORS_THRESHOLD_MEDIUM = "#f39c12"
COLORS_THRESHOLD_HIGH = "#e74c3c"

COLORS_MARKER_SOURCE = "red"
COLORS_MARKER_PROBE = "#3498db"
COLORS_MARKER_REFERENCE = "#95a5a6"
COLORS_VECTOR = "white"

COLORS_CODE_DEFAULT = "#3498db"
COLORS_CODE_CALIBRATED = "#e67e22"

THRESHOLD_COLORS_ORDERED = [
    COLORS_THRESHOLD,
    COLORS_THRESHOLD_MEDIUM,
    COLORS_THRESHOLD_HIGH,
]


def apply_style() -> None:
    plt.rcParams.update(_RCPARAMS)


apply_style()
