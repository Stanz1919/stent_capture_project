"""
figures.style
=============
Dissertation-quality matplotlib rcParams applied as a context manager or
directly at import time.

Usage (at import):
    import stent_capture.figures.style  # applies rcParams immediately

Usage (selective):
    from stent_capture.figures.style import apply_style
    apply_style()
"""

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — safe for headless/batch use
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


def apply_style() -> None:
    """Apply dissertation rcParams to the current matplotlib session."""
    plt.rcParams.update(_RCPARAMS)


# Apply automatically when this module is imported.
apply_style()
