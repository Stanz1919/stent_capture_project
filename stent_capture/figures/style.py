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

# Consistent color palette for all figures
COLORS_CONCENTRATION = "viridis"  # Primary colormap for concentration fields

# Unified threshold colors (magnetic field & therapeutic)
COLORS_THRESHOLD = "#2ecc71"  # Green for therapeutic thresholds (5 ng/mL)
COLORS_THRESHOLD_MEDIUM = "#f39c12"  # Orange for medium thresholds (100 T/m, 15-25 ng/mL)
COLORS_THRESHOLD_HIGH = "#e74c3c"  # Red for high thresholds (300 T/m)

# Markers and annotations
COLORS_MARKER_SOURCE = "red"  # Red for cell/source positions
COLORS_MARKER_PROBE = "#3498db"  # Blue for probe/observation points
COLORS_MARKER_REFERENCE = "#95a5a6"  # Gray for reference data (COMSOL, etc.)
COLORS_VECTOR = "white"  # White for velocity vectors/arrows

# Code comparison colors (consistent across figures)
COLORS_CODE_DEFAULT = "#3498db"  # Blue for default/standard model
COLORS_CODE_CALIBRATED = "#e67e22"  # Dark orange for calibrated/adjusted model

# Discrete threshold palette (replaces multicolor approach)
# Use these consistently across all threshold indicators
THRESHOLD_COLORS_ORDERED = [
    COLORS_THRESHOLD,  # Green - lowest/safe threshold
    COLORS_THRESHOLD_MEDIUM,  # Orange - intermediate
    COLORS_THRESHOLD_HIGH,  # Red - highest/critical
]


def apply_style() -> None:
    """Apply dissertation rcParams to the current matplotlib session."""
    plt.rcParams.update(_RCPARAMS)


# Apply automatically when this module is imported.
apply_style()
