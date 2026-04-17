# Unified Color Scheme Guide

## Rationale

Dissertation figures previously used scattered colors (red, blue, yellow, green, purple, etc.) making it difficult to match meaning across figures. This unified scheme reduces the palette to **3 primary threshold colors** and **4 marker colors** that are reused consistently.

## Primary Threshold Colors (Green → Orange → Red)

These three colors represent increasing severity/threshold levels across all figures:

### Green: Low/Safe Threshold
```
Color: #2ecc71 (Emerald green)
Usage: 
  - 40 T/m (capture threshold - easiest to achieve)
  - 5 ng/mL VEGF (therapeutic minimum - basal EC secretion)
  - Lower-risk scenarios
```

### Orange: Medium Threshold
```
Color: #f39c12 (Burnt orange)
Usage:
  - 100 T/m (standard clinical capture threshold)
  - 15-25 ng/mL VEGF (therapeutic target range - optimized secretion)
  - Intermediate difficulty/impact
```

### Red: High/Critical Threshold
```
Color: #e74c3c (Bright red)
Usage:
  - 300 T/m (maximum capture threshold - most difficult)
  - >100 ng/mL VEGF (aberrant/excessive growth)
  - High-risk or critical scenarios
```

**Visual progression:** Green → Orange → Red mirrors traffic light semantics (safe → caution → critical)

---

## Marker/Data Colors (Purpose-Based)

### Blue: Code/Simulation Data
```
Color: #3498db (Sky blue)
Usage:
  - Default code model (standard parameters)
  - Simulation results (baseline)
  - Probe/observation points
  - Reusable for plots where code is primary subject
```

### Dark Orange: Calibrated/Adjusted Data
```
Color: #e67e22 (Dark orange)
Usage:
  - Calibrated model (tuned to match reference)
  - Model with adjusted parameters
  - Code results after optimization
  - Distinctly different from medium orange (#f39c12)
```

### Gray: Reference Data (COMSOL, Literature, etc.)
```
Color: #95a5a6 (Medium gray)
Usage:
  - COMSOL FEM results
  - Experimental data from literature
  - External/third-party references
  - Comparison benchmarks
```

### Red: Cell/Source Positions
```
Color: Red (pure)
Usage:
  - Captured cell locations
  - VEGF source positions
  - Any origin/starting point markers
```

### White: Vector/Flow Annotations
```
Color: White
Usage:
  - Velocity vectors
  - Flow arrows
  - Directional annotations (overlaid on colored backgrounds)
```

---

## Implementation

### Importing in Figure Scripts

```python
from stent_capture.figures.style import (
    COLORS_THRESHOLD,          # Green (low)
    COLORS_THRESHOLD_MEDIUM,   # Orange (medium)
    COLORS_THRESHOLD_HIGH,     # Red (high)
    COLORS_CODE_DEFAULT,       # Blue (standard code)
    COLORS_CODE_CALIBRATED,    # Dark orange (adjusted)
    COLORS_MARKER_REFERENCE,   # Gray (reference)
    COLORS_MARKER_SOURCE,      # Red (sources/cells)
    COLORS_VECTOR,             # White (vectors)
    THRESHOLD_COLORS_ORDERED,  # List of [green, orange, red]
)
```

### Example Usage Patterns

#### Pattern 1: Threshold Lines (Ordered Low→High)
```python
import matplotlib.pyplot as plt
from stent_capture.figures.style import THRESHOLD_COLORS_ORDERED

thresholds = [40, 100, 300]  # Ordered from low to high
for thresh, color in zip(thresholds, THRESHOLD_COLORS_ORDERED):
    ax.axhline(thresh, color=color, ls='--', lw=1.5, label=f'{thresh} T/m')
```

#### Pattern 2: Model Comparison
```python
from stent_capture.figures.style import COLORS_CODE_DEFAULT, COLORS_CODE_CALIBRATED

ax.plot(x, y_default, color=COLORS_CODE_DEFAULT, label='Default')
ax.plot(x, y_calibrated, color=COLORS_CODE_CALIBRATED, label='Calibrated')
ax.scatter(x_ref, y_ref, color=COLORS_MARKER_REFERENCE, label='COMSOL')
```

#### Pattern 3: Concentration Fields
```python
from stent_capture.figures.style import COLORS_CONCENTRATION, COLORS_THRESHOLD

im = ax.contourf(..., cmap=COLORS_CONCENTRATION)  # viridis
ax.contour(..., colors=[COLORS_THRESHOLD], ...)   # Green threshold overlay
```

---

## Color Specifications

| Purpose | Hex | RGB | Usage |
|---------|-----|-----|-------|
| Threshold (Low) | #2ecc71 | (46, 204, 113) | 40 T/m, 5 ng/mL |
| Threshold (Medium) | #f39c12 | (243, 156, 18) | 100 T/m, 15-25 ng/mL |
| Threshold (High) | #e74c3c | (231, 76, 60) | 300 T/m, >100 ng/mL |
| Code Default | #3498db | (52, 152, 219) | Standard model |
| Code Calibrated | #e67e22 | (230, 126, 34) | Adjusted model |
| Reference Data | #95a5a6 | (149, 165, 166) | COMSOL, literature |
| Cell Source | Red | (255, 0, 0) | Positions, origins |
| Vector/Flow | White | (255, 255, 255) | Arrows, annotations |

---

## Backwards Compatibility

### Old Approach (Multi-Color)
Figure 3 previously used `TH_COLORS = ["#2ecc71", "#e74c3c", "#3498db"]` which mixed threshold meaning with model types.

### New Approach (Semantic)
- **Thresholds only** use green→orange→red (ascending severity)
- **Models/Data** use blue/orange/gray (type-specific colors)
- **Much fewer total colors** (8 instead of 12+)

### Migration Path
1. Replace `TH_COLORS` usage with `THRESHOLD_COLORS_ORDERED`
2. Replace hardcoded colors ("b", "r", "darkorange") with named constants
3. Verify color usage matches semantic meaning (low/medium/high or default/calibrated/reference)

---

## Accessibility Considerations

This color scheme is designed to be:
- **Colorblind-friendly:** Viridis (concentration) and ordered severity (green→orange→red) both work for most colorblindness types
- **Print-safe:** No pure yellow or light greens; all colors have sufficient contrast
- **High-contrast:** Each color has distinct hue and brightness
- **Semantically meaningful:** Color progression (green safe → red critical) is intuitive

---

## Validation Checklist

For each figure:
- [ ] All threshold lines use THRESHOLD_COLORS_ORDERED (green, orange, red)
- [ ] Model comparisons use COLORS_CODE_DEFAULT (blue) and COLORS_CODE_CALIBRATED (orange)
- [ ] Reference data is COLORS_MARKER_REFERENCE (gray)
- [ ] Cell/source positions are COLORS_MARKER_SOURCE (red)
- [ ] No other colors used except when semantically required (e.g., concentration colormap)
- [ ] Color imports from `style.py`, not hardcoded
- [ ] Document legend clearly explains color meaning

---

## Examples in Current Figures

### ✓ Fig 3 (Gradient vs Distance)
- Blue line: Code default
- Dark orange line: Code calibrated
- Green/orange/red lines: Thresholds (40/100/300 T/m)
- Gray circles: COMSOL reference

### ✓ Fig 22-27 (Paracrine)
- Viridis: VEGF concentration heatmap
- Green contour: Therapeutic threshold (5 ng/mL)
- Red markers: Captured cells

### ✓ Fig 24 (Time to Threshold)
- Green/orange/red bars: Thresholds (ascending distance)
- Yellow background: Info box (light, neutral)

---

## Future Standardization

All remaining figures should migrate to this palette by:
1. Reviewing color usage in each figure
2. Mapping old colors to new semantic colors
3. Updating imports to use `style.py` constants
4. Testing output to confirm visual consistency
