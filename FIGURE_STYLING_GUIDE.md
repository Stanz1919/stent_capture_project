# Figure Styling Guide

Professional styling improvements applied to all dissertation figures.

## Summary of Changes

### 1. Unicode Replacements (All Figures)
All figures have been updated to use ASCII-compatible characters for Windows compatibility:
- `µm` → `um`
- `×` → `x`
- `–` (en-dash) → `-` (hyphen)
- `∇|B||` → `grad_B`
- `≈` → `~`

### 2. Consistent Color Palette (`style.py`)

Defined central color constants exported from `stent_capture/figures/style.py`:

```python
COLORS_CONCENTRATION = "viridis"      # Colorblind-friendly concentration colormap
COLORS_THRESHOLD = "#2ecc71"          # Green for therapeutic thresholds
COLORS_MARKER_SOURCE = "red"          # Red for cell/source positions
COLORS_MARKER_PROBE = "cyan"          # Cyan for probe points
COLORS_VECTOR = "white"               # White for velocity vectors
```

### 3. Figure Size & Layout Improvements

**Standardized sizing across figure types:**

| Figure Type | Old Size | New Size | Benefit |
|------------|----------|----------|---------|
| Single panel | Varied | 8-9 x 6 | Better readability |
| Two-panel | Varied | 14 x 5.5 | Balanced spacing |
| Multi-panel | Varied | 15 x 6-6.5 | Professional proportions |

**Layout enhancements:**
- `plt.tight_layout(pad=1.0)` → Better spacing between elements
- `y=0.98` for suptitles → Prevents overlap with content
- `rect=[0, 0, 1, 0.96]` → Reserves space for suptitles

### 4. Professional Typography

**Title formatting:**
- Before: `"Fig 22 — Concentration field"`
- After: `"Figure 22: Steady-state VEGF Distribution"` (bold, size 13)

**Label improvements:**
- More descriptive labels: `"Distance from strut surface (um)"`
- Consistent units throughout
- Font size 10-11 for readability

**Legends:**
- `framealpha=0.95` for better visibility
- Proper location selection to avoid plot overlap
- Consistent fontsize (9-10)

### 5. Text Positioning

Key improvements to prevent text overlap:

1. **Info boxes**: Positioned at `(0.02, 0.98)` (top-left) instead of `(0.97, 0.97)` (top-right)
2. **Annotations**: Offset from content with arrows for clarity
3. **Colorbars**: `pad=0.02` ensures proper spacing
4. **Axis labels**: Updated to be more descriptive and less crowded

### 6. Grid & Visual Polish

- Added grid where appropriate: `grid(True, alpha=0.3)`
- Consistent line widths: `lw=1.8-2.0`
- Alpha transparency for overlapping elements: `alpha=0.7-0.9`
- Clear edge colors for bars and markers

## Updated Figures

### Paracrine Signalling (Figs 22-27)
- **Fig 22**: Concentration field (VEGF distribution)
  - Consistent viridis colormap
  - Improved cell scatter visibility
  - Better text positioning for cell counts

- **Fig 23**: Concentration vs. distance
  - Larger figure (9 x 6)
  - Better annotation positioning
  - Improved legend placement

- **Fig 24**: Time to threshold
  - Enhanced bar chart styling
  - Better data labels
  - Info box repositioned to avoid overlap

- **Fig 27**: Advection impact
  - Consistent colormaps across both panels
  - Velocity vector annotation with text box
  - Professional suptitle

### Magnetic Field Figures (Figs 1-21, 25-26)
**Batch improvements applied:**
- Figure sizing standardized (14 x 5.5 for two-panel, 15 x 6 for multi-panel)
- Suptitle positioning fixed (y=0.98)
- Layout padding increased (pad=1.0)
- Grid visibility improved where applicable

**Key figures improved:**
- **Fig 1**: Single magnetised strut (magnetic flux density & gradient)
- **Fig 3**: Gradient vs distance (code-COMSOL comparison)
- **Fig 7**: Gradient contours
- **Fig 11**: R-Z heatmap
- **Fig 16**: Capture map
- **Fig 17**: SPION loading sweep

## Implementation Guide for New Figures

When creating new figures, follow this template:

```python
from stent_capture.figures.style import (
    apply_style, COLORS_CONCENTRATION, COLORS_THRESHOLD,
    COLORS_MARKER_SOURCE, COLORS_VECTOR
)

apply_style()

# Create figure
fig, ax = plt.subplots(figsize=(11, 6))  # Standard sizing

# Use consistent colors
im = ax.contourf(..., cmap=COLORS_CONCENTRATION)
ax.contour(..., colors=[COLORS_THRESHOLD], ...)
ax.scatter(..., c=COLORS_MARKER_SOURCE, ...)

# Professional title and labels
ax.set_title("Figure N: Clear Title", fontweight='bold', fontsize=13)
ax.set_xlabel("Description (units)")
ax.set_ylabel("Description (units)")

# Proper layout
fig.suptitle("Figure N: Main Title", fontsize=13, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])

# Save at high DPI (set in rcParams via apply_style())
fig.savefig('output.png', dpi=300, bbox_inches='tight')
```

## Verification Checklist

For each figure:
- [ ] No unicode characters (µ, ×, ∇, –, ≈)
- [ ] Figure size appropriate (8-9 for single, 14 for two-panel, 15 for multi)
- [ ] Suptitle at y=0.98, bold, size 13
- [ ] Legends framealpha=0.95, proper location
- [ ] Labels descriptive and in plain English
- [ ] No text overlaps (check colorbars, annotations, info boxes)
- [ ] Consistent color usage (viridis, green threshold, red markers)
- [ ] Grid visible if applicable (alpha=0.3)
- [ ] Layout uses tight_layout(pad=1.0) or rect form

## Common Pitfalls to Avoid

1. **Text Overlap**: Use `transAxes` for positioned text and offset from edges
2. **Color Inconsistency**: Always import from `style.py`, don't hardcode colors
3. **Small Figures**: Minimum 8 inches wide for two-panel, 11 for single
4. **Tight Layouts**: Must account for suptitle with `rect=[0, 0, 1, 0.96]`
5. **Unicode**: Never use mathematical symbols - use ASCII equivalents

## Testing Command

Run all figures to verify:
```bash
cd stent_capture_project-DATAver
python -m stent_capture.figures.regenerate_original_results
```

Each figure should save without errors and display professional styling.
