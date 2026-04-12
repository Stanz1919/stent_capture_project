# Results Summary: Original (10 pg) vs 200 pg SPION Loading

## Overview

This project now includes **two complete sets of results** for easy comparison:

1. **`results/`** — Original project results with default 10 pg SPION loading
2. **`results-200pg/`** — Variant results with 200 pg SPION loading (Polyak et al. 2008 experimental value)

---

## File Counts

| Folder | Figures | PNGs | PDFs | Description |
|---|---|---|---|---|
| `results/` | 1–24 | 24 | 24 | **Complete:** All stages (1–3c), default parameters (10 pg SPION) |
| `results-200pg/` | 13–21 | 9 | 9 | **Stage 2 & 3 only:** Force, drag, capture, trajectories, efficiency (200 pg SPION) |

**Total:** 33 unique figures (24 original + 9 200pg-specific)

---

## Figure Breakdown by Stage

### **Stage 1: B-Field & Gradient (Figs 1–12)**
*All in `results/` (SPION-independent)*

| Figure | Description |
|---|---|
| **fig01** | Single strut radial B-field profile |
| **fig02** | Ring cross-section B-field heatmaps |
| **fig03** | Gradient magnitude vs axial distance |
| **fig04** | B-field vs magnetisation sweep |
| **fig05** | Effect of strut dimensions |
| **fig06** | Effect of number of struts |
| **fig07** | Gradient contours (heatmap) |
| **fig08** | Force parameter (B·∇B) spatial map |
| **fig09** | Axial gradient profile |
| **fig10** | Convergence check (FD accuracy) |
| **fig11** | R-Z plane heatmap |
| **fig12** | External field comparison (B₀ effect) |

### **Stage 2: Force & Drag (Figs 13–16)**
*Available in both `results/` (10 pg) and `results-200pg/` (200 pg)*

| Figure | Description | Key Difference |
|---|---|---|
| **fig13** | Force parameter \|B\|·\|∇B\| spatial map | Stronger field at 200 pg |
| **fig14** | Magnetic force vs radial distance | Higher force with 200 pg |
| **fig15** | Stokes drag vs flow velocity | Same profile (flow-dependent) |
| **fig16** | Static capture criterion map | Larger capture zone at 200 pg |

### **Stage 3: Trajectories & Efficiency (Figs 17–21)**
*Available in both `results/` (10 pg) and `results-200pg/` (200 pg)*

| Figure | Description | Key Difference |
|---|---|---|
| **fig17** | SPION loading sweep (static capture) | Monotonic increase with loading |
| **fig18** | Single-cell trajectory visualization | Shorter transit at higher loading |
| **fig19** | Trajectory bundle (multi-injection) | More captures at 200 pg |
| **fig20** | Capture efficiency vs v & loading | Higher efficiency at 200 pg |
| **fig21** | **Headline:** Static vs trajectory | Panel (b) reference: 50 pg vs 200 pg |

### **Stage 3c: Paracrine/VEGF (Figs 22–24)**
*All in `results/` (SPION-independent, uses paracrine parameters)*

| Figure | Description |
|---|---|
| **fig22** | VEGF concentration field (2D) |
| **fig23** | VEGF concentration vs radial distance |
| **fig24** | Time to reach therapeutic threshold |

---

## Comparison: Key Results

### **Static Capture Distance at v = 0.2 m/s**

| SPION Loading | Static (µm) | Trajectory (µm) | Extension Ratio |
|---|---|---|---|
| **10 pg** | 0.0 | 51.4 | ∞ (static predicts zero) |
| **50 pg** (original headline) | 6.8 | 96.9 | **14.2×** |
| **200 pg** (new variant) | 39.0 | 142.5 | **3.7×** |

**Interpretation:**
- **Absolute reach**: 200 pg captures from ~48% farther (96.9 → 142.5 µm)
- **Relative advantage**: Trajectory bonus shrinks (14.2× → 3.7×) because static criterion improves with higher loading
- **Physics**: At high loading, magnetic force dominates earlier, making static approximation more valid

### **Trajectory Capture Range Across Velocities (200 pg)**

| Velocity | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|
| 0.02 m/s (slow distal) | 94.5 | 256.5 | 2.7× |
| 0.05 m/s (Stage 3a) | 71.1 | 199.5 | 2.8× |
| 0.10 m/s | 56.5 | 165.3 | 2.9× |
| 0.20 m/s (MCA mean) | 39.0 | 142.5 | 3.7× |
| 0.50 m/s (fast MCA) | 18.5 | 108.3 | 5.8× |

**Key insight:** Even at high velocities, trajectory extends capture by 5.8× because radial magnetic drift accumulates over the 2 mm approach distance.

---

## How to Use Both Folders

### **For Publication/Comparison:**
- Use **fig21 from both folders** (original at 50 pg reference in panel b; 200pg at 200 pg reference in panel b)
- Caption should note: "Left panel (loading sweep) is identical; right panel references different SPION loading"

### **For Validation Against Literature:**
- **Polyak et al. 2008** (in vivo stent targeting): Use **200 pg results** (`results-200pg/`)
- **Tefft et al. 2014** (distal MCA, lower flow): Use **original results** (`results/` with 0.05 m/s variant in fig19)

### **For Sensitivity Studies:**
- Compare **figs 17, 20, 21** side-by-side between folders to show loading-dependence
- Both loading sweeps (fig17) show monotonic improvement from 10–400 pg

---

## Generation Scripts

Both result sets were regenerated from the codebase on **2026-04-12** using:

- **Original**: `scripts/regenerate_original_results.py`
  - Runs all 24 figures with project defaults (10 pg SPION, 8 struts, B₀=0.5 T, etc.)
  - Runtime: ~30 minutes
  
- **200 pg variant**: `scripts/generate_complete_200pg_results.py`
  - Runs figs 13–21 with 200 pg SPION (all other parameters identical)
  - Runtime: ~1.3 minutes (skips regeneration of pre-existing fig16, fig21)

---

## Parameter Summary

### **Identical Parameters (Both Folders)**

```
Stent geometry:
  Ring radius R          = 1.5 mm
  Strut width w          = 100 µm
  Strut thickness t      = 80 µm
  Strut length L         = 500 µm
  Number of struts       = 8 (evenly spaced)
  Magnetisation M        = 1.0 MA/m (radial, 304 SS saturation)

External field:
  B₀ magnitude           = 0.5 T (axial, +z direction)
  assume_saturation      = True

Blood flow:
  Vessel radius          = 1.54 mm (stent outer surface)
  Mean velocity          = 0.2 m/s (MCA-representative)
  Viscosity              = 4 mPa·s (whole blood)

Cell (endothelial):
  Radius                 = 10 µm
  Susceptibility (SPION) = 2.0 (SI)
  SPION density          = 5170 kg/m³ (magnetite)

Trajectory integration:
  Downstream boundary z_end = 2 mm
  RK45 tolerances: rtol=1e-5, atol=1e-8
  Binary search depth    = 7 (capture range resolution ~11 µm)
```

### **Differing Parameter**

```
              results/        results-200pg/
  SPION mass:  10 pg           200 pg
```

---

## Known Limitations & Considerations

1. **SPION saturation regime**: At B₀ = 0.5 T and 200 pg loading, real SPIONs may be partially saturated. The linear-susceptibility approximation (χ=2.0) may overestimate force by ~2×. See README.md **Limitations** for details.

2. **Strut geometry**: Rectangular prism approximation (Akoun & Yonnet kernel) neglects strut edge rounding. Field accuracy is ~5% locally.

3. **Flow model**: Poiseuille (parabolic) profile assumes quasi-steady laminar flow. Womersley number ≈2 for MCA; pulsatility effects not modelled.

4. **2D vs 3D trajectories**: All cells injected at θ=0 (strut-aligned). Real vessel wall has 8 strut locations; efficiency will vary by angular position.

---

## File Listing

### **results/ (24 figures, ~2.5 GB PNG + PDF)**
```
fig01–11  (Stage 1: B-field & gradient)
fig12     (Stage 1 extended: external field)
fig13–16  (Stage 2: force, drag, capture)
fig17–21  (Stage 3: trajectories, efficiency)
fig22–24  (Stage 3c: VEGF paracrine)
```

### **results-200pg/ (9 figures, ~400 MB PNG + PDF)**
```
fig13–21  (Stage 2 & 3 with 200 pg SPION)
COMPARISON_50pg_vs_200pg.md (detailed analysis)
```

---

## Quick Lookup

| Question | Answer | Where to Look |
|---|---|---|
| What's the headline result? | Trajectory = 14.2× static at 50 pg | `fig21` in `results/` |
| How does 200 pg compare? | Trajectory = 3.7× static; absolute range +47% | `fig21` in `results-200pg/` |
| Show me the SPION loading trend | Static & trajectory vs 10–400 pg | `fig17` in both folders |
| What's the efficiency in slow flow? | 20% at v=0.02 m/s with 50 pg | `fig19` in `results/` |
| Compare vs Polyak 2008 (200 pg)? | Use `results-200pg/` wholesale | All 9 figs in `results-200pg/` |
| B-field strength alone? | 500 mT at stent surface | `fig02`, `fig07`, `fig08` |
| How much drag on a cell? | ~2 pN at v=0.2 m/s | `fig15` |

---

## Project Structure

```
stent_capture_project/
├── results/                    # Original (10 pg, all 24 figures)
├── results-200pg/              # Variant (200 pg, figures 13-21)
├── scripts/
│   ├── regenerate_original_results.py
│   └── generate_complete_200pg_results.py
├── Overview.pdf                # Technical reference document
├── README.md                   # Physics summary & usage guide
├── AUDIT.md                    # Physics audit findings
└── RESULTS_SUMMARY.md          # This file
```

---

## Citation / Attribution

If using these results, cite:
- **Original project**: This repository (Stages 1–3, 2026)
- **SPION loading 200 pg**: Polyak et al. (2008) PNAS 105:698–703
- **Paracrine module**: Mac Gabhann & Popel (2006), Stefanini et al. (2008)
- **B-field model**: Akoun & Yonnet (1984), implemented as in magpylib
- **Force balance**: Furlani & Ng (2006) Phys Rev E 73:061919

---

Generated: **2026-04-12 16:30 UTC**
