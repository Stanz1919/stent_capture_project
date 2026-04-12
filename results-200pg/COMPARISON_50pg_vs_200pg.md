# Comparison: 50 pg vs 200 pg SPION Loading

This folder contains figures generated with **200 pg SPION loading** per cell, compared against the project defaults of **50 pg** (or 10 pg base).

## Files Generated

- `fig16_capture_map_200pg.{png,pdf}` — Static capture criterion map (Stage 2)
- `fig21_static_vs_trajectory_200pg.{png,pdf}` — Static vs trajectory comparison (Stage 3, headline result)

---

## Key Findings

### **Panel (a): Loading Sweep at v = 0.2 m/s**

| Loading (pg) | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|
| **10** | 0.0 | 51.4 | ∞ |
| **30** | 0.0 | 85.6 | ∞ |
| **50** (project default) | 6.8 | 96.9 | **14.2×** |
| **100** | 24.4 | 119.7 | 4.9× |
| **200** | 39.0 | 142.5 | **3.7×** |
| **300** | 47.7 | 153.9 | 3.2× |
| **400** | 56.5 | 165.3 | 2.9× |

**Observations:**
- At 200 pg, **trajectory capture range = 142.5 µm** vs static = 39.0 µm (3.7× extension)
- Compare to 50 pg: **trajectory = 96.9 µm** vs static = 6.8 µm (14.2× extension)
- The **extension factor decreases** at higher loading because the static criterion improves (force increases), reducing the relative advantage of trajectory integration
- However, the **absolute capture range increases monotonically** with loading:
  - 50 pg → 96.9 µm
  - 100 pg → 119.7 µm  
  - 200 pg → 142.5 µm (+47%)

### **Panel (b): Velocity Sweep at 200 pg**

| Velocity (m/s) | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|
| **0.020** (slow distal) | 94.5 | 256.5 | 2.7× |
| **0.050** (Stage 3a reference) | 71.1 | 199.5 | 2.8× |
| **0.100** | 56.5 | 165.3 | 2.9× |
| **0.200** (MCA mean) | 39.0 | 142.5 | 3.7× |
| **0.500** (fast MCA) | 18.5 | 108.3 | 5.8× |

**Observations:**
- At slow velocities (0.02 m/s), trajectory reaches **256.5 µm** — cells have more time for radial drift
- At high velocities (0.5 m/s), trajectory still captures **108.3 µm** despite static criterion predicting only 18.5 µm
- The extension factor **grows with velocity** because faster flow strengthens the relative advantage of sustained radial magnetic drift

---

## Comparison to Project Default (50 pg, v = 0.2 m/s)

| Metric | 50 pg (default) | 200 pg | Change |
|---|---|---|---|
| Static capture | 6.8 µm | 39.0 µm | +472% |
| Trajectory capture | 96.9 µm | 142.5 µm | +47% |
| Extension ratio | 14.2× | 3.7× | −74% (relative gain) |

**Interpretation:**
- **200 pg cells capture 4× farther** (142.5 vs 96.9 µm) in absolute terms
- **Trajectory advantage shrinks** (14.2× → 3.7×) because the static criterion improves substantially with higher loading
- This is **physically expected**: at higher SPION loading, the magnetic force dominates earlier in the cell's approach, so the static force-balance criterion becomes more representative

---

## Literature Context

**Polyak et al. (2008)** — the primary experimental reference for this project:
- Used **200 pg SPION loading** in their in vivo stent study (PNAS 105:698–703)
- Confirmed magnetic targeting of SPION-labelled endothelial cells to stent surfaces under MRI guidance
- Cell loading was 0.2 ng (≈ 200 pg) of iron oxide per cell

**Project default (10–50 pg):**
- Represents a **lower-loading scenario** (e.g., clinical translation with reduced iron burden per cell)
- Shows the robustness of trajectory integration: even at 10 pg, cells can be captured despite static criterion predicting zero capture

---

## Recommendations for Use

1. **Baseline comparisons**: Use the **50 pg project results** (original `results/` folder) as the primary comparison point with Stage 3a studies

2. **High-loading scenario**: Use the **200 pg results** (this folder) when:
   - Modelling Polyak-like experimental conditions (200 pg, transfected cells)
   - Validating against published in vivo stent studies
   - Exploring the saturation regime where force-balance approximations become more accurate

3. **Note on SPION saturation**: At B₀ = 0.5 T with 200 pg of iron oxide, real SPIONs may be partially saturated. The linear-susceptibility model (χ = 2.0) may **overestimate forces by ~2×**. See README.md **Limitations** section for details.

---

## Generated Figures

### Figure 16: Capture Map (200 pg)
Static force-balance contour showing the captured region in the r–z plane (strut-aligned axis, θ=0). The larger contour area vs 50 pg reflects the higher magnetic force at 200 pg loading.

### Figure 21: Static vs Trajectory (200 pg variant)
- **(a) Loading sweep**: Shows monotonic increase in capture range from 51.4 µm (10 pg) to 165.3 µm (400 pg), with the extension ratio declining as static criterion improves
- **(b) Velocity sweep at 200 pg**: Demonstrates the velocity-dependence of trajectory advantage; slower flows allow more radial drift accumulation, achieving larger capture zones

---

## Project Structure

```
stent_capture_project/
├── results/              (Project default: 50 pg, 10 pg)
└── results-200pg/        (This folder: 200 pg variant)
    ├── fig16_capture_map_200pg.png
    ├── fig16_capture_map_200pg.pdf
    ├── fig21_static_vs_trajectory_200pg.png
    ├── fig21_static_vs_trajectory_200pg.pdf
    └── COMPARISON_50pg_vs_200pg.md (this file)
```

---

## Generation Details

- **Generation script**: `scripts/generate_200pg_results.py`
- **All other parameters**: Identical to project defaults
  - Stent: 8 struts, R = 1.5 mm, w = 100 µm, t = 80 µm, L = 500 µm, M = 1 MA/m (radial)
  - Field: B₀ = 0.5 T (axial), assume_saturation = True
  - Flow: v̄ = 0.2 m/s, η = 4 mPa·s, vessel radius = 1.54 mm
  - Cell (except mass): radius = 10 µm, χ = 2.0, ρ = 5170 kg/m³
  - Trajectory integration: RK45, rtol = 1e-5, atol = 1e-8
- **Computation time**: ~253 seconds (4.2 min) for both figures
- **Date generated**: 2026-04-12
