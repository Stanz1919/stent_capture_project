# MSC Version Summary
**Generated:** 2026-04-12  
**Branch:** MSCand200pg  
**Script:** [scripts/generate_msc_results.py](scripts/generate_msc_results.py)  
**Output:** [results-msc/](results-msc/)

---

## 1. What This Variant Does

The MSC variant re-runs Stage 2 & 3 figures (13–21) using **Mesenchymal Stem Cells (MSCs)** instead of endothelial cells. All stent geometry, B-field, and blood flow parameters are held identical to the original project defaults; only the cell parameters change.

The core scientific question being addressed: **can a magnetic stent capture MSCs under physiological flow conditions, and does trajectory integration change the answer compared to the static force-balance criterion?**

---

## 2. Cell Parameters

| Parameter | MSC (this work) | Endothelial (original) | Source |
|---|---|---|---|
| Cell radius | **12.5 µm** (25 µm diameter) | 10 µm | Conservative culture estimate |
| SPION loading (default) | **25 pg Fe/cell** | 10 pg (lower regime) | Centre of 10–30 pg TA-assisted range |
| Loading sweep range | 5–100 pg | 10–400 pg | MSC-relevant clinical range |
| Cell density | 1050 kg/m³ | — | **Documented only; not modelled in physics** |
| SPION susceptibility χ | 2.0 (SI) | 2.0 (SI) | Magnetite — unchanged |
| SPION density | 5170 kg/m³ | 5170 kg/m³ | Magnetite Fe₃O₄ — unchanged |

**Key literature basis for 25 pg default:**
- Arbab et al. (2003): 30.1 ± 3.7 pg with transfection agent
- Farrell et al. (2008): ~22 pg with poly-L-lysine
- Saldanha et al. (2021): 1.5–5 pg (passive incubation — lower bound)

---

## 3. Physics Differences vs Endothelial Variant

| Effect | Change | Direction |
|---|---|---|
| Stokes drag | +25% (radius 10 → 12.5 µm) | Harder to capture |
| Magnetic force (at 25 pg vs 10 pg EC) | +2.5× (loading ratio) | Easier to capture |
| Magnetic force (at 25 pg vs 200 pg EC) | −8× (loading ratio) | Much harder to capture |
| Net vs 10 pg EC | Force +2.5×, Drag +1.25× | **MSC advantaged** vs low-loading EC |
| Net vs 200 pg EC | Force −8×, Drag +1.25× | **MSC severely disadvantaged** vs high-loading EC |

The increased cell radius enters only through Stokes drag (F = 6πηRv). The SPION mass enters both the magnetic force (F ∝ V_SPION ∝ m_SPION) and the Stokes drag term (no dependence). This means doubling the SPION loading always helps capture; doubling the cell radius always hurts.

---

## 4. Figure-by-Figure Summary

### Fig 13 — Force Parameter |B|·|∇B| Map
**File:** [results-msc/fig13_force_parameter_msc.png](results-msc/fig13_force_parameter_msc.png)

This is a **purely field-dependent quantity** — it does not depend on cell type or SPION loading. The map is **identical to the original fig13** in `results/`. The figure is included for completeness and cross-reference but contains no MSC-specific information. The force parameter peaks at the stent strut surface (~88,000 T²/m near the strut tip) and drops exponentially into the lumen.

> **Flag:** This figure is redundant with `results/fig13_force_parameter.png`. It is not wrong, but it consumes space in the MSC results folder without adding new information. Any caption should make clear the map is cell-type-independent.

---

### Fig 14 — Magnetic Force vs Radial Distance
**File:** [results-msc/fig14_force_vs_distance_msc.png](results-msc/fig14_force_vs_distance_msc.png)

Compares the magnetic force on a 25 pg MSC (green) against a 10 pg endothelial cell (grey dashed). The MSC force is ~2.5× higher throughout, consistent with the 25/10 = 2.5 loading ratio (force scales linearly with SPION mass). Both curves peak at the stent inner surface (~1.46 mm) and fall steeply into the lumen on a log scale.

> **Flag:** The endothelial reference used here is the **10 pg lower-loading baseline**, not the 200 pg Polyak working dose. A reader comparing these figures may underestimate how much harder the MSC is to capture than a heavily labelled endothelial cell. The comparison is not wrong, but a note or additional reference line at 200 pg endothelial would be useful for context.

---

### Fig 15 — Stokes Drag vs Flow Velocity
**File:** [results-msc/fig15_drag_vs_velocity_msc.png](results-msc/fig15_drag_vs_velocity_msc.png)

Shows MSC drag (r = 12.5 µm) vs endothelial drag (r = 10 µm) across 0.01–1.0 m/s. The drag ratio is exactly 1.25× across all velocities (Stokes drag scales linearly with radius). The annotation in the figure correctly states "Drag ratio MSC/EPC = 1.25×". At v = 0.2 m/s and r = 1.5 mm (mid-lumen), MSC drag ≈ 3.3 pN.

No issues.

---

### Fig 16 — Static Capture Criterion Map (3-panel)
**File:** [results-msc/fig16_capture_map_msc.png](results-msc/fig16_capture_map_msc.png)

Three cross-sectional panels at v = 0.05, 0.20, 0.50 m/s. Colour shows log₁₀(F_mag/F_drag). All panels are predominantly red (drag dominates). The **maximum ratio across all panels is 1.04**, meaning the magnetic force barely exceeds drag only in a thin shell right at the strut tips, and only at the slowest velocity (0.05 m/s).

**Key result:** At 25 pg MSC and v ≥ 0.2 m/s, the static force-balance criterion predicts **zero capture everywhere in the lumen**. This motivates the trajectory analysis (figs 18–21).

> **Flag — minor cosmetic issue:** With max ratio = 1.04 at the strut surface (0.05 m/s panel), the green "capture possible" zone is essentially invisible at the resolution printed. A separate close-up inset near the strut tip would make this barely-positive region visible, otherwise the figure reads as "uniformly no capture" across all three panels, missing the nuance at slow flow.

---

### Fig 17 — SPION Loading Sweep (Static Criterion)
**File:** [results-msc/fig17_spion_loading_sweep_msc.png](results-msc/fig17_spion_loading_sweep_msc.png)

Log-log plot of static capture distance vs SPION loading (5–100 pg) at v = 0.2 m/s. Zero capture is explicitly annotated as "0 µm (no static capture)" for loadings 5–50 pg. The first non-zero values appear at 75 pg (12.7 µm) and 100 pg (18.5 µm). The 25 pg default is marked with a green dashed line — it sits firmly in the zero-capture zone.

**Key result:** The static criterion predicts zero capture for MSCs at clinical SPION loadings (≤ 50 pg) under MCA-representative flow, making trajectory integration not just useful but **essential** for any quantitative prediction.

No issues with this figure.

---

### Fig 18 — Single-Cell Trajectory
**File:** [results-msc/fig18_single_trajectory_msc.png](results-msc/fig18_single_trajectory_msc.png)

Trajectory of a single MSC injected at 50 µm from the stent inner surface (z = −2 mm upstream). Left panel: R-Z projection showing the cell drifting inward as it advects downstream. Right panel: radial position vs time — the cell reaches the stent surface (~1.46 mm) at t ≈ 28 ms. Status: **captured**.

The figure demonstrates that despite the static criterion predicting zero capture at this condition, the cell IS captured due to radial magnetic drift accumulated over the 2 mm upstream approach. This is the central mechanistic insight of the MSC variant.

No issues.

---

### Fig 19 — Trajectory Bundle
**File:** [results-msc/fig19_trajectory_bundle_msc.png](results-msc/fig19_trajectory_bundle_msc.png)

Bundle of 6 trajectories at fixed injection distances of 10, 20, 40, 60, 80, and 100 µm. The legend shows outcomes: 10, 20, 40, 60, and 80 µm all marked **OK (captured)**; 100 µm marked **X (escaped)**.

> **Flag — numerical inconsistency with fig21:** The trajectory bundle shows 80 µm as captured, but the binary search in fig21 reports a maximum capture range of **65 µm** at the same conditions (25 pg, v = 0.2 m/s). These results are contradictory. The stated binary search resolution is ~10.8 µm (7 iterations over r_inner × 0.95 ≈ 1.39 mm), giving an expected precision of ±11 µm. A result of 65 µm with 80 µm actually captured falls 15 µm outside this tolerance — slightly exceeding the expected precision band. The true capture boundary is likely in the **75–82 µm** range, and the binary search has underestimated it by ~10–17 µm. This is the most significant numerical inconsistency in the MSC results. **Recommendation:** Increase binary search depth from 7 to 10 iterations (resolution improves to ~1.4 µm) or cross-validate the fig21 binary search result against the direct trajectory test at 80 µm.

---

### Fig 20 — Capture Efficiency Trends
**File:** [results-msc/fig20_capture_efficiency_msc.png](results-msc/fig20_capture_efficiency_msc.png)

Two-panel log-log figure:
- **(a) Loading sweep** at v = 0.2 m/s: static capture distance vs 5–100 pg. Zero values are floored to 0.01 µm for log-axis plotting.
- **(b) Velocity sweep** at 25 pg: static capture distance vs 0.02–0.5 m/s. Shows a sharp falloff from ~200 µm at 0.02 m/s to near-zero above 0.1 m/s.

> **Flag — log-axis floor:** Panel (a) plots zero-capture cases at 0.01 µm (the `np.maximum(..., 0.01)` floor). Unlike fig17, these floored points are not annotated with "0 µm" labels, so a reader inspecting only fig20 might interpret the low-loading regime as having finite static capture (of order 0.01 µm) rather than genuinely zero capture. This is a cosmetic inconsistency with fig17's clearer treatment of zero values.

> **Note on panel (b):** The velocity sweep at 25 pg shows non-zero static capture only at 0.02 m/s (very slow distal flow). This is internally consistent with the MSC_NOTES table (18.5 µm static at v = 0.05 m/s, 39.0 µm at v = 0.02 m/s), but the figure's log-log scale compresses the physiologically relevant range (0.1–0.5 m/s) into a narrow band where all values appear near zero.

---

### Fig 21 — Static vs Trajectory Comparison (Headline Figure)
**File:** [results-msc/fig21_static_vs_trajectory_msc.png](results-msc/fig21_static_vs_trajectory_msc.png)

The headline comparison figure, two panels:

**Panel (a) — Loading sweep at v = 0.2 m/s:**

| SPION loading (pg) | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|
| 5 | 0.0 | 32.5 | ∞ |
| 10 | 0.0 | 43.3 | ∞ |
| 15 | 0.0 | 54.2 | ∞ |
| **25 (default)** | **0.0** | **65.0*** | **∞** |
| 50 | 0.0 | 86.7 | ∞ |
| 75 | 12.7 | 97.5 | 7.7× |
| 100 | 18.5 | 97.5 | 5.3× |

*See fig19 discrepancy note — true value may be ~75–82 µm.

**Panel (b) — Velocity sweep at 25 pg:**

| Velocity (m/s) | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|
| 0.020 | 39.0 | 130.0 | 3.3× |
| 0.050 | 18.5 | 97.5 | 5.3× |
| 0.100 | 0.0 | 86.7 | ∞ |
| **0.200 (MCA)** | **0.0** | **65.0** | **∞** |
| 0.500 | 0.0 | 43.3 | ∞ |

**Key finding:** At all MCA-relevant velocities (≥ 0.1 m/s) and at the clinical 25 pg loading, the static force-balance criterion predicts **zero capture**. Trajectory integration predicts 43–86 µm capture range over the same velocity range. The trajectory approach is not optional — it is the only method that gives a non-zero prediction.

> **Flag — suptitle "N/A":** The figure suptitle reads "At 25 pg, v = 0.2 m/s: trajectory extends static by **N/A**." This is because static = 0, making the ratio undefined (division by zero). The code outputs "N/A" as a fallback. For a dissertation figure, this should read "**∞** (static predicts zero capture)" to convey the physical meaning. The current "N/A" label looks like a formatting error rather than a meaningful physical statement.

---

## 5. Consolidated Results at Reference Conditions (v = 0.2 m/s, B₀ = 0.5 T)

| Cell type | Loading | Radius | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|---|---|
| Endothelial | 10 pg | 10 µm | 0.0 | 51.4 | ∞ |
| Endothelial | 50 pg | 10 µm | 6.8 | 96.9 | 14.2× |
| Endothelial | 200 pg | 10 µm | 39.0 | 142.5 | 3.7× |
| **MSC** | **25 pg** | **12.5 µm** | **0.0** | **65.0** | **∞** |
| **MSC** | **75 pg** | **12.5 µm** | **12.7** | **97.5** | **7.7×** |

**Interpretation:** A 25 pg MSC has a trajectory capture range (65 µm) intermediate between a 10 pg EC (51 µm) and a 50 pg EC (97 µm), despite its larger radius — because the 2.5× loading advantage outweighs the 1.25× drag disadvantage. The MSC is significantly harder to capture than a 200 pg EC (the Polyak 2008 experimental condition).

---

## 6. Flags and Inconsistencies — Summary

| # | Severity | Figure(s) | Issue | Recommendation |
|---|---|---|---|---|
| 1 | **Medium** | fig19 vs fig21 | Trajectory capture range discrepancy: fig19 shows 80 µm captured at 25 pg, v = 0.2 m/s; fig21 binary search reports 65 µm. Discrepancy (~15 µm) slightly exceeds the binary search resolution (~11 µm). | Increase binary search to 10 iterations (resolution ~1.4 µm) and regenerate fig21, or cross-validate by re-running a direct trajectory at 65, 70, 75, 80 µm. |
| 2 | **Low–Medium** | fig21 | Suptitle shows "N/A" for the trajectory-to-static ratio when static = 0. Looks like a formatting error rather than a physics result. | Change code to output "∞ (static = 0)" when `ref_static == 0`. |
| 3 | **Low** | fig13 | Force parameter map is cell-type-independent and identical to `results/fig13_force_parameter.png`. Inclusion in `results-msc/` is redundant. | Add a note in the figure caption or README: "Fig 13 is field-geometry-only and is reproduced here for completeness; it is identical to the original variant." |
| 4 | **Low** | fig20 | Panel (a) plots zero static capture values at 0.01 µm (log-axis floor) without annotation. Inconsistent with fig17, which explicitly labels them "0 µm (no static capture)". | Add annotations for floored-to-floor points, matching fig17's treatment. |
| 5 | **Low** | fig14 | Endothelial reference is 10 pg — the minimum loading regime — not the Polyak working dose (200 pg). | Add a second reference curve at 200 pg, or add a note clarifying the comparison is against the lower loading baseline. |
| 6 | **Low** | fig16 | The barely-positive capture zone at v = 0.05 m/s (ratio 1.04) is visually invisible in the 3-panel plot at the printed scale. | Add a fourth panel or inset at v = 0.05 m/s with tighter colorbar limits near ratio = 1 to make the thin capture zone visible. |
| 7 | **Info** | all | Cell density (1050 kg/m³) is documented in the script but not modelled in the physics. Buoyancy force ≈ 0.08 fN vs ~2.6 pN drag — safely negligible, as confirmed in MSC_NOTES. | No action needed; already documented in [MSC_NOTES.md](results-msc/MSC_NOTES.md). |
| 8 | **Inherited** | all | Aaslid 0.2 m/s attribution (partially fixed in docstrings post-audit, but some figure annotations may still read "MCA mean"). | Verify figure axis labels and reference-line labels read "distal/diseased MCA" not "healthy MCA (Aaslid 1982)". |
| 9 | **Inherited** | all | Linear-χ SPION force model overestimates F_mag by ~2× at B₀ = 0.5 T (SPION saturation). | Already documented in AUDIT.md and README Limitations. Confirm the MSC results section of the dissertation also states this limitation. |

---

## 7. Modelling Limitations (MSC-Specific)

These are additional limitations beyond those already documented in [AUDIT.md](AUDIT.md) for the original variant:

1. **Non-spherical morphology** — MSCs in culture are often spindle-shaped (aspect ratio ~3:1). Stokes drag for a prolate spheroid is orientation-dependent; the spherical approximation overestimates drag for cells aligned with flow by ~10–20%.

2. **Cell deformability** — MSCs are mechanically softer than endothelial cells. A deformable cell has an effective drag that is lower than the rigid-sphere Stokes value (Fung 1981), potentially by 5–15%.

3. **SPION clustering** — The model assumes SPION mass is distributed uniformly within the cell volume. Real MSCs accumulate SPIONs in endosomes, creating a localised high-permeability inclusion. The effective dipole approximation remains valid (SPION cluster << cell), but the polarisability tensor becomes anisotropic.

4. **Loading range** — The MSC sweep extends to 100 pg. The endothelial sweep goes to 400 pg. The cross-variant comparison (fig21 vs the equivalent 200 pg endothelial result) cannot be made directly from the sweep panels because the x-axis ranges differ. A combined plot including both cell types on the same axes would strengthen the comparison.

---

## 8. Files Generated

| File | Type | Content |
|---|---|---|
| [fig13_force_parameter_msc.png/pdf](results-msc/fig13_force_parameter_msc.png) | Stage 1 field | Force parameter |B|·|∇B| (cell-independent) |
| [fig14_force_vs_distance_msc.png/pdf](results-msc/fig14_force_vs_distance_msc.png) | Stage 2 force | F_mag vs radial distance; 25 pg MSC vs 10 pg EC |
| [fig15_drag_vs_velocity_msc.png/pdf](results-msc/fig15_drag_vs_velocity_msc.png) | Stage 2 drag | Stokes drag vs flow velocity; MSC vs EC |
| [fig16_capture_map_msc.png/pdf](results-msc/fig16_capture_map_msc.png) | Stage 2 map | Static capture criterion; 3-panel velocity sweep |
| [fig17_spion_loading_sweep_msc.png/pdf](results-msc/fig17_spion_loading_sweep_msc.png) | Stage 3 static | Static capture distance vs loading; 5–100 pg |
| [fig18_single_trajectory_msc.png/pdf](results-msc/fig18_single_trajectory_msc.png) | Stage 3 traj | Single MSC trajectory at 50 µm injection distance |
| [fig19_trajectory_bundle_msc.png/pdf](results-msc/fig19_trajectory_bundle_msc.png) | Stage 3 traj | Bundle of 6 trajectories at 10–100 µm injection |
| [fig20_capture_efficiency_msc.png/pdf](results-msc/fig20_capture_efficiency_msc.png) | Stage 3 static | Static capture efficiency vs loading and velocity |
| [fig21_static_vs_trajectory_msc.png/pdf](results-msc/fig21_static_vs_trajectory_msc.png) | Stage 3 headline | Static vs trajectory: loading sweep + velocity sweep |
| [MSC_NOTES.md](results-msc/MSC_NOTES.md) | Reference | Parameter sources, literature citations, physics notes |

**Not included in MSC results:** Stage 1 figures (figs 1–12, field-only, cell-independent) and Stage 3c paracrine figures (figs 22–24, VEGF model, cell-type-independent).

---

## 9. How to Reproduce

```bash
python scripts/generate_msc_results.py
```

Runtime: ~4 minutes (dominated by fig21 trajectory integration, 7×7 = 49 binary search trajectories).

Requirements: Anaconda 3.12, scipy 1.13.1, numpy, matplotlib (same as base project).
