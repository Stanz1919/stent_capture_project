# MSC Variant: Parameters, Literature, and Key Results

## Cell Parameters

| Parameter | Value | Source |
|---|---|---|
| Cell radius | **12.5 µm** (diameter 25 µm) | User estimate; conservative culture value |
| SPION loading (default) | **25 pg Fe/cell** | Centre of 10–30 pg standard range (see below) |
| SPION loading (sweep) | 5–100 pg | MSC-relevant range |
| Cell density | 1050 kg/m³ | User-specified; **not currently in physics model** |
| SPION susceptibility | 2.0 (SI) | Magnetite — unchanged from endothelial variant |
| SPION material density | 5170 kg/m³ | Magnetite Fe₃O₄ — unchanged |

---

## Literature: SPION Labelling of MSCs

### Key References

1. **Arbab et al. (2003)** *Radiology* 229(3):838–846  
   PubMed: 12819345  
   - Method: Ferumoxides + transfection agent (TA)  
   - **Loading: 30.1 ± 3.7 pg Fe/cell**  
   - Cell viability: 103–123% of control at 9 days post-labelling  
   - ~100% labelling efficiency  
   - *Upper end of practical range; transfection-assisted protocol*

2. **Farrell et al. (2008)** *(poly-L-lysine assisted)*  
   PubMed: 23577035 (related; see also Walczak 2008)  
   - Method: Iron oxide nanoparticles + poly-L-lysine  
   - **Loading: ~22 pg Fe/cell**  
   - Quantified by inductively-coupled mass spectrometry  
   - No adverse effects on MSC viability or differentiation  

3. **Saldanha et al. (2021)** *PubMed* 33747602  
   *Biodistribution of poly clustered SPION-labelled MSCs*  
   - Method: Resovist and CMF particles, passive incubation  
   - **Loading: 1.5–5.1 pg Fe/cell** (at 21-hour incubation, 10 µg USPIO/10⁵ cells)  
   - *Lower end — passive incubation without enhancement*

### Summary of Loading Range

| Protocol | Fe loading (pg/cell) | Notes |
|---|---|---|
| Passive incubation | 1.5–5 pg | Low uptake |
| Poly-L-lysine assisted | 10–22 pg | No cytotoxicity |
| Magnetofection (TA) | 20–30 pg | ~100% efficiency |
| **Model default (this work)** | **25 pg** | Centre of TA range |
| Upper range used in sweep | 100 pg | Theoretical maximum for sensitivity |

---

## Physics Implications vs Endothelial Cells

### At Equal SPION Loading (25 pg):

| Property | Endothelial (10 µm radius) | MSC (12.5 µm radius) | Ratio |
|---|---|---|---|
| Stokes drag | F_drag(EC) | 1.25 × F_drag(EC) | **+25%** |
| Magnetic force | F_mag | F_mag (same SPION mass) | 1× |
| Net: harder to capture | — | 25% more drag at equal loading | MSC disadvantaged |

### At Clinically-Comparable Conditions (25 pg MSC vs 10 pg EC):

| Property | Endothelial (10 pg) | MSC (25 pg) |
|---|---|---|
| Magnetic force ratio | 1× | ~2.5× |
| Drag ratio | 1× | 1.25× |
| Static capture at v=0.2 m/s | 0.0 µm | 0.0 µm |
| Trajectory capture at v=0.2 m/s | 51.4 µm | **65.0 µm** |

---

## Key Results (from fig21)

### Loading Sweep at v = 0.2 m/s (MSC, r = 12.5 µm)

| Loading (pg) | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|
| 5 | 0.0 | 32.5 | ∞ |
| 10 | 0.0 | 43.3 | ∞ |
| 15 | 0.0 | 54.2 | ∞ |
| **25** (default) | **0.0** | **65.0** | **∞** |
| 50 | 0.0 | 86.7 | ∞ |
| 75 | 12.7 | 97.5 | 7.7× |
| 100 | 18.5 | 97.5 | 5.3× |

**Key finding:** Static criterion predicts zero capture for MSCs below ~75 pg loading. Trajectory integration shows cells ARE captured (32–86 µm range), demonstrating the trajectory advantage is **critical** for MSCs.

### Velocity Sweep at 25 pg (MSC default)

| Velocity (m/s) | Static (µm) | Trajectory (µm) | Ratio |
|---|---|---|---|
| 0.020 (slow distal) | 39.0 | 130.0 | 3.3× |
| 0.050 | 18.5 | 97.5 | 5.3× |
| 0.100 | 0.0 | 86.7 | ∞ |
| **0.200** (MCA mean) | **0.0** | **65.0** | **∞** |
| 0.500 (fast MCA) | 0.0 | 43.3 | ∞ |

**Key finding:** At all MCA-relevant velocities (≥0.1 m/s), the static criterion fails entirely for 25 pg MSCs. Trajectory integration is the only valid predictor.

---

## Comparison: All Three Variants at v = 0.2 m/s

| Cell type | Loading | Radius | Static (µm) | Trajectory (µm) |
|---|---|---|---|---|
| Endothelial | 10 pg | 10 µm | 0.0 | 51.4 |
| Endothelial | 50 pg | 10 µm | 6.8 | 96.9 |
| Endothelial | 200 pg | 10 µm | 39.0 | 142.5 |
| **MSC** | **25 pg** | **12.5 µm** | **0.0** | **65.0** |
| **MSC** | **75 pg** | **12.5 µm** | **12.7** | **97.5** |

---

## Modelling Limitations for MSCs

1. **Cell density not modelled**: MSC density (1050 kg/m³) is slightly lower than blood (1060 kg/m³).
   Buoyancy force ≈ (1060 − 1050) × 9.81 × (4/3)π(12.5×10⁻⁶)³ ≈ **0.08 fN**.
   This is ~10⁵× smaller than Stokes drag at v=0.2 m/s (~2.6 pN) — safely negligible.

2. **Cell deformability**: MSCs deform under shear stress; rigid-sphere Stokes drag may
   overestimate drag by ~5–15% depending on deformability index.

3. **Non-spherical morphology**: MSCs in culture are often spindle-shaped (aspect ratio ~3).
   Effective drag depends on orientation; spherical model gives order-of-magnitude estimate.

4. **SPION distribution**: Model assumes uniform intracellular SPION distribution.
   Real MSCs may cluster SPIONs in endosomes, affecting local force.

---

## Generation Details

- Script: `scripts/generate_msc_results.py`
- Date: 2026-04-12
- Runtime: ~4 min (including fig21 trajectory sweeps at 7-step binary search resolution)
- Python: Anaconda 3.12 + scipy 1.13.1
