# Polyak (2008) Comparison Report — B₀ = 0.1 T

**Generated:** 2026-04-12
**Cell:** EC, radius = 10 µm (BAEC representative)
**B₀ (this run):** 0.1 T  (1,000 G — Polyak 2008 experimental)
**B₀ (code standard):** 0.5 T  (5× stronger)
**Binary search depth:** 10 iterations (~1.4 µm precision)

---

## Loading Sweep at v = 0.2 m/s

| Loading (pg) | Static 0.1T (µm) | Static 0.5T (µm) | Traj 0.1T (µm) | Traj 0.5T (µm) | Traj ratio 0.5T/0.1T |
|---|---|---|---|---|---|
| 10 | 0.0 | 0.0 | 51.5 | 51.4 | 1.00× |
| 50 | 6.8 | 6.8 | 84.0 | 96.9 | 1.15× |
| 100 | 24.4 | 24.4 | 101.6 | 119.7 | 1.18× |
| 200 | 39.0 | 39.0 | 119.2 | 142.5 | 1.20× |
| 300 | 47.7 | 47.7 | 131.4 | 153.9 | 1.17× |

**Key finding:** Code (0.5 T) overestimates trajectory capture distance relative to Polyak's experimental field (0.1 T).

---

## Velocity Sweep at 200 pg (Polyak dose)

| Velocity (m/s) | Static 0.1T (µm) | Static 0.5T (µm) | Traj 0.1T (µm) | Traj 0.5T (µm) | Traj ratio 0.5T/0.1T |
|---|---|---|---|---|---|
| 0.020 | 94.5 | 94.5 | 199.1 | 256.5 | 1.29× |
| 0.050 | 71.1 | 71.1 | 162.5 | 199.5 | 1.23× |
| 0.100 | 53.6 | 56.5 | 139.5 | 165.3 | 1.18× |
| 0.200 | 39.0 | 39.0 | 119.2 | 142.5 | 1.20× |
| 0.500 | 18.5 | 18.5 | 96.2 | 108.3 | 1.13× |

---

## Interpretation

### Why the results differ

The magnetic force is:
```
F = (V_spion × χ / μ₀) × |B_total| × |∇|B_total||
```

At distances far from the strut, |B_total| ≈ B₀. The 5× increase in B₀
(0.1 → 0.5 T) increases |B_total| by ~5×. The gradient |∇|B_total|| also
changes because B_total direction rotates. The net effect on capture distance
is visible in the tables above.

### Stent saturation note

Both runs use `assume_saturation = True` (M = 1.0×10⁶ A/m). At 0.1 T,
Polyak reports ~80–90% saturation for 304 SS — so the stent gradient field
is slightly overestimated (~10–20%) in the 0.1 T run. This means the 0.1 T
results here are a slight overestimate of Polyak's true experimental forces.

### SPION saturation note

χ = 2.0 is used at both B₀ values. Polyak reports ~80–90% SPION saturation
at 0.1 T. At 0.5 T (code standard), SPIONs are fully saturated, making the
linear χ approximation an overestimate — so the code's 0.5 T force is
doubly overestimated (higher B₀ AND fixed χ rather than saturated M_sat).

### Capture efficiency (Polyak: 20%)

Polyak measured 20% capture of 2.5×10⁶ cells at 0.1 T, 200 pg, 30 ml/min.
This code computes capture DISTANCE (µm), not efficiency (%). These cannot
be directly compared without the cell spatial distribution and stent area.

---

## Files Generated

- `figA_force_comparison.{png,pdf}` — Force vs distance at both field strengths
- `figB_loading_sweep_polyak_vs_code.{png,pdf}` — Loading sweep comparison
- `figC_velocity_sweep_polyak_vs_code.{png,pdf}` — Velocity sweep comparison