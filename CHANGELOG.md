# Changelog

All notable changes are documented here.

---

## [0.2.1] — 2026-04-09  *(Stage 2.5 — SPION loading sweep)*

### Added

**Figure (`stent_capture/figures/`)**
- `fig17_spion_loading_sweep.py` — 2×1: (a) capture distance from stent inner
  surface vs SPION loading (1–300 pg, log-spaced), three velocities, B0 = 0.5 T
  axial; overlays for inter-strut half-distance (589 µm), lumen radius (1460 µm),
  literature benchmark loadings (10/50/200 pg), and experimental range band
  (30–100 pg). (b) force ratio |F_mag|/|F_drag| at 5 µm from stent inner surface
  vs loading, log–log, with ratio = 1 threshold line.
  Field and gradient are precomputed once along the radial line; loading/velocity
  loops only scale the precomputed force-parameter array.

**Tests (`stent_capture/tests/test_capture_criterion.py`)**
- `TestLoadingSweep::test_capture_distance_monotonic_in_loading` — capture
  distance (inward sweep) is non-decreasing over 10→30→100→300 pg at v=0.05 m/s.
- `TestLoadingSweep::test_capture_distance_positive_at_100pg` — capture distance
  > 50 µm at 100 pg, v=0.05 m/s, B0=0.5T (validates non-trivial capture exists
  in upper experimental loading regime; actual value ≈ 56.5 µm).

Total tests: **60** (all passing).

### Changed

**`stent_capture/physics/capture_criterion.py`**
- `capture_distance()` gains `direction='inward'` option: sweeps from the stent
  inner surface (R − t/2) toward the vessel centre, which is the physically
  correct lumen geometry for capture-distance analysis.  The legacy `'radial'`
  mode (default, sweeps outward from stent outer surface) is preserved for
  backward compatibility.

### Physics notes — Stage 2.5 results

Capture distance table (µm from stent inner surface, B0 = 0.5 T axial):

| Loading | v = 0.05 m/s | v = 0.20 m/s | v = 0.50 m/s |
|---------|-------------|-------------|-------------|
| 10 pg   | 0 µm        | 0 µm        | 0 µm        |
| 50 pg   | 38.5 µm     | 5.7 µm      | 0 µm        |
| 200 pg  | 71.4 µm     | 38.5 µm     | 20.3 µm     |

- Static capture is only predicted above ~19 pg at the slowest flow (0.05 m/s).
- A capture distance of 100 µm is not reached within the 1–300 pg sweep range
  at any of the three physiological velocities.
- The 8-strut stent geometry concentrates capture within ~70 µm of the stent
  inner surface even at 200 pg; inter-strut half-distance (589 µm) and lumen
  radius (1460 µm) are far beyond the static capture envelope.
- These constraints motivate Stage 3: trajectory-based capture can be
  substantially larger than the static criterion predicts, as cells accumulate
  radial displacement over their transit time through the stent.

---

## [0.2.0] — 2026-04-09  *(Stage 2)*

### Added

**Physics modules (`stent_capture/physics/`)**
- `magnetic_force.py`:
  - `SPIONLabelledCell` — model cell with `radius`, `spion_mass`, `spion_susceptibility`,
    `spion_density`; `volume` and `spion_volume` properties.  Default: 10 µm radius,
    10 pg iron oxide (magnetite), χ = 2.0.
  - `magnetic_force(cell, total_field, points)` — returns (N, 3) force vectors in Newtons
    using the scalar-gradient approximation F = (V_spion · χ / µ₀) · |B| · ∇|B|
    (Furlani & Ng 2006).
- `hydrodynamics.py`:
  - `BloodFlow` — Poiseuille flow in a cylindrical vessel; `velocity_at`, `shear_rate_at`,
    `wall_shear_stress`.  Default: R_vessel = 1.54 mm, v_mean = 0.2 m/s, η = 4 mPa·s.
  - `stokes_drag(cell, blood_flow, points)` — returns (N, 3) drag vectors via
    F_drag = 6π η R_cell (v_blood − v_cell).
- `capture_criterion.py`:
  - `capture_map(cell, total_field, blood_flow, points)` — returns dict with F_mag_vec,
    F_drag_vec, F_mag, F_drag, margin, captured (Furlani & Ng scalar criterion).
  - `capture_distance(cell, total_field, blood_flow)` — outermost radius at which
    |F_mag| ≥ |F_drag|.

**Core gradient (`stent_capture/core/gradient.py`)**
- `compute_gradient_vector` — new function returning (N, 3) vector gradient ∇|B|;
  `compute_gradient_magnitude` now delegates to it.

**Figures (`stent_capture/figures/`)**
- `fig14_force_vs_distance.py` — 2×2: (a/b) |F_mag| vs distance for axial/transverse
  B0; (c) polar plot of |F_mag| vs angle at 200 µm showing 8-strut periodicity;
  (d) |F_mag| at 200 µm vs SPION load 1–100 pg.
- `fig15_drag_vs_velocity.py` — 2×2: (a) Stokes drag vs distance from vessel wall;
  (b) Poiseuille velocity profile; (c) drag vs velocity with F_mag reference lines;
  (d) wall shear stress vs velocity with physiological reference bands.
- `fig16_capture_map.py` — 1×3: 2D cross-section force ratio |F_mag|/|F_drag|
  maps (log₁₀ scale) for v_mean = 0.05, 0.2, 0.5 m/s; B0 = 0.5 T axial.
  Diverging colormap (red = drag dominates, green = mag dominates, white = balance).
  Static force balance yields no capture anywhere in the lumen at physiological
  coronary flow for 10 pg SPION-loaded cells.

**Tests (`stent_capture/tests/`)**
- `test_magnetic_force.py` — 15 tests: zero gradient → zero force, force direction
  toward strut, χ and V_spion linearity, cell properties, order-of-magnitude check.
- `test_hydrodynamics.py` — 12 tests: Poiseuille profile (centreline, wall, mean),
  shear rate, Stokes formula and scaling.
- `test_capture_criterion.py` — 10 tests: dict structure, capture at surface, no
  capture at centreline, range decreases with velocity, range increases with B0;
  new test documenting that max |F_mag|/|F_drag| < 1.0 in lumen at v_mean=0.05 m/s
  (regression guard for the no-static-capture finding).

Total tests: **58** (all passing).

### Physics notes

- Force uses only SPION volume (not cell volume) × χ_spion: consistent with
  Furlani & Ng 2006 Eq. 2 and Polyak 2008 experimental regime.
- Default spion_mass_per_cell = 10e-15 kg (10 pg); force at 100 µm with B0=0.5T ≈ 340 pN.
- Capture criterion is scalar |F_mag| > |F_drag| (conservative; Stage 3 will
  use directional trajectory integration).
- Wall shear stress range 0.5–5 Pa for coronary conditions (v_mean 0.05–0.5 m/s).
- **Static force balance yields no capture at physiological flow**: at v_mean =
  0.05–0.5 m/s, Stokes drag (nN scale) exceeds magnetic force (pN–nN) throughout
  the lumen for 10 pg SPION-loaded cells. Maximum force ratio at v_mean = 0.05 m/s
  with B0 = 0.5 T is ≈ 0.52 (< 1.0 everywhere in the lumen). This motivates Stage 3
  trajectory integration, where cells accumulate radial drift over their transit
  time even when the instantaneous magnetic force is weaker than drag.

### Not included (future stages)

- Cell trajectory ODE integration (Stage 3)
- Capture efficiency maps (Stage 3)
- Geometry optimisation (Stage 4)

---

## [0.1.0] — 2026-04-08  *(Stage 1)*

### Added

**Package structure**
- `stent_capture/` Python package with `core/`, `physics/`, `figures/`, `tests/`
  sub-packages replacing the original monolithic `stent_analysis.py`.
- `legacy_stent_analysis.py` copied to project root for reference.

**Core field model (`stent_capture/core/`)**
- `field_model.py`: `StentRing` class — exact 3-D B-field using the
  Akoun & Yonnet (1984) / Furlani (2001) closed-form analytical expressions
  for uniformly magnetised rectangular prisms.
  - `field_at(points: (N,3))` primary interface returning `(N, 3)` B-vectors.
  - `B_field`, `B_magnitude`, `grad_B` convenience wrappers.
  - `assume_saturation: bool` flag to signal that M is the saturation value.
  - `mag_mode`: `'radial'`, `'circumferential'`, or `'axial'`.
- `gradient.py`: `compute_gradient_magnitude(field_func, points, dx)` —
  field-model-agnostic central finite-difference gradient; works with any
  callable returning `(N, 3)` B-vectors.

**External field physics (`stent_capture/physics/`)**
- `external_field.py`:
  - `UniformExternalField(B0_vector)` — spatially uniform B0 field with
    `field_at`, `magnitude`, and broadcast support.
  - `TotalField(stent_ring, external_field)` — composes stent + external
    field into a unified `B_total = B_stent + B0` object exposing the same
    API (`field_at`, `B_field`, `B_magnitude`, `grad_B`) for use in future
    force and trajectory stages.

**Figures (`stent_capture/figures/`)**
- `style.py`: dissertation-quality `rcParams` (serif font, 300 dpi, grid).
- `common.py`: `DEFAULTS`, `THRESHOLDS`, `make_ring()`, `save_fig()`.
- `fig01_single_strut.py` — single strut |B| and |∇B| vs distance; 3-D
  model at default L plus long-strut (L = 5 mm) for comparison.
- `fig02_ring_heatmaps.py` — ring cross-section |B| and |∇B| heatmaps.
- `fig03_gradient_vs_distance.py` — through-strut and between-struts profiles.
- `fig04_magnetisation_sweep.py` — M sweep and capture distance.
- `fig05_strut_dimensions.py` — thickness and width sweep.
- `fig06_n_struts.py` — number of struts and angular uniformity.
- `fig07_gradient_contours.py` — gradient heatmap with capture threshold contours.
- `fig08_force_parameter.py` — B·∇|B| force parameter.
- `fig09_axial_profile.py` — axial gradient profile for multiple L values.
- `fig10_convergence.py` — 3-D convergence vs strut length (replaces the
  former 2D-vs-3D comparison now that the 2-D model is removed).
- `fig11_rz_heatmap.py` — gradient in the r-z (axial cross-section) plane.
- `fig12_external_field_comparison.py` *(new)* — 2×2 panel:
  (a) B0 = 0 reference profile;
  (b) B0 = 0.5 T axial;
  (c) B0 = 0.5 T transverse (parallel to M);
  (d) capture distance vs B0 magnitude (0–1 T, axial direction).

**Tests (`stent_capture/tests/`)**
- `test_external_field.py`: 17 pytest tests covering:
  1. Regression — `TotalField(ring, None)` reproduces `StentRing.grad_B` at
     20 reference points to rtol = 1e-6.
  2. Superposition — single-strut |B_total| = |B_stent| + B0 when B0 ∥ M,
     verified to rtol = 1e-10.
  3. Rotation invariance — 45° symmetry of 8-strut ring conserved by
     `TotalField.grad_B` to rtol = 1e-3 (finite-difference precision limit).
  4. Far-field limit — |B_total| → |B0| and |∇B_total| → 0 at r = 50 mm.
  5. `UniformExternalField` — broadcast shape, magnitude, invalid input.

### Fixed

- **Akoun & Yonnet kernel sign** (`core/field_model.py`): the raw corner sum
  in `_akoun_yonnet_local` accumulated a global sign inversion relative to the
  `H = -∇φ_m` convention.  The negation was latent in the original
  `stent_analysis.py` but harmless there because all outputs used `|B|` and
  `|∇|B||`.  The sign is now corrected so that `B_stent + B0` is physically
  correct for external-field superposition.  All figure outputs are unchanged.

### Changed

- The 2-D surface-charge model (`B_rect_charge_2D`, `StentRing2D`) is removed.
  Figures 1–8 now use the 3-D Akoun & Yonnet model evaluated at z = 0 for the
  cross-section (midplane) slice.  Values differ slightly from the original
  2-D figures because the 3-D model accounts for the finite axial length
  (L = 500 µm) of the struts.
- `fig10` repurposed from "2-D vs 3-D comparison" to "3-D convergence check"
  showing the model approaches the infinite-length limit as L increases.

### Not included (future stages)

- Magnetic force on SPION-labelled cells (Stage 2)
- M(H) hysteresis curve beyond the `assume_saturation` flag (Stage 2)
- Cell trajectory integration and capture efficiency (Stage 3)
- Geometry optimisation (Stage 4)
