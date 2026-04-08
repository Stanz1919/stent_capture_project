# Changelog

All notable changes are documented here.

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
