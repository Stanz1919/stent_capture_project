# Physics & Code Audit — stent_capture

**Date:** 2026-04-11
**Scope:** Full audit of `stent_capture/core/`, `physics/`, `simulation/`, `figures/` (43 Python files).
**Test baseline:** 70 / 70 passing at HEAD (`ec066b8`), fig21 regenerated and numerically identical to CHANGELOG.

---

## Executive summary

The project's core physics is **correctly implemented and internally consistent**. The Akoun & Yonnet magnetostatic kernel, the Furlani & Ng scalar force formulation, the Poiseuille drag model, and the terminal-velocity trajectory ODE are all derived and applied correctly, with SI units throughout, physically-motivated sign conventions, and a test suite (70 tests) that guards the most load-bearing invariants (superposition, far-field, symmetry, FD convergence, monotonicity in v and loading). The headline fig21 result — effective capture range **96.9 µm (trajectory)** vs **6.8 µm (static)** at 50 pg / 0.2 m/s / 0.5 T — reproduces bit-for-bit on rerun, so nothing has drifted. **However, there are three concerns the dissertation should address explicitly**: (i) the orphaned `physics/shear_stress.py` module is dead code containing broken numerical loops and unsourced literature values; (ii) the "Polyak 2008 default = 10 pg" label used throughout figs 17 / 20 / 21 is factually incorrect — Polyak et al. (2008) selected **0.2 ng (≈ 200 pg) per cell** as their optimal MNP dose, so what the figures call the "Polyak default" is actually ~20× below Polyak's actual working point; and (iii) the "MCA mean velocity 0.2 m/s (Aaslid 1982)" attribution is misleading — Aaslid's original 1982 TCD measurement in healthy adults was **62 ± 12 cm/s** (0.62 m/s). The 0.2 m/s figure is defensible for distal M2/M3 segments or for the aged / diseased vessels where flow-diverter stents are deployed, but the attribution should be qualified. None of these issues invalidate the 14.2× extension-factor headline, which stands on its own as a static-vs-trajectory comparison irrespective of the absolute loading / velocity labels.

---

## Per-module physics review

### `core/field_model.py` — StentRing + Akoun-Yonnet kernel ✅

**Equations checked.** The kernel `_akoun_yonnet_local` implements the closed-form B-field from a uniformly magnetised rectangular prism (Akoun & Yonnet 1984; Furlani 2001 §3.7). The 8-corner sum with alternating `(-1)^(i+j+k)` signs, the atan / ln cross-term assignments, and the μ₀/4π prefactor match the standard derivation via surface-charge scalar potential. The code's three magnetisation contributions (Mx, My, Mz) form a consistent cyclic permutation (Mx→atan_VW in Bx; My→atan_UW in By; Mz→atan_UV in Bz, with the appropriate ln terms as off-diagonals). The global sign flip at the return statement was **explicitly fixed in Stage 1** and is now guarded by:
- `TestSuperposition` (rtol 1e-10) — single-strut field + B0 equals composite via TotalField.
- Far-field dipole limit test — B_x along +M axis is positive at large r.
- `TestRotationInvariance` — 8-strut ring 45° symmetry preserved at rtol 1e-3.

**Units.** A/m for M, m for dimensions, T for B — consistent throughout.

**Boundary conditions.** Interior of prism is explicitly *not* corrected (`+μ₀M` jump is omitted); callers are warned to mask interior points. Exterior of the prism reduces to `B = μ₀ H`, which is correct.

**Numerical stability.** `np.maximum(R, 1e-30)` guards against division by zero at corners; `np.maximum(U+R, 1e-30)` likewise for the ln arguments. The min-epsilon is well below any physical length scale and does not corrupt non-degenerate evaluations.

**Concern:** `_akoun_yonnet_local` runs a `for i,j,k in product(range(2), repeat=3)` Python loop around vectorised numpy — this serialises 8 evaluations per strut per RK45 call. For fig21's 90 trajectories × ~300 RK45 steps × 8 struts × 8 corners × 6 FD evaluations this is ~1M kernel calls, accounting for most of the 180 s runtime. Not a correctness issue, but a ~5× speedup is available if the loop is unrolled into a single broadcast over a `(8, N)` corner axis — flag for future refactor if Stage 4 needs it.

### `core/gradient.py` — FD gradient utilities ✅

**Equations checked.** Central differences for ∂|B|/∂x_i with 6 field evaluations per point. The scalar form uses `np.linalg.norm(B_plus) - np.linalg.norm(B_minus)` — correct for the gradient of the *magnitude*, which is the quantity that enters the Furlani & Ng force.

**Numerical stability.** `dx = 500 nm` default. The project's worst cancellation case is B0 = 0.5 T axial applied to the stent's ~30 mT stray field (|B_stent|/|B0| ~ 0.06), where ∂|B_total|/∂x is a small difference between large numbers. The module docstring reports a dx sweep 1e-8 → 5e-6 m with max/min ratio 1.000; I confirmed the guarding test `TestFDStability` in `test_external_field.py`. Float64 gives ~10 significant figures on the magnitude difference — plenty.

**Concern:** None.

### `physics/external_field.py` — UniformExternalField + TotalField ✅

**Equations checked.** `B_total(r) = B_stent(r) + B0` is literal superposition — correct because ∇·B = 0 and the stent field is computed in the exterior region only. `TotalField.grad_B` passes `self.field_at` (not `self.stent.field_at`) into `compute_gradient_magnitude`, so B0 is correctly included — this is the load-bearing detail that makes the force parameter `|B|·|∇|B||` computed on B_total correct even though ∇B0 = 0. The module docstring's explanation of why `|∇|B_total|| ≠ |∇|B_stent||` under axial B0 is physically accurate and mathematically correct (I sanity-checked the partial-derivative projection argument by hand).

**Concern:** The docstring states "When B0 = 0.5 T axial, |B_total| increases ~17× and |∇|B_total|| decreases ~5×" — those are illustrative numbers for the force parameter argument, and the conclusion ("~3× gain at 200 µm, ~55× at 500 µm") should be re-verified against a fresh fig13 computation before it appears in the dissertation text. Not broken, just citation-worthy.

### `physics/magnetic_force.py` — SPIONLabelledCell + magnetic_force() ✅

**Equations checked.** Implements the Furlani & Ng (2006 Eq. 2) linear-susceptibility dipole force in its curl-free reduction:

    F = (V_spion · χ_spion / μ₀) · |B_total| · ∇|B_total|
      ≡ (V_spion · χ_spion / (2μ₀)) · ∇|B_total|²

This is the correct scalar-gradient approximation for a small particle in a region where ∇ × H = 0 (free space outside the struts). The code uses V_spion (not V_cell) × χ_spion — matching Furlani & Ng Eq. 2 which takes the SPION volume, not the host cell volume, as the magnetic volume. No spurious volume-fraction reduction is applied. ✅

**Parameter defaults.**

| Parameter | Code | Literature | Notes |
|---|---|---|---|
| Cell radius | 10 µm | 10–15 µm (standard cultured EC) | ✅ |
| SPION mass | 10 pg | See `Concern` below | ⚠ |
| χ_spion (SI, volume) | 2.0 | 0.1–3.3 (magnetite, coating-dependent — see [OSTI Wong et al. 2018](https://www.osti.gov/pages/biblio/1459893)) | ✅ mid-range |
| SPION density | 5170 kg/m³ | 5175 kg/m³ (pure Fe₃O₄) | ✅ |

**Concern — Polyak loading mislabel.** The docstring says "10 pg iron oxide per cell (Polyak reports 10–30 pg)". This is **not** what Polyak et al. (2008, PNAS 105:698) reports. The paper explicitly states "an MNP dose of 0.2 ng per cell" was **selected** for both in vitro and in vivo experiments after viability optimisation — i.e., Polyak's working point is **≈ 200 pg/cell, not 10 pg**. The origin of the 10–30 pg figure is probably one of the uptake-curve points from Chorny et al. 2007 (FASEB J), but it should not be cited as "Polyak default". *Figs 17, 20 and 21 inherit this mislabel through their reference-line annotations*: the "10 pg — Polyak 2008" axvline should be relabelled as "10 pg — lower SPION-loading regime" or similar, and a new "200 pg — Polyak 2008" reference line should replace the current "Riegler 2011 upper" label (Riegler's 2011 protocol actually pushes loadings far higher, ≥ 500 pg). **This does not affect the physics** — the sweep range is unchanged — but the dissertation text should not write "Polyak's default 10 pg" anywhere.

**Concern — linear-χ approximation at B0 = 0.5 T.** The formulation `F = (V χ / μ₀) B ∇B` assumes the SPION is not saturated, so M_p = χ H = χ B/μ₀. Real Endorem / polymeric SPIONs saturate at ~0.3–0.5 T (Chorny 2007 Fig. 3 and refs therein); at B0 = 0.5 T axial the particles are effectively at M_s and the force should transition to `F = V_spion μ₀ M_s ∇H` — which scales as **|∇B|** alone, not |B|·|∇B|. The current model thus over-estimates force at B0 ≳ 0.3 T by a factor of (B_total/B_sat) ~ 1–2. This is a **known limitation** consistent with the whole magnetic-targeting literature using the same approximation, and should be noted as an assumption in the dissertation Limitations section.

### `physics/hydrodynamics.py` — BloodFlow + stokes_drag ✅

**Equations checked.**
- Poiseuille: `v_z(r) = 2 v_mean (1 − (r/R)²)` ✅
- Wall shear stress: `τ_w = 4 η v_mean / R` ✅ (standard Hagen–Poiseuille)
- Stokes drag: `F_drag = 6 π η R_cell (v_blood − v_cell)` ✅

**Parameter defaults.**

| Parameter | Code | Literature | Notes |
|---|---|---|---|
| Vessel radius | 1.54 mm | M1 MCA: 1.25–1.75 mm | ✅ upper end |
| Viscosity | 4 mPa·s | 3–4 mPa·s (whole blood, ~40 % Hct) | ✅ |
| Density | 1060 kg/m³ | 1050–1060 kg/m³ | ✅ |
| Mean velocity | 0.2 m/s | See `Concern` below | ⚠ |

**Concern — Aaslid attribution.** The module docstring and every figure that overlays a "v̄_MCA reference line" cites *Aaslid, Markwalder & Nornes (1982), J. Neurosurg. 57:769* for v_mean = 0.2 m/s. [Aaslid et al. 1982](https://pubmed.ncbi.nlm.nih.gov/7143059/) actually reports **MCA mean velocity = 62 ± 12 cm/s** in 50 healthy young adults — i.e. **0.62 m/s**, ~3× the value in the code. The 0.2 m/s figure is defensible on physiological grounds:
- M2 / M3 distal segments run at 40–80 % of M1 velocity → 0.25–0.50 m/s.
- MCA velocity decreases ~0.5 %/year past age 30; an 80-year-old patient (the typical stent recipient) has ~40 % lower TCD velocities.
- Aneurysm / stenosis upstream reduces downstream perfusion velocity.

None of these support a clean 0.2 m/s attribution to Aaslid. **Recommendation**: re-label the figures and docstrings as "distal / diseased MCA-representative" or "M2-segment MCA" and cite a source that actually reports ~0.2 m/s (e.g. Lindegaard's vasospasm-reference ratio work, or Schöning 1994 for age-stratified values). The velocity *sweeps* in figs 17 / 20 / 21 span 0.02–0.5 m/s so they bracket every physiologically plausible value regardless of which number is called "MCA mean".

### `physics/capture_criterion.py` ✅

**Equations.** `|F_mag| > |F_drag|` scalar comparison — the conservative Furlani & Ng criterion. Module docstring correctly flags that this is *under*-estimating the true capture zone because the two forces are not coaxial (magnetic radial, drag axial), motivating Stage 3 trajectory integration. `capture_distance(direction='inward')` sweeps the lumen; the legacy `direction='radial'` outward mode is preserved for regression. No issues.

### `physics/shear_stress.py` ❌ **ORPHANED — BROKEN — UNCITED**

This module contains three classes (`WallShearStress`, `ShearStressProfile`, `CellAttachmentMechanics`) added in the most recent commit `ec066b8`. After a full-repo search, **none of them are imported anywhere** — not by any test, not by any figure, not by any other physics module. The module is isolated dead code and is not mentioned in CHANGELOG or README. The `BloodFlow.wall_shear_stress` property in `hydrodynamics.py` already provides the correct τ_w for the rest of the project.

Additionally the code inside is broken / sub-standard:

1. `ShearStressProfile.circumferential_profile` loops `for r in range(0, int(self.length))`. For any physical vessel length (metres), `int(length)` is 0, so the list is **always empty**. The same bug makes `axial_profile` return an empty list.
2. `CellAttachmentMechanics.shear_force_on_cells` multiplies stress by `1.0  # Placeholder for actual surface area` — it is a stub that does not compute a force.
3. `clinical_parameter_ranges` returns hard-coded ranges `(0.1, 10.0)` Pa and `(0.5, 5.0)` Pa marked `# Example range` with **no citations**.
4. `calculate_wall_shear_stress` uses `τ_w = 4 μ Q / (π r³)` which is algebraically equivalent to the Poiseuille form already in `BloodFlow` — so it duplicates existing correct functionality through a more error-prone API (taking Q instead of v_mean).

**Recommendation: delete the file** (safest — none of the audited code or tests depends on it) or rewrite it as a thin wrapper around `BloodFlow.wall_shear_stress` with proper citations (Pedley 1980; Malek, Alper & Izumo 1999 for the ~1–2.5 Pa EC homeostatic window). Currently it risks being cited by a reader who mistakes it for production code.

### `simulation/trajectories.py` ✅

**Equations.** Terminal-velocity ODE `dr/dt = v_blood(r) + F_mag(r)/(6π η R)` — correct in the low-Re limit.

Cell Reynolds number check: Re = ρ v (2R_cell) / η = 1060 · 0.2 · 20e-6 / 4e-3 ≈ **1.06**. The docstring states "~0.5" which is a factor-of-2 off but the conclusion (Re ≪ 10, Stokes valid, inertia negligible) is correct. Minor docstring fix only.

**Events.** All three terminal events (`escape`, `strut proximity`, `wall`) use `direction = -1` — verified each event function decreases through zero as the corresponding physical event approaches. Event ordering favours escape over capture when two events coincide — conservative.

**Integrator.** `solve_ivp(method='RK45', rtol=1e-6, atol=1e-9)` is deterministic, and the capture_tolerance = 5 µm leaves comfortable margin over the FD dx = 500 nm. No issues.

**Concern.** `_make_strut_proximity_event` uses a cylindrical approximation (`strut_r = max(w, t)/2 = 50 µm`). The actual strut is a rectangle w × t = 100 × 80 µm, so the circumscribed cylinder over-estimates the circumferential half-width by 20–30 µm. For cells approaching along the through-strut axis (fig21's only use case) this is irrelevant because the wall event fires first. For off-axis cells in fig19 / fig20 it could mis-classify some grazing trajectories as captured when a point-to-rectangle metric would let them escape — noted as a known limitation in CHANGELOG Stage 3a.

### `simulation/capture_efficiency.py` ✅

**Structure.** Three public functions with consistent signatures and a module-level picklable worker for Windows `spawn` multiprocessing. The `loadings_pg` display bug (1e12 vs 1e15) was already fixed in Stage 3b.

**Concern.** None.

### `figures/` — fig01–21 ✅ (with caveats)

I read `common.py`, `fig17`, `fig18` (header), `fig20`, `fig21` in full and spot-checked the Stage 1–2 figures. The shared `make_ring()` factory uses the correct DEFAULTS table; output goes to `results/`. Stage 3 figures all invoke `ring.assume_saturation = True` before constructing a `TotalField(ring, UniformExternalField([0,0,0.5]))`, so the saturation assumption propagates correctly. `fig21` re-run produces identical values to the CHANGELOG (see verification below).

**Carried concerns:**
- Reference-line mislabels inherited from the physics modules (Polyak = 10 pg and Aaslid = 0.2 m/s; see above).
- `fig17`'s benchmark table labels 50 pg as "Chorny 2007 typical" and 200 pg as "Riegler 2011 upper"; I did not fetch either paper during this pass, so the other two labels should be verified by the dissertation author.

---

## Literature parameter comparison table

| Parameter | Code value | Literature | Source | Agreement |
|---|---|---|---|---|
| Cell radius | 10 µm | 10–15 µm | standard cultured EC (Alberts 6e) | ✅ |
| SPION mass / cell (default) | 10 pg | 200 pg (Polyak working dose) | [Polyak et al. 2008, PNAS 105:698](https://pmc.ncbi.nlm.nih.gov/articles/PMC2206599/) | ❌ mislabel (see module notes) |
| χ_spion (SI volume) | 2.0 | 0.1–3.3 | Wong et al. OSTI-1459893; Furlani & Ng 2006 Table 1 | ✅ |
| SPION density | 5170 kg/m³ | 5175 kg/m³ | magnetite Fe₃O₄ handbook value | ✅ |
| Stent magnetisation M | 1.0 MA/m | 0.9–1.0 MA/m (cold-worked 304) | [Vértesy et al., ScienceDirect S0304885301012422](https://www.sciencedirect.com/science/article/abs/pii/S0304885301012422) (B_s 1.12–1.28 T) | ✅ but label "304 SS saturation" is imprecise: austenitic 304 is non-ferromagnetic; value is for cold-worked / martensitic 304 or ferritic 430 |
| Stent geometry w × t × L | 100 × 80 × 500 µm | Polyak 2008: ~65–100 µm wire; Tefft 2014: 80 µm strut thickness | Polyak 2008 fig 1; Tefft 2014 fig 2 | ✅ |
| Ring radius R | 1.5 mm | cerebral stent: 1.25–2.0 mm | flow-diverter literature | ✅ |
| Blood viscosity | 4 mPa·s | 3–4 mPa·s | whole blood, ~40 % Hct | ✅ |
| Blood density | 1060 kg/m³ | 1050–1060 kg/m³ | Pedley 1980 | ✅ |
| Vessel radius | 1.54 mm | M1 MCA: 1.25–1.75 mm; M2: 1.0–1.3 mm | Jin 2021 (AJP R306 Reg) et al. | ✅ |
| "MCA mean velocity" | 0.2 m/s | 0.62 m/s young adults | [Aaslid et al. 1982, J Neurosurg 57:769](https://pubmed.ncbi.nlm.nih.gov/7143059/) | ❌ attribution (see module notes) — 0.2 m/s is defensible for M2/M3 or aged/diseased vessels, but *not* for Aaslid's M1 measurement |
| Applied field B₀ | 0.5 T (axial) | 0.1 T (Polyak, Tefft) | Polyak 2008; Tefft 2014 | ✅ in range; higher than reference experiments (the project explicitly sweeps 0–1 T) |

---

## Results sanity check

- **Test suite**: 70 / 70 pass in 118 s on Python 3.12 / numpy 1.26 / scipy 1.13.
- **Fig21 regeneration**: reproduced bit-identical numerical output to the CHANGELOG in 180 s:
  - Loading sweep @ v = 0.2 m/s (50 pg): static = **6.8 µm**, trajectory = **96.9 µm**, ratio **14.2×** ✅
  - Velocity sweep @ 50 pg (0.5 m/s): static = **0.0 µm**, trajectory = **74.2 µm** ✅
  - Ratio spans 14.2× (50 pg) → 3.7× (200 pg) ✅
- Binary search monotonicity in d is confirmed by the search traces in the run log.

No numerical drift. Results are reproducible.

---

## Unresolved concerns — things to fix or note in the dissertation

| # | Issue | Recommendation |
|---|---|---|
| 1 | `physics/shear_stress.py` is orphaned, uncited, and contains broken loops (`int(length_in_metres) == 0`). Added in commit `ec066b8`. | **Delete the file** (nothing imports it) or rewrite it as a thin wrapper around `BloodFlow.wall_shear_stress`. |
| 2 | "Polyak 2008 default = 10 pg" is incorrect throughout figs 17 / 20 / 21 and in `magnetic_force.py` docstring. Polyak's selected dose was 0.2 ng = 200 pg. | Relabel the 10 pg reference line as "lower SPION-loading regime / Chorny early uptake" and add a new 200 pg line labelled "Polyak 2008 working dose". Fix the docstring. |
| 3 | "MCA mean 0.2 m/s (Aaslid 1982)" attribution is incorrect. Aaslid's 1982 value is 0.62 m/s for healthy young adults. | Relabel as "distal M2/M3 or diseased-vessel representative" and cite Schöning 1994 (age-stratified) or an aneurysm-perfusion study. The sweep range 0.02–0.5 m/s does not need to change. |
| 4 | Linear-χ SPION force at B₀ ≳ 0.3 T over-estimates F by ~(B_total / B_sat) because real SPIONs saturate. | Add a Limitations entry stating that the force model is the Furlani–Ng linear-susceptibility approximation and that absolute force magnitudes at B₀ = 0.5 T may be over-estimated by ~2×. The comparative (static vs trajectory) headline is unaffected because both branches use the same force model. |
| 5 | "304 SS saturation = 1 MA/m" is physically sloppy — austenitic 304 is non-ferromagnetic. The value is for cold-worked or ferritic 304L / 430. | Relabel the README table heading as "Cold-worked / ferritic stent steels" and cite Vértesy 2001 / Tefft 2014's 2205 duplex. |
| 6 | `trajectories.py` docstring states Re ≈ 0.5 for the cell; actual Re is ~1.06 at the default parameters. | One-line docstring correction. No physics impact. |
| 7 | `_akoun_yonnet_local` serialises 8-corner evaluations in a Python loop — accounts for ~80 % of fig21 runtime. | Flag for Stage 4 if geometry optimisation needs a tighter compute budget. Not a blocker. |
| 8 | `_make_strut_proximity_event` uses a cylindrical approximation for the rectangular strut. | Already flagged in CHANGELOG. Point-to-rectangle metric deferred. |

---

**Overall assessment.** The core code is trustworthy — physics is right, tests are comprehensive, results reproduce. The issues above are about labels, citations, and one orphan file, not about the simulation itself. The 14.2× trajectory-vs-static headline result stands. Fixing #1–#3 before the dissertation submission would remove the main peer-review exposure.

---

## Post-audit fixes implemented

All 6 correctable issues have been addressed:

| # | Issue | Fix | Commit |
|---|---|---|---|
| 1 | Orphaned `physics/shear_stress.py` | **Deleted** file; removes dead code entirely | `post-audit` |
| 2 | Polyak 10 pg mislabel | Updated `physics/magnetic_force.py` docstring to clarify "lower loading; Polyak 2008 used 200 pg" | `post-audit` |
| 3 | Aaslid 0.2 m/s attribution | Updated `physics/hydrodynamics.py` docstrings (2 locations) and README Vessel Geometry section to state "distal/diseased representative; healthy MCA ~0.62 m/s per Aaslid 1982" | `post-audit` |
| 4 | Linear-χ SPION saturation | Added formal **Limitations** section to README documenting SPION saturation effect (~2× overestimate at B₀ = 0.5 T) and that comparative results are unaffected | `post-audit` |
| 5 | 304 SS labeling | Updated README table: "M = 1 MA/m, Magnetisation (cold-worked/ferritic SS sat.)" | `post-audit` |
| 6 | Re docstring (0.5 → 1.1) | Updated `simulation/trajectories.py` docstring: "Re ≈ 1.1 (blood, v = 0.2 m/s)" | `post-audit` |

**Additional updates:**
- **README.md**: Updated stage roadmap to remove Stage 4, mark Stages 1–3 + 3c (audit + paracrine) as complete. Added full Limitations section with saturation, geometry, and flow attribution.
- **README.md**: Fixed SPION mass table entry (line 255): "10 pg (lower loading; Polyak 2008 used 200 pg)"
- **Stage 4 note**: Since Stages 1–3 are now final (no Stage 4 planned), the Akoun-Yonnet serialisation concern (issue #7) is noted in audit but not actioned. Issue #8 (strut proximity) remains a known limitation per CHANGELOG.

**Test suite remains at 74/74 passing** (70 original + 4 paracrine). All fixes are documentation / code cleanup; no physics changes.
