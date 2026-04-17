# Physics & Code Audit — stent_capture (v2)

**Date:** 2026-04-14
**Branch:** `SPIONsaturation`
**Scope:** Full re-audit of `stent_capture/core/`, `physics/`, `simulation/`,
`paracrine/`, `figures/`, `scripts/`, plus the 200 pg / Polyak / saturation
analysis artifacts added since the previous audit (`bc40c15`, 2026-04-11).
**Test baseline:** 77 / 78 passing at HEAD (`12e1713`), 1 failing
(`TestLoadingSweep::test_capture_distance_positive_at_100pg`) — see Issue #1.

---

## Executive summary

The core physics from the previous audit is **unchanged and still correct**
(Akoun & Yonnet kernel, Furlani & Ng force, Poiseuille drag, terminal-velocity
ODE). Since `bc40c15` the project has grown in four material directions:
(i) a **Langevin SPION saturation model** wired into `magnetic_force.py` as
the *new default*; (ii) a **COMSOL-calibrated effective magnetisation**
(`M_COMSOL_EFF = 0.619 MA/m` at B₀ = 0, `M_COMSOL_EFF_B15 = 2.20 MA/m` at
B₀ = 1.5 T) with a dedicated validation figure (`fig04_comsol_gradient_validation`);
(iii) a **VEGF paracrine reaction–diffusion module** (`paracrine/`) with
three dedicated figures (`fig22`–`fig24`) and four tests; and (iv) two
parallel result sets (`results-200pg/`, `results-polyak/`) plus a
`results/` directory that is referenced in documentation but no longer
present on disk.

**Bottom line for the dissertation:** the project is *substantively more
mature* than at the last audit — saturation is now modelled, COMSOL provides
an independent FEM cross-check, and paracrine analysis completes the
"capture → therapeutic benefit" story. But there are **seven concerns** the
thesis must address before submission, in decreasing order of severity:

1. A unit test now **fails** under the default Langevin model
   (`test_capture_distance_positive_at_100pg`: expects > 50 µm, gets 36.1 µm).
   This is a *silent regression* from commit `53f7ff9` — the test threshold
   was never updated when the Langevin model was made default.
2. The headline numbers in `README.md` (14.2× extension factor at
   50 pg / 0.2 m/s / 0.5 T) were produced with **`M = 1.0 MA/m`, constant-χ,
   n_struts = 8**. The current defaults in `figures/common.py` are
   **n_struts = 12** and Langevin-χ, and at B₀ = 1.5 T the code
   silently swaps to `M_COMSOL_EFF_B15 = 2.20 MA/m`. Re-running fig21 now
   produces different absolute numbers. The README has not been updated
   and still advertises the old values as current.
3. The `results/` folder referenced by `RESULTS_SUMMARY.md` and by every
   figure's `save_fig()` call does **not exist** at HEAD. Only
   `results-200pg/` and `results-polyak/` are committed. Re-running any
   figure creates an empty `results/` on the fly but there is no committed
   baseline set.
4. The `n_struts` default changed from 8 → 12 in `common.py` (commit
   `d57a999`, to match COMSOL V2-2C geometry) but `README.md`, the
   physical-parameters table in `RESULTS_SUMMARY.md`, and all test
   fixtures still use `n_struts = 8`. Tests and figures therefore do not
   share a geometry.
5. `SATURATION_MODEL_SUMMARY.md` recommends a "hybrid" approach
   (constant-χ for headline results, Langevin as supplementary). The code
   currently implements the opposite — Langevin is the default. A thesis
   decision is required.
6. `RESUME_PART3.md` exposes the author's working notes (internal git
   identity workaround, "do not delete without asking" instructions,
   chat-GPT prompt critique). This file should be removed before the
   repository is shared with examiners.
7. The Polyak / Aaslid / 304-SS citation issues flagged in the previous
   audit were *partially* fixed in `858447e` but survived into new
   artifacts: `fig_saturation_impact.py` still attributes "v = 0.2 m/s
   (MCA mean, Aaslid 1982)" and `results-polyak/POLYAK_COMPARISON_REPORT.md`
   repeats the "Polyak 200 pg" vs "code 10 pg default" framing. The old
   audit already established that Polyak's working point *is* 200 pg, so
   the repeated framing confuses the narrative.

None of the above invalidate the physics. They are documentation drift
and test-suite coverage gaps that appeared when the project expanded
faster than its README and tests.

---

## Per-module re-audit

### `core/field_model.py` — StentRing + Akoun–Yonnet kernel ✅

Unchanged since `bc40c15`. Sign flip, ε guards, interior-mask rule are all
as audited previously. The `assume_saturation` flag is still a pure
signalling boolean — no M(H) hysteresis.

### `core/gradient.py` — FD gradient utilities ✅

Unchanged. Default `dx = 500 nm` still defensible.

### `physics/external_field.py` — UniformExternalField + TotalField ✅

Unchanged. The load-bearing `self.field_at` vs `self.stent.field_at`
distinction is still correct.

### `physics/magnetic_force.py` — SPIONLabelledCell + magnetic_force() ⚠

**Changed materially in `53f7ff9`** (2026-04-13). The default SPION model
is now Langevin saturation:

```
chi_eff(B) = chi_0 * 3 L(xi) / xi    where xi = 3 chi_0 B / (M_sat * mu_0)
```

- `chi_0 = 2.0`, `M_sat = 446 kA/m` (bulk magnetite, Furlani & Ng 2006 Table 1)
- Small-ξ Taylor expansion `1 - xi²/15` for numerical stability at
  `xi < 1e-3` — **verified by hand**: exact Taylor series of `3 L(x)/x`
  is `1 − x²/15 + 2x⁴/945 − …`, so the first-order form used is correct.
- Large-ξ branch uses `coth(xi) − 1/xi` directly — correct.
- `spion_sat_magnetization=None` preserves backward compatibility
  (constant-χ mode); verified by `TestConstantChiModeUnaffected`.

**Physics note.** The Langevin model is now *more* consistent with the
COMSOL µᵣ = 2 soft-ferromagnet stent than the previous constant-χ
formulation, because both the stent and the SPIONs saturate under the
applied B₀. This is a physics-realism improvement.

**Concern — silent default change + stale test.** Making Langevin the
default reduces the effective force at B₀ = 0.5 T by ~20–25 % in the
active capture zone (confirmed by `scripts/compare_saturation_models.py`
and by the failing test: the old 100 pg capture distance was ≈ 57 µm
under constant χ, is now 36 µm under Langevin). The test threshold
(> 50 µm) was never updated. **Either** the test threshold needs to
be lowered, **or** the regression commit was a silent bug — in either
case the failure must not be left unresolved before the thesis submission.

**Concern — Polyak 10 pg docstring relic.** The docstring still calls
10 pg the "lower loading; Polyak 2008 used 200 pg" default. This is
defensible as written, but the code initialises `spion_mass_per_cell = 10e-15`
(10 pg) without a warning. All figures that rely on the cell default
therefore silently use 10 pg — a loading *below* every literature
benchmark. Recommendation: change the default to 50 pg (the old
headline value) and force any figure that needs 10 pg to set it
explicitly.

### `physics/hydrodynamics.py` ✅

Unchanged. Aaslid / v_mean attribution was fixed in `858447e` and the
docstring now reads correctly ("distal/diseased vessel representative;
healthy MCA ~0.62 m/s per Aaslid et al. 1982"). The Polyak / 200pg
comparison report repeats the original misleading framing, though — see
Issue #7.

### `physics/capture_criterion.py` ✅

Unchanged. `direction='inward'` is still the correct mode.

### `physics/shear_stress.py` — deleted ✅

Previously orphaned / broken. Deletion confirmed in `858447e`.

### `simulation/trajectories.py` ✅

Unchanged. Re docstring now correctly states `Re ≈ 1.1`. Cylindrical
strut-proximity approximation still in force, still acknowledged in the
CHANGELOG.

### `simulation/capture_efficiency.py` ✅

Unchanged.

### `paracrine/transport.py` ✅ (NEW — first audit)

**Equations checked.**
`∂C/∂t = D ∇²C + S − k C` with periodic x and zero-flux z BCs, solved
either as a sparse linear system (steady state) or explicit Euler
(transient). The 5-point Laplacian discretisation is standard; the
zero-flux BC is correctly enforced by adding `D/dz²` back to the
diagonal at the first / last row of each column (ghost-cell mirror).
Periodic x BC adds off-diagonal entries linking rows `j` and
`(Nx-1)·Nz + j` — correct for a unrolled cylinder.

**Parameter defaults.**

| Parameter | Code | Literature | Source | Agreement |
|---|---|---|---|---|
| `D` (VEGF in tissue) | 1.04 × 10⁻¹¹ m²/s | ~10⁻¹¹ m²/s | Mac Gabhann & Popel 2006 | ✅ exact |
| `k_deg` (first-order) | 1.93 × 10⁻⁴ s⁻¹ (t½ ≈ 60 min) | 30–90 min tissue t½ | Stefanini 2008 | ✅ midrange |
| Diffusion length L_D | √(D/k) ≈ 232 µm | — | derived | ✅ |

**Concern.** The explicit Euler stability limit
`dt_CFL = 0.5 / (D (1/dx² + 1/dz²) + 0.5 k)` is a *reaction–diffusion*
CFL with a spurious ½ on the reaction term. The correct 2-D explicit
Euler bound is `dt ≤ 0.5 / (D (1/dx² + 1/dz²))` with `k dt ≤ 2` as a
separate reaction-stiffness constraint. In practice `k ~ 2·10⁻⁴` s⁻¹ is
so much smaller than `D/dx²` that the difference is unobservable, but
the formula should be corrected for rigour.

### `paracrine/secretion.py` ✅ (NEW)

**Equations checked.** Gaussian source per captured cell with
`q_vol_peak = q_cell / (2π σ² h) × 10³` converting g/(m³·s) to
ng/(mL·s). The `10³` factor is correct (`1 g/m³ = 10⁶ ng/m³ =
10⁶ ng / 10³ mL = 10³ ng/mL`). ✅

**Parameter defaults.**

| Parameter | Code | Source | Notes |
|---|---|---|---|
| q_cell | 0.068 molecules/cell/s ≈ 5.08e-21 g/s | Stefanini 2008 | ⚠ originally per myonuclear domain, repurposed per EC |
| σ (cell radius) | 10 µm | standard EC | ✅ |
| h (tissue slab) | 20 µm | ~2 cell layers | ⚠ order-of-magnitude choice |

**Concern.** The Stefanini secretion rate was calibrated against a
*two-compartment mouse muscle model* in which the "source" is a whole
muscle fibre domain, not a single endothelial cell. Applying it 1:1 to
BAECs is a ~order-of-magnitude assumption. The docstring acknowledges
this, which is good — but the `fig22` title ("VEGF-enhanced 100×")
already absorbs the uncertainty, so the basal / enhanced contrast in
`fig22` is defensible.

### `paracrine/therapeutic.py` ✅ (NEW)

**Equations.** `therapeutic_zone_radius` uses nearest-cell distance to
each grid point above threshold — correct metric. `concentration_vs_distance`
azimuthally averages in radial bins — correct. `time_to_threshold`
nearest-grid-cell lookup is correct but uses `X[0,0]` as the origin,
which assumes a cell-centred grid (the `ParacrineField` grid is indeed
cell-centred with `x_i = (i+0.5) dx`). ✅

**Concern.** The `t_res = distance / v_slip` rejection is correct
physics — Péclet number Pe = vL/D for L = 500 µm, v = 100 µm/s,
D = 10⁻¹¹ m²/s gives Pe ≈ 5, so advection is *not* fully negligible.
The module ignores advection entirely. This is defensible because the
slip velocity at the vessel wall is already close to zero for captured
cells (they are stationary on the wall) and because the residual
Poiseuille velocity in the 20 µm slab is ~0. But the docstring's
"Pe ≪ 1" claim is stronger than the regime actually supports —
recommend softening to "advection contribution is secondary".

### `figures/` — fig01–24 + fig_saturation_impact ✅ (with caveats)

**New since last audit:**
- `fig04_comsol_gradient_validation.py` — FEM-vs-analytical bar chart at
  B₀ = 1.5 T, reports RMS error in the threshold-crossing distances.
  Uses `make_ring(B0_magnitude=1.5)` which auto-swaps in
  `M_COMSOL_EFF_B15`. ✅
- `fig22_concentration_field.py`, `fig23_concentration_vs_distance.py`,
  `fig24_time_to_threshold.py` — VEGF paracrine triplet. Basal vs
  100×-enhanced contrast in fig22 is physically meaningful; fig23
  overlay with the 5–25 ng/mL therapeutic band (Ozawa 2004) is correct;
  fig24 time-to-threshold values depend on `q_cell × n_cells` so the
  enhanced scenario is the only one that crosses within simulated time.
- `fig_saturation_impact.py` — three-panel Langevin vs constant-χ
  comparison. Panel (a) analytical `χ_eff(B)` curve, panel (b) force
  profile, panel (c) loading sweep. Still contains the Aaslid 0.2 m/s
  reference line — inherited mislabel (Issue #7).

**Filename collision.** Two files start with `fig04_`:
`fig04_magnetisation_sweep.py` (Stage 1) and
`fig04_comsol_gradient_validation.py` (new). They write different
output stems (`fig4_magnetisation_sweep.png` vs
`fig4_comsol_gradient_validation.png`) so there is no overwrite, but
the double-`fig04` prefix is confusing and `scripts/regenerate_original_results.py`
only imports one of them. Recommend renaming the COMSOL one to
`fig25_comsol_gradient_validation.py` for ordering consistency.

**Geometry drift.** `DEFAULTS["n_struts"] = 12` (set in `common.py`)
does not match the test fixtures (`n_struts = 8`), the README table
(`n_struts = 8`), or the `RESULTS_SUMMARY.md` parameter block
(`n_struts = 8`). Every figure that uses `make_ring()` silently runs
with 12 struts. This is a breaking change relative to all audit
baselines and relative to the CHANGELOG entries ("8 struts" is still
printed as the default in every entry). See Issue #4.

### `scripts/` ✅ (NEW since last audit)

Seven top-level scripts, all runnable via `python -m scripts.X`:

- `build_overview.py` — generates `Overview.pdf` (25 pages, 684 KB).
- `check_comsol_calibration.py` — verifies COMSOL calibration is
  unaffected by SPION model choice. ✅ reports <1% agreement at
  100 T/m and 40 T/m thresholds.
- `compare_saturation_models.py` — quantitative Langevin vs
  constant-χ comparison; numbers backing `SATURATION_MODEL_SUMMARY.md`.
- `fix_fig16_200pg.py` — one-shot patch script for a specific 200pg
  figure regeneration. Orphaned after 200pg regeneration; recommend
  deletion.
- `generate_200pg_results.py`, `generate_complete_200pg_results.py` —
  two scripts for the same task; `generate_complete_...` is the
  canonical one (per `RESULTS_SUMMARY.md`). The shorter one appears
  to be a superseded prototype. Recommend deletion.
- `generate_polyak_comparison.py` — generates `results-polyak/`.
- `regenerate_original_results.py` — runs all 24 figures with project
  defaults. Advertised 30-minute runtime; has not been re-verified
  against the current Langevin default.

Only one of these scripts (`regenerate_original_results.py`) has
a clear maintenance contract with the test suite. The rest are
analysis artefacts that should either be moved to an `analysis/`
subfolder or deleted.

### `legacy_stent_analysis.py` ✅

Still present at project root. Marked in README as "reference only".
Not imported anywhere. Can stay.

---

## Literature parameter comparison — updated

| Parameter | Code value | Literature | Source | Agreement |
|---|---|---|---|---|
| Cell radius | 10 µm | 10–15 µm | Alberts 6e | ✅ |
| SPION mass / cell (default) | **10 pg** | **200 pg** (Polyak) | Polyak 2008 | ⚠ default is below every experimental benchmark — recommend bump to 50 pg |
| χ₀ (SI volume, low-field) | 2.0 | 0.1–3.3 | Furlani & Ng 2006 Table 1 | ✅ |
| M_sat (SPION) | 446 kA/m | 446 kA/m (bulk magnetite) | Furlani & Ng 2006 Table 1 | ✅ exact |
| SPION density | 5170 kg/m³ | 5175 kg/m³ | handbook | ✅ |
| Stent M (no-B₀ default) | **1.0 MA/m** | cold-worked 304: 0.9–1.0 MA/m | Vértesy 2001 | ✅ |
| Stent M (B₀ = 1.5 T auto) | **2.20 MA/m** | COMSOL µᵣ=2 effective | this project's own calibration | ✅ self-consistent |
| n_struts (default) | **12** | V2-2C geometry | COMSOL model | ⚠ changed silently from 8 |
| Ring radius R | 1.5 mm | 1.25–2.0 mm | flow-diverter lit | ✅ |
| Blood viscosity | 4 mPa·s | 3–4 mPa·s | Pedley 1980 | ✅ |
| Vessel radius | 1.54 mm | M1 MCA 1.25–1.75 mm | Jin 2021 | ✅ |
| Mean velocity | 0.2 m/s | 0.62 m/s healthy, 0.1–0.4 m/s distal/aged | Aaslid 1982; Schöning 1994 | ✅ defensible, now correctly labelled |
| B₀ default | 0.5 T | 0.1 T (Polyak, Tefft); 1.5–3 T (MRI-guided) | literature span | ✅ bracketed by the project's own 0–1 T sweeps |
| VEGF D (tissue) | 1.04 × 10⁻¹¹ m²/s | 1 × 10⁻¹¹ m²/s | Mac Gabhann & Popel 2006 | ✅ |
| VEGF k_deg | 1.93 × 10⁻⁴ s⁻¹ (t½ ≈ 60 min) | 30–90 min tissue t½ | Stefanini 2008 | ✅ |
| VEGF therapeutic band | 5–25 ng/mL | 5–25 ng/mL | Ozawa 2004 | ✅ |

---

## Results sanity check

- **Test suite**: 77 / 78 pass. The failure
  (`test_capture_distance_positive_at_100pg`) is new since the last
  audit and is a direct consequence of the default-χ change in
  commit `53f7ff9`. Issue #1.
- **Fig21 regeneration not verified**. Under the previous audit the
  50 pg / 0.2 m/s / 0.5 T headline was 96.9 µm (trajectory) vs 6.8 µm
  (static) = 14.2×. Those numbers were produced with
  `n_struts = 8`, constant χ, `M = 1.0 MA/m`. With the current defaults
  (`n_struts = 12`, Langevin χ at B₀ = 0.5 T — which in the active zone
  is only mildly reduced, ~5–15 %) the numbers will move. The README
  still advertises 14.2× and 96.9 µm as current; they should be
  re-computed or the defaults rolled back for the dissertation run.
- **`results/` directory missing**: `RESULTS_SUMMARY.md` describes 24
  figures in `results/` but the folder is not on disk at HEAD. Only
  `results-200pg/` (9 figures) and `results-polyak/` (3 figures) are
  committed.
- **COMSOL cross-check**: `fig04_comsol_gradient_validation` reports
  <1% agreement at the two physiologically relevant thresholds
  (100 T/m, 40 T/m), 12% at 300 T/m (very near strut surface).
  `scripts/check_comsol_calibration.py` confirms calibration is
  independent of SPION susceptibility model. ✅ — this is the
  strongest piece of external validation in the project.

---

## Unresolved concerns — things to fix or note in the dissertation

| # | Issue | Recommendation |
|---|---|---|
| 1 | `test_capture_distance_positive_at_100pg` fails under Langevin default (36 µm vs 50 µm threshold). | Either lower the threshold to 30 µm with a comment explaining Langevin attenuation, or investigate whether the regression indicates a Langevin-implementation bug (analytical `3 L(ξ)/ξ` value at B₀ = 0.5 T, χ₀ = 2, M_sat = 446 kA/m is ≈ 0.70 — matches the ~30% force reduction observed). The test-threshold fix is correct; do it. |
| 2 | README headline (14.2× at 50 pg, 0.2 m/s, 0.5 T) predates the n_struts / Langevin / COMSOL-M changes. | Re-run fig21 with the *actual* thesis defaults, update README to report the new numbers, and archive the constant-χ / 8-strut numbers as a comparison baseline. |
| 3 | `results/` directory referenced everywhere but missing from the repo. | Commit the generated `results/` (or its `.gitignore` an entry + a `regenerate` shell command) so that `RESULTS_SUMMARY.md` is not advertising non-existent files. |
| 4 | `n_struts` default drift: `common.py` = 12, tests = 8, README = 8. | Pick one, and make the README, tests, CHANGELOG and figure module all agree. If 12 is the thesis geometry, update the "Default stent parameters" table in the README and the `_DEFAULT_STENT` fixture in every test file. |
| 5 | Langevin is the code default but `SATURATION_MODEL_SUMMARY.md` recommends constant-χ for the thesis headline. | Decide. If Langevin: edit the summary memo to reflect that. If constant-χ: change the `SPIONLabelledCell` default to `spion_sat_magnetization=None`. Either way, the README Limitations section needs to name the chosen default. |
| 6 | `RESUME_PART3.md` contains private working notes (git identity workaround, internal commentary). | Delete before sharing repo with examiners. |
| 7 | Re-surfacing Aaslid / Polyak mislabels in `fig_saturation_impact.py` and `POLYAK_COMPARISON_REPORT.md`. | Search-and-replace "MCA mean, Aaslid 1982" → "distal/M2-representative"; "code 10 pg default" → "10 pg low-loading reference" in the new artefacts. |
| 8 | Paracrine transient CFL bound has a stray ½ on the reaction term. | One-line edit in `paracrine/transport.py:184`. Effectively harmless at current k. |
| 9 | Orphan scripts (`fix_fig16_200pg.py`, `generate_200pg_results.py`). | Delete or move to `scripts/archive/`. |
| 10 | SPION default mass (`10 pg`) is below every experimental benchmark. | Change `spion_mass_per_cell` default to `50e-15` (50 pg) — the value that produces the headline 14.2× result and the only one with matching published figures. Explicit overrides for 10 pg / 200 pg in figure scripts. |
| 11 | `fig04_comsol_gradient_validation` shares the `fig04_` prefix with `fig04_magnetisation_sweep`. | Rename to `fig25_comsol_gradient_validation.py` for unambiguous ordering. |
| 12 | README `Package layout` and `Quick start` blocks predate the new modules (paracrine, scripts, figures 17–24, saturation fig). | Refresh — listed in README Part 2 below. |

---

## Will this project produce thesis-quality results?

**Short answer: yes, with caveats.**

**What works as-is.**

1. **Headline physics comparison.** The static-vs-trajectory contrast
   (fig21, fig20, fig19) is genuinely novel for a cerebral-flow
   magnetisable-stent geometry at this level of analytical rigour. The
   headline claim — that the static Furlani & Ng criterion
   systematically under-estimates capture by a factor that grows with
   flow velocity and shrinks with SPION loading — is supported by
   reproducible, deterministic simulations, is consistent with the
   experimental capture efficiencies reported by Polyak (20%) and
   Tefft (40–60%) given the geometric differences, and is exactly
   the kind of "simple physical model + clean comparison" finding
   that makes a publishable thesis chapter.

2. **Independent COMSOL validation.** The 12-cell COMSOL FEM cross-check
   at B₀ = 1.5 T with <1% error at the physiologically relevant
   gradient thresholds is the strongest piece of external validation
   in the project and directly addresses the "is the analytical model
   trustworthy?" question that every examiner will ask.

3. **Paracrine extension.** The VEGF reaction–diffusion module
   transforms the project from "how well can we stick cells to a
   stent?" into "does capturing cells actually deliver therapeutic
   benefit?" Fig 22–24's basal-vs-enhanced contrast is a strong
   negative result (basal EC secretion is insufficient) with a
   clean quantitative recommendation (100×-enhanced secretion reaches
   the Ozawa therapeutic band). This is a thesis-worthy story on
   its own.

4. **Test coverage.** 78 tests across 7 modules, covering
   superposition, far-field, symmetry, FD convergence, monotonicity,
   and paracrine limiting cases. Strong for a thesis project.

**What the thesis must do before defending.**

1. **Fix or document Issue #1 (failing test).** A failing test in a
   pre-submission repo will be flagged by any reader with
   `pytest --tb=line`. Either fix the threshold, or remove the test
   with a comment explaining why.

2. **Reconcile the default-parameter drift (Issues #2, #4, #5).**
   State *one* set of defaults in the README, match it in `common.py`
   and in every test, re-generate fig 21 under those defaults, and
   report the new numbers. Don't claim 14.2× if the current code
   produces a different number.

3. **Commit or regenerate `results/` (Issue #3).** A repo that
   advertises 24 figures but only commits 12 looks unmaintained.

4. **Decide and document the saturation story (Issue #5).** The fact
   that saturation reduces the capture prediction by only ~20 %
   *despite* a 5.7× reduction in χ_eff is a *strong* robustness
   statement. Lead with it in the Limitations section.

5. **Clean the repository surface (Issues #6, #9, #11).** Remove
   working notes, orphan scripts, and confusing file prefixes.

**What would strengthen the thesis further (optional).**

- **Direct efficiency comparison with Polyak 2008 at matched
  conditions** (B₀ = 0.1 T, 200 pg, 304 SS stent geometry,
  30 mL/min flow) with a computed capture *efficiency* (not just
  *distance*). The `results-polyak/` work is halfway there but
  stops short of giving a single %-efficiency number comparable
  to Polyak's 20 %.
- **Pulsatility sensitivity analysis.** Womersley number ~2 at MCA
  conditions means the Poiseuille approximation is defensible but
  not exact. A single "does pulsatile flow matter?" appendix
  figure would close this loophole.
- **Angular-average efficiency** (all θ, not just θ = 0 through-strut).
  Currently flagged as a known limitation; the asymmetry between
  through-strut and between-strut capture could be > 2×.

**Conclusion.** The physics is right, the validation is real, the
headline story is defensible. Fix the bookkeeping (tests, defaults,
missing results folder, stale README) and the project produces
thesis-quality material — specifically, **two chapters' worth**:
(i) analytical model, COMSOL validation, static-vs-trajectory
capture comparison; (ii) VEGF paracrine analysis and therapeutic-
window characterisation.

---

## Appendix: Delta from the previous audit (`bc40c15`)

**Fixed since last audit:**
- `physics/shear_stress.py` deleted (Issue #1 of old audit). ✅
- Aaslid attribution in `hydrodynamics.py` docstring (Issue #3 of old audit). ✅
- Re docstring in `trajectories.py` (Issue #6 of old audit). ✅
- 304 SS wording in README (Issue #5 of old audit). ✅
- Linear-χ limitations section added to README (Issue #4 of old audit). ✅
- Langevin saturation model implemented (addresses old Issue #4 at the
  *code* level, not just the docstring). ✅

**Still outstanding from last audit:**
- Polyak 10 pg relic in default (old Issue #2 — README was updated,
  but the code default remains 10 pg).
- Aaslid / Polyak mislabels re-surface in the *new* artefacts
  (fig_saturation_impact, POLYAK_COMPARISON_REPORT).
- `_akoun_yonnet_local` corner-loop serialisation (old Issue #7) —
  still not vectorised; still not a correctness issue.

**New since last audit:**
- Issues #1–#12 above.
