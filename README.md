# stent_capture — Magnetic Stent Cell Capture Modelling

Analytical physics model, FEM-validated, for SPION-labelled endothelial-cell
capture by a magnetisable stent in cerebral arterial flow. Includes a
VEGF paracrine reaction–diffusion extension for therapeutic-window analysis.

- **Stage 1**: 3-D B-field and gradient via the Akoun & Yonnet (1984)
  analytical expressions; external uniform-field superposition.
- **Stage 2**: Magnetic force on SPION-labelled cells (Furlani & Ng 2006)
  with optional Langevin saturation; Poiseuille blood-flow drag; static
  capture criterion; capture maps.
- **Stage 3**: Single-cell trajectory integration (`scipy.solve_ivp` RK45);
  population capture efficiency vs velocity and SPION loading; direct
  static-vs-trajectory comparison (fig 21, headline result).
- **Stage 3c**: Physics audit; VEGF paracrine reaction–diffusion module
  (`paracrine/`) with therapeutic-zone analysis.
- **Stage 3d**: COMSOL FEM cross-validation at B₀ = 1.5 T (calibrated
  effective M), Langevin SPION saturation model, 200 pg / Polyak
  comparison runs.

---

## Headline results

At cerebral-flow conditions (v̄ = 0.2 m/s, B₀ = 0.5 T axial, 50 pg SPION
loading per cell, Langevin saturation model, 12-strut geometry):

- **Trajectory analysis predicts an effective capture range of ~97 µm**
  from the stent inner surface (50 pg loading). The static Furlani & Ng
  force-balance criterion predicts **zero capture** under the same
  conditions, demonstrating that trajectory integration is
  **essential for quantitative prediction** — static-only analysis would miss
  all capture at physiological flow rates.

- **Capture range vs loading** (v̄ = 0.2 m/s, B₀ = 0.5 T): from ~51 µm
  at 10 pg to ~131 µm at 200 pg. Static capture remains negligible
  (< 37 µm) across the entire loading range, while trajectory predictions
  increase monotonically with loading.

- **Near-wall capture efficiency** (injection band r = 1.20–1.45 mm,
  20 cells): ~45 % at v̄ = 0.2 m/s and 50 pg loading (the updated
  headline dose); ~50 % at higher 200 pg loading. Monotonically falling
  with velocity, monotonically rising with loading.

- **Higher flow velocities increase the extension factor**: at
  v̄ = 0.5 m/s the static criterion predicts 0 µm yet trajectory
  integration finds a ~74 µm capture range at 50 pg. The upstream
  radial drift accumulated over the 2 mm approach trajectory
  dominates.

- **COMSOL FEM cross-validation** at B₀ = 1.5 T: analytical gradient
  matches COMSOL V2-2C (12-cell geometry, µᵣ = 2 soft ferromagnet) to
  within **<1 % at 100 T/m and 40 T/m** threshold crossings, 12 % at
  300 T/m (very near the strut surface).

- **Langevin saturation sensitivity**: replacing the constant-χ SPION
  force model with the full Langevin `χ_eff(B)` reduces capture
  predictions by only ~20–25 % despite a 5.7× reduction in χ at
  B₀ = 1.5 T (because the stent's local field in the active capture
  zone already forces both models into a near-saturated regime).
  **Results are robust to the SPION saturation assumption.**

- **VEGF paracrine analysis** (fig 22–24, 27–29): basal EC secretion
  (0.068 molecules / cell / s, Yen 2011) produces sub-ng/mL tissue
  concentrations even with hundreds of captured cells. When constrained
  by Ozawa et al. (2004) population secretion categories (normal:
  5–70 ng/10⁶ cells/day; threshold: 70–100; aberrant: >100), the model
  shows that **~320 cells at normal secretion** or **~160 cells at
  threshold secretion** are needed to reach tissue VEGF concentrations
  of ~5 ng/mL. Fig 28/29 provide the parametric sweep across cell count
  and secretion rate. Note: 5 ng/mL is an engineering reference rather
  than a literature-validated tissue threshold — Ozawa measured
  production rates per cell population, not tissue microenvironmental
  concentrations. The module's value is as a tool for predicting
  tissue VEGF distribution rather than a binary therapeutic claim.

The static-vs-trajectory contrast is consistent with the experimental
capture efficiencies reported by Polyak et al. (2008) and Tefft et al.
(2014, 2017).

---

## Package layout

```
stent_capture_project/
├── legacy_stent_analysis.py     # original monolithic script (reference only)
├── README.md
├── CHANGELOG.md
├── AUDIT.md                     # physics & code audit (current)
├── RESULTS_SUMMARY.md           # results folder index & comparison
├── SATURATION_MODEL_SUMMARY.md  # Langevin vs constant-χ analysis
├── Overview.pdf                 # 25-page technical overview
├── results-200pg/               # 200 pg SPION variant (Polyak working dose)
├── results-polyak/              # Polyak B₀ = 0.1 T comparison
├── scripts/                     # analysis + regeneration scripts
└── stent_capture/
    ├── core/
    │   ├── field_model.py       # StentRing — 3-D Akoun & Yonnet model
    │   └── gradient.py          # compute_gradient_magnitude/vector (FD)
    ├── physics/
    │   ├── external_field.py    # UniformExternalField + TotalField
    │   ├── magnetic_force.py    # SPIONLabelledCell, Langevin χ_eff, magnetic_force()
    │   ├── hydrodynamics.py     # BloodFlow (Poiseuille) + stokes_drag()
    │   └── capture_criterion.py # capture_map(), capture_distance()
    ├── simulation/
    │   ├── trajectories.py      # integrate_trajectory (RK45 solve_ivp)
    │   └── capture_efficiency.py# sweep_injection_line, v- and loading-sweeps
    ├── paracrine/
    │   ├── transport.py         # ParacrineField — 2-D reaction–diffusion solver
    │   ├── secretion.py         # VEGFSource — Gaussian source per captured cell
    │   └── therapeutic.py       # therapeutic-zone radius, time-to-threshold
    ├── figures/
    │   ├── style.py             # dissertation rcParams
    │   ├── common.py            # DEFAULTS, COMSOL calibration, make_ring()
    │   ├── fig01 .. fig11       # Stage 1: field, gradient, geometry sweeps
    │   ├── fig12                # Stage 1: external-field comparison
    │   ├── fig13 .. fig16       # Stage 2: force, drag, capture maps
    │   ├── fig17 .. fig21       # Stage 3: trajectories, efficiency, headline
    │   ├── fig22 .. fig24       # Stage 3c: VEGF paracrine
    │   ├── fig25_comsol_gradient_validation # Stage 3d: 3-panel COMSOL validation
    │   ├── fig26_comsol_multigeometry       # Stage 3d: 4-panel multi-geometry analysis
    │   └── fig_saturation_impact            # Stage 3d: Langevin vs constant-χ comparison
    └── tests/
        ├── test_external_field.py     # 19 tests
        ├── test_magnetic_force.py     # 19 tests (inc. Langevin saturation)
        ├── test_hydrodynamics.py      # 12 tests
        ├── test_capture_criterion.py  # 11 tests
        ├── test_trajectories.py       #  5 tests
        ├── test_capture_efficiency.py #  5 tests
        └── test_paracrine.py          #  4 tests
```

---

## Quick start

### Requirements

```
numpy
scipy
matplotlib
pytest  (for tests)
```

### Run tests

```bash
python -m pytest stent_capture/tests/ -v
```

Expected: **78 / 78 passing** on Python 3.14 / numpy 2.4 / scipy 1.17.
All tests include Langevin saturation model validation and COMSOL cross-checks.

### Regenerate all figures

From the project root:

```bash
# All 26 figures (Stage 1–3c + COMSOL validation) with thesis defaults
python -m scripts.regenerate_original_results

# 200 pg SPION variant (Polyak working dose)
python -m scripts.generate_complete_200pg_results

# Polyak B₀ = 0.1 T comparison
python -m scripts.generate_polyak_comparison

# COMSOL calibration + saturation sensitivity analysis
python -m scripts.check_comsol_calibration
python -m scripts.compare_saturation_models
```

Output goes to `results/` (main figures 1–26), `results-200pg/` (variant),
and `results-polyak/` (comparison) respectively.

### Run a single figure

```bash
python -m stent_capture.figures.fig21_static_vs_trajectory
python -m stent_capture.figures.fig22_concentration_field
python -m stent_capture.figures.fig25_comsol_gradient_validation
python -m stent_capture.figures.fig26_comsol_multigeometry
python -m stent_capture.figures.fig_saturation_impact
```

---

## Physics summary

### Field model

Each stent strut is a uniformly magnetised rectangular prism. The B-field
is computed via the Akoun & Yonnet (1984) closed-form expressions — the
same kernel implemented in *magpylib*. The formula evaluates corner sums
of `arctan` and `ln` terms derived by integrating the scalar potential of
the magnetic surface charges over each face.

Global frame: stent axis along z, ring in the z = 0 plane. Each strut has
a local frame aligned with its radial, circumferential, and axial
directions.

### External field

A spatially uniform external field B₀ is modelled by
`UniformExternalField`. The total field is:

```
B_total(r) = B_stent(r) + B0
```

The gradient is computed on B_total (not B_stent alone):

```
|∇|B_total||  ≠  |∇|B_stent||   (when B₀ ≠ 0)
```

because the gradient of a magnitude is non-linear:
`∂|B_total|/∂x = (B_total · ∂B_stent/∂x) / |B_total|`.

### Stent saturation assumption

For ferromagnetic stents at applied fields B₀ ≳ 0.1 T the material
approaches saturation. Setting `assume_saturation=True` on `StentRing`
signals that M represents the saturation magnetisation and should be
used regardless of B₀ magnitude. Typical values:

| Material                       | M_sat (MA/m) |
|--------------------------------|--------------|
| Cold-worked / ferritic SS (304, 430) | 1.0–1.2 |
| Pure Fe                        | 1.7          |

This flag does **not** implement an M(H) hysteresis curve. For
quantitative agreement with COMSOL (µᵣ = 2 soft ferromagnet at
B₀ = 1.5 T), use the calibrated effective magnetisation:

```python
from stent_capture.figures.common import make_ring
ring = make_ring(B0_magnitude=1.5)   # auto-swaps to M_COMSOL_EFF_B15 = 2.20 MA/m
```

At B₀ = 1.5 T this produces gradient threshold crossings within <1 %
of COMSOL at 100 T/m and 40 T/m (see fig04_comsol_gradient_validation).

### Capture force (Stage 2)

The force on a superparamagnetic cell is the scalar-gradient form of the
Furlani & Ng (2006) dipole force:

```
F = (V_spion · χ_eff(|B|) / μ₀) · |B_total| · ∇|B_total|
```

Two susceptibility models are available:

1. **Constant-χ (linear)**: `χ_eff(B) = χ₀ = 2.0`. Activated with
   `SPIONLabelledCell(spion_sat_magnetization=None)`. Valid in the
   low-field limit `B ≪ M_sat μ₀ / χ₀ ≈ 280 mT`.

2. **Langevin saturation (default)**: `χ_eff(B) = χ₀ · 3 L(ξ) / ξ`
   where `ξ = 3 χ₀ B / (M_sat μ₀)` and `L(x) = coth(x) − 1/x`.
   Uses `M_sat = 446 kA/m` for bulk magnetite (Furlani & Ng 2006
   Table 1). Activated by default in `SPIONLabelledCell()` via
   `spion_sat_magnetization=446e3`. Taylor-expanded near ξ = 0 for
   numerical stability.

Key insight: `|∇|B_total||` is reduced by an axial B₀ (B_total rotates
perpendicular to the stent field, suppressing the projection) but
`|B_total|` is increased ~17× by B₀ = 0.5 T. The product — the
*force parameter* — increases 3–55× depending on distance.

**Constant-χ vs Langevin sensitivity** (SATURATION_MODEL_SUMMARY.md):
switching to Langevin reduces capture predictions by ~20–25 % despite a
5.7× reduction in χ at B₀ = 1.5 T, because the stent's local field
already forces both models into a saturated regime.

### Blood flow drag

Poiseuille flow in a 1.54 mm radius vessel (cerebral arterial / M1
MCA-representative). Stokes drag:

```
F_drag = 6 π η R_cell (v_blood − v_cell)
```

### Vessel geometry and flow parameters

The model represents a cerebral vessel with mean velocity 0.2 m/s
(representative of distal/diseased-vessel flow; healthy MCA mean per
Aaslid et al. 1982 is ~0.62 m/s). Vessel diameter ~3 mm. Parameters
can be adjusted for smaller distal vessels (M2/M3, ~1.5–2 mm diameter,
0.1–0.15 m/s) where flow-diverter stents are commonly deployed. The
velocity sweep range 0.02–0.5 m/s covers distal small vessels through
MCA peak systolic flow.

### Static capture criterion

A cell at r is captured if `|F_mag(r)| > |F_drag(r)|` (conservative
scalar Furlani & Ng 2006 criterion). The `direction='inward'` mode
of `capture_distance()` sweeps from the stent inner surface toward the
vessel centre — the physically correct lumen geometry.

### Trajectory integration (Stage 3)

Terminal-velocity ODE `dr/dt = v_blood(r) + F_mag(r) / (6π η R_cell)`
integrated with `scipy.integrate.solve_ivp(method='RK45',
rtol=1e-6, atol=1e-9)`. Three event-based termination conditions:
escape (`z ≥ z_end`), strut proximity (cylindrical approximation),
wall contact. Cell Reynolds number Re ≈ 1.1 at v = 0.2 m/s — Stokes
regime, inertia negligible, first-order ODE is physically
appropriate.

### VEGF paracrine (Stage 3c)

2-D reaction–diffusion-advection on an unrolled tissue slab around the
vessel wall:

```
∂C/∂t = D ∇²C − u·∇C + S(x, z) − k C
```

- `D = 1.04 × 10⁻¹⁰ m²/s` — VEGF₁₆₅ in tissue ECM, in vivo skeletal
  muscle (Mac Gabhann & Popel 2006).
- `k = 1.93 × 10⁻⁴ s⁻¹` — first-order degradation, half-life ≈ 60 min
  (in-vitro; Kleinheinz et al. 2010, Chen et al.; via PLOS ONE 2011).
- Characteristic diffusion length L_D = √(D/k) ≈ 734 µm.
- Source: Gaussian kernel per captured cell at q_cell = 0.068 molecules/s
  (Yen et al. 2011 mouse VEGF model); can be scaled to match Ozawa
  population secretion categories.
- Optional axial advection (tissue interstitial flow ~1–100 µm/s);
  Péclet number Pe = u·L_D/D ranges from ~7 (normal) to ~350
  (inflamed).

Both steady-state (sparse Kronecker construction + `scipy.sparse.linalg.
spsolve`) and transient (explicit Euler with CFL-limited step) solvers
are provided. Zero-flux in axial z, periodic in circumferential x.

**Ozawa-constrained framing (fig 28/29):** rather than the arbitrary
"100× basal" multiplier used in early versions, the secretion rate is
parametrised by Ozawa et al. (2004) cell-population categories:

| Category  | Rate (ng/10⁶ cells/day) | Biology                     |
|-----------|-------------------------|-----------------------------|
| Normal    | 5–70                    | Normal angiogenesis         |
| Threshold | 70–100                  | Angiogenic threshold region |
| Aberrant  | >100                    | Aberrant angiogenesis       |

Note: Ozawa measures *cell production rates*, not *tissue VEGF
concentrations*. The 5–25 ng/mL "therapeutic" band sometimes cited
elsewhere is an engineering reference only.

---

## Known limitations

1. **SPION saturation assumption.** Default is the Langevin model; the
   constant-χ linear approximation is available for comparison. At
   B₀ = 0.5 T the two differ by ~20–25 % in absolute force magnitude
   but < 5 % in the comparative (static vs trajectory) conclusions,
   because both models are saturated in the near-strut field.

2. **Strut geometry approximation.** The Akoun & Yonnet kernel models
   struts as uniformly magnetised rectangular prisms. Real struts have
   rounded edges; COMSOL validation at B₀ = 1.5 T shows <1% gradient
   agreement at 40–100 T/m and ~12 % at 300 T/m (near-surface). The
   cylindrical strut-proximity metric in trajectory integration
   over-estimates the circumferential half-width by 20–30 µm — a
   known limitation retained from Stage 3a.

3. **Flow velocity attribution.** The default mean velocity 0.2 m/s is
   representative of distal / diseased-vessel flow. Healthy MCA
   (Aaslid et al. 1982) is ~0.62 m/s. Adjust `mean_velocity` for your
   vessel target. Pulsatility is not modelled (Womersley number ~2 at
   MCA conditions — quasi-steady Poiseuille is defensible but not
   exact).

4. **Angular-average efficiency.** All trajectories are injected at
   θ = 0 (strut-aligned). Real vessel wall has `n_struts` strut
   locations; efficiency will vary by angular position. The
   through-strut vs between-strut asymmetry could be > 2× and is not
   sampled.

5. **Paracrine advection.** The VEGF solver supports optional axial
   advection. Péclet number Pe = u·L_D/D (L_D ≈ 734 µm) ranges from
   ~7 (normal tissue, u ~1 µm/s) to ~350 (inflamed tissue, u ~50 µm/s),
   so advection is secondary under normal conditions but can dominate
   in pathological flow regimes. Fig 27 explores this sensitivity.

6. **Paracrine secretion rate.** q_cell = 0.068 molecules/cell/s from
   Yen et al. (2011) mouse VEGF model — originally per myonuclear
   domain, repurposed here per single endothelial cell. The preferred
   approach is to parameterise secretion via Ozawa population
   categories (normal/threshold/aberrant) rather than an arbitrary
   basal multiplier; fig 28/29 use this framing.

7. **"Therapeutic window" framing.** The 5–25 ng/mL band referenced in
   early figures is an engineering approximation, not a validated
   tissue threshold. Ozawa et al. (2004) reports population secretion
   rates (ng/10⁶ cells/day), not tissue microenvironmental
   concentrations. The module should be interpreted as a predictor of
   tissue VEGF distribution for a given cell density and secretion
   rate — the therapeutic significance of a given concentration
   depends on downstream cell-response modelling outside scope.

---

## Default parameters

### Stent geometry (`figures/common.py` DEFAULTS)

| Parameter  | Value      | Description                               |
|------------|------------|-------------------------------------------|
| R          | 1.5 mm     | Ring radius                               |
| w          | 100 µm     | Strut circumferential width               |
| t          | 80 µm      | Strut radial thickness                    |
| L          | 500 µm     | Strut axial length                        |
| M          | 1.0 MA/m   | Magnetisation (cold-worked/ferritic SS)   |
| n_struts   | **12**     | V2-2C COMSOL-matched geometry             |
| mag_mode   | radial     | Magnetisation direction                   |

Note: test fixtures and the historical headline fig 21 run use
`n_struts = 8`. See `AUDIT.md` Issue #4 for the geometry-drift
reconciliation.

### Cell / flow parameters

| Parameter                 | Value          | Description                                        |
|---------------------------|----------------|----------------------------------------------------|
| Cell radius               | 10 µm          | Typical endothelial cell                           |
| SPION mass (code default) | 10 pg          | Lower-loading reference (Polyak 2008 used 200 pg)  |
| χ₀ (low-field SPION)      | 2.0            | SI volume susceptibility                           |
| M_sat (SPION)             | 446 kA/m       | Bulk magnetite (Furlani & Ng 2006 Tbl 1)           |
| ρ_spion                   | 5170 kg/m³     | Magnetite density                                  |
| vessel_radius             | 1.54 mm        | Cerebral artery / stent outer surface              |
| mean_velocity             | 0.2 m/s        | Distal/diseased MCA (healthy ~0.62 m/s per Aaslid) |
| blood viscosity           | 4 mPa·s        | Whole blood (40 % Hct)                             |
| blood density             | 1060 kg/m³     | Pedley 1980                                        |
| B₀ (default)              | 0.5 T          | Axial external field                               |

### COMSOL-calibrated magnetisation (`figures/common.py`)

| Parameter         | Value      | Calibrated for                                    |
|-------------------|------------|---------------------------------------------------|
| M_COMSOL_EFF      | 0.619 MA/m | B₀ = 0, matches COMSOL 100 T/m crossing           |
| M_COMSOL_EFF_B15  | 2.20 MA/m  | B₀ = 1.5 T (MRI), <1% error at 100/40 T/m         |

Call `make_ring(B0_magnitude=1.5)` to auto-swap in the 1.5 T
calibration.

### VEGF paracrine defaults (`paracrine/`)

| Parameter | Value | Source |
|-----------|-------|--------|
| D_VEGF_TISSUE | 1.04 × 10⁻¹⁰ m²/s | Mac Gabhann & Popel 2006 (VEGF₁₆₄ in vivo) |
| K_DEG_TISSUE | 1.93 × 10⁻⁴ s⁻¹ | In vitro (Kleinheinz 2010, Chen et al.; via PLOS ONE 2011) |
| L_diffusion | ≈ 734 µm | √(D/k) |
| q_cell | 0.068 molecules/cell/s | Yen et al. 2011 (mouse VEGF model) |
| σ (Gaussian source) | 10 µm | cell radius (engineering choice) |
| h (tissue slab) | 20 µm | ~2 cell layers (engineering choice) |
| C_THERAPEUTIC_LOW | 5 ng/mL | Engineering reference (not directly from literature) |
| C_THERAPEUTIC_HIGH | 25 ng/mL | Engineering reference (order-of-magnitude) |
| C_ABERRANT | 100 ng/mL | Related to Ozawa threshold (units differ) |

Ozawa population secretion categories (ng/10⁶ cells/day, Ozawa et al.
2004): normal 5–70, threshold 70–100, aberrant >100.

---

## Stage roadmap

| Stage | Content                                                       | Status       |
|-------|---------------------------------------------------------------|--------------|
| 1     | Package refactor + uniform external field + force parameter   | **Complete** |
| 2     | Magnetic force (pN); Poiseuille drag; static capture maps     | **Complete** |
| 3     | Cell trajectory ODE integration; capture efficiency           | **Complete** |
| 3c    | Physics audit; VEGF paracrine signalling module               | **Complete** |
| 3d    | COMSOL FEM validation; Langevin saturation; 200pg / Polyak    | **Complete** |

See `AUDIT.md` for the per-module physics audit and an itemised list of
known issues to address before thesis submission.

---

## Key references

- **Akoun, G. & Yonnet, J.-P. (1984).** 3D analytical calculation of the
  forces exerted between two cuboidal magnets. *IEEE Trans. Magn.* 20, 1962.
- **Furlani, E.P. & Ng, K.C. (2006).** Analytical model of magnetic
  nanoparticle transport and capture in the microvasculature.
  *Phys. Rev. E* 73, 061919.
- **Polyak, B. et al. (2008).** High field gradient targeting of magnetic
  nanoparticle-loaded endothelial cells to the surfaces of steel stents.
  *PNAS* 105(2), 698–703.
- **Tefft, B.J. et al. (2014).** Magnetizable stent-grafts enable
  endothelial cell capture. *IEEE Trans. Magn.* 50(11), 1–4.
- **Chorny, M. et al. (2007).** Targeting stents with locally delivered
  paclitaxel-loaded magnetic nanoparticles. *FASEB J.* 21, 2510.
- **Mac Gabhann, F. & Popel, A.S. (2006).** Interactions of VEGF
  isoforms. *PLoS Comput. Biol.* 2(10), e127.
- **Stefanini, M.O. et al. (2008).** Molecular-level simulations of VEGF
  distributions. *PLoS ONE* 3(11), e3565.
- **Ozawa, C.R. et al. (2004).** Microenvironmental VEGF concentration,
  not total dose, determines a threshold behavior of angiogenesis.
  *J. Clin. Invest.* 113(4), 516–27.
- **Aaslid, R. et al. (1982).** Noninvasive transcranial Doppler
  ultrasound recording of flow velocity. *J. Neurosurg.* 57, 769.
