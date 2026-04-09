# stent_capture — Magnetic Stent Cell Capture Modelling

Stages 1–2 of a 4-stage build-out for modelling magnetic SPION-labelled
endothelial cell capture in vascular stents.

- **Stage 1**: 3-D B-field and gradient magnitude using the Akoun & Yonnet
  (1984) analytical expressions; external uniform field superposition.
- **Stage 2**: Magnetic force on SPION-labelled cells; Poiseuille blood flow
  drag; static capture criterion |F_mag| > |F_drag|; capture maps.

Roadmap: Stage 3 = trajectory integration, Stage 4 = geometry optimisation.

---

## Package layout

```
stent_capture_project/
├── legacy_stent_analysis.py     # original monolithic script (reference only)
├── README.md
├── CHANGELOG.md
└── stent_capture/
    ├── core/
    │   ├── field_model.py       # StentRing — 3-D Akoun & Yonnet model
    │   └── gradient.py          # compute_gradient_magnitude/vector (FD)
    ├── physics/
    │   ├── external_field.py    # UniformExternalField + TotalField
    │   ├── magnetic_force.py    # SPIONLabelledCell + magnetic_force()
    │   ├── hydrodynamics.py     # BloodFlow (Poiseuille) + stokes_drag()
    │   └── capture_criterion.py # capture_map() + capture_distance()
    ├── figures/
    │   ├── style.py             # shared rcParams (dissertation quality)
    │   ├── common.py            # DEFAULTS, THRESHOLDS, make_ring()
    │   ├── fig01_single_strut.py  .. fig11_rz_heatmap.py   # Stage 1
    │   ├── fig12_external_field_comparison.py               # Stage 1
    │   ├── fig13_force_parameter.py                         # Stage 1
    │   ├── fig14_force_vs_distance.py                       # Stage 2
    │   ├── fig15_drag_vs_velocity.py                        # Stage 2
    │   └── fig16_capture_map.py                             # Stage 2
    └── tests/
        ├── test_external_field.py    # 19 tests — Stage 1
        ├── test_magnetic_force.py    # 15 tests — Stage 2
        ├── test_hydrodynamics.py     # 12 tests — Stage 2
        └── test_capture_criterion.py #  9 tests — Stage 2
```

---

## Quick start

### Requirements

```
numpy
matplotlib
pytest  (for tests only)
```

### Run all figures

From the project root:

```bash
# All 11 original figures
python -c "
from stent_capture.figures import (
    fig01_single_strut, fig02_ring_heatmaps, fig03_gradient_vs_distance,
    fig04_magnetisation_sweep, fig05_strut_dimensions, fig06_n_struts,
    fig07_gradient_contours, fig08_force_parameter, fig09_axial_profile,
    fig10_convergence, fig11_rz_heatmap,
)
for mod in [fig01_single_strut, fig02_ring_heatmaps, fig03_gradient_vs_distance,
            fig04_magnetisation_sweep, fig05_strut_dimensions, fig06_n_struts,
            fig07_gradient_contours, fig08_force_parameter, fig09_axial_profile,
            fig10_convergence, fig11_rz_heatmap]:
    mod.main()
"

# New Stage 1 external-field comparison figure
python -m stent_capture.figures.fig12_external_field_comparison
```

Output goes to `results/` in the project root (PNG + PDF).

### Run a single figure

```bash
python -m stent_capture.figures.fig03_gradient_vs_distance
```

### Run tests

```bash
python -m pytest stent_capture/tests/ -v
```

---

## Physics summary

### Field model

Each stent strut is a uniformly magnetised rectangular prism.  The B-field is
computed using the Akoun & Yonnet (1984) closed-form expressions — the same
kernel implemented in *magpylib*.  The formula evaluates corner sums of arctan
and ln terms derived from integrating the scalar potential of the magnetic
surface charges over each face.

Global frame: stent axis along z, ring in the z = 0 plane.  Each strut has a
local frame aligned with its radial, circumferential, and axial directions.

### External field

A spatially uniform external field B0 is modelled by `UniformExternalField`.
The total field is:

```
B_total(r) = B_stent(r) + B0
```

The gradient is computed on B_total (not B_stent alone):

```
|∇|B_total|| ≠ |∇|B_stent||  (when B0 ≠ 0)
```

because the gradient of a magnitude is non-linear:
`∂|B_total|/∂x = (B_total · ∂B_stent/∂x) / |B_total|`

### Saturation assumption

For ferromagnetic stents at applied fields B0 ≳ 0.1 T, the material approaches
saturation.  Setting `assume_saturation=True` on `StentRing` signals that M
represents the saturation magnetisation and should be used directly regardless
of B0 magnitude.  Typical values:

| Material | M_sat (MA/m) |
|----------|-------------|
| 304 SS   | ~1.0        |
| 430 SS   | ~1.2        |
| Pure Fe  | ~1.7        |

This flag does **not** implement an M(H) hysteresis curve — that is deferred.

### Capture force (Stage 2)

The force on a superparamagnetic cell is (Furlani & Ng 2006):

```
F = (V_spion * chi_spion / mu_0) * |B_total| * nabla|B_total|
```

Key insight: `|nabla|B_total||` is **reduced** by an axial B0 (B_total rotates
perpendicular to the stent field, suppressing the projection).  But
`|B_total|` is increased ~17× by B0 = 0.5 T.  The product — the **force
parameter** — increases 3–55× depending on distance.  Use `B_magnitude * grad_B`
not `grad_B` alone to assess capture benefit.

### Blood flow drag

Poiseuille flow in a 1.54 mm radius vessel (cerebral arterial / M1-segment
MCA-representative); Stokes drag:

```
F_drag = 6 * pi * eta * R_cell * v_blood
```

### Vessel geometry and flow parameters

The model represents an M1-segment cerebral artery (vessel diameter ~3 mm,
mean velocity 0.2 m/s; Aaslid et al. 1982).  Parameters can be adjusted for
smaller distal vessels (M2/M3, ~1.5–2 mm diameter, 0.1–0.15 m/s) where
flow-diverter stents are commonly deployed.

The velocity range 0.05–0.5 m/s covers distal small vessels through MCA peak
systolic flow.

### Capture criterion

A cell is captured at position r if:

```
|F_mag(r)| > |F_drag(r)|
```

This is the conservative Furlani & Ng (2006) scalar criterion.  Directional
analysis (Stage 3 trajectories) will extend the capture zone.

---

## Default stent parameters

| Parameter | Value  | Description                    |
|-----------|--------|--------------------------------|
| R         | 1.5 mm | Ring radius                    |
| w         | 100 µm | Strut circumferential width    |
| t         | 80 µm  | Strut radial thickness         |
| L         | 500 µm | Strut axial length             |
| M         | 1 MA/m | Magnetisation (~304 SS sat.)   |
| n_struts  | 8      | Number of struts               |
| mag_mode  | radial | Magnetisation direction        |

---

## Stage roadmap

| Stage | Content                                                      | Status       |
|-------|--------------------------------------------------------------|--------------|
| 1     | Package refactor + uniform external field + force parameter  | **Complete** |
| 2     | Magnetic force (pN); Poiseuille drag; static capture maps    | **Complete** |
| 3     | Cell trajectory ODE integration; capture efficiency          | Planned      |
| 4     | Geometry optimisation (strut shape, n_struts, R)             | Planned      |

### Default cell / flow parameters (Stage 2)

| Parameter            | Value          | Description                           |
|----------------------|----------------|---------------------------------------|
| Cell radius          | 10 µm          | Typical endothelial cell              |
| SPION mass           | 10 pg          | Iron oxide per cell (Polyak 2008)     |
| chi_spion            | 2.0            | SPION material susceptibility         |
| rho_spion            | 5170 kg/m³     | Magnetite density                     |
| vessel_radius        | 1.54 mm        | Cerebral artery (MCA) / stent outer surface |
| mean_velocity        | 0.2 m/s        | MCA mean flow (Aaslid et al. 1982)    |
| blood_viscosity      | 4 mPa·s        | Whole blood                           |
