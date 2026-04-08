# stent_capture — Magnetic Stent Field Gradient Analysis

Stage 1 of a 4-stage build-out for modelling magnetic cell capture in
vascular stents.  Computes the 3-D B-field and gradient magnitude around
magnetised stent struts using the Akoun & Yonnet (1984) / Furlani (2001)
closed-form analytical expressions for uniformly magnetised rectangular
prisms.

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
    │   └── gradient.py          # compute_gradient_magnitude (FD, any field)
    ├── physics/
    │   └── external_field.py    # UniformExternalField + TotalField composer
    ├── figures/
    │   ├── style.py             # shared rcParams (dissertation quality)
    │   ├── common.py            # DEFAULTS, THRESHOLDS, make_ring(), save_fig()
    │   ├── fig01_single_strut.py
    │   ├── fig02_ring_heatmaps.py
    │   ├── fig03_gradient_vs_distance.py
    │   ├── fig04_magnetisation_sweep.py
    │   ├── fig05_strut_dimensions.py
    │   ├── fig06_n_struts.py
    │   ├── fig07_gradient_contours.py
    │   ├── fig08_force_parameter.py
    │   ├── fig09_axial_profile.py
    │   ├── fig10_convergence.py
    │   ├── fig11_rz_heatmap.py
    │   └── fig12_external_field_comparison.py
    └── tests/
        └── test_external_field.py
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

This flag does **not** implement an M(H) hysteresis curve — that is Stage 2.

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

| Stage | Content                                              | Status      |
|-------|------------------------------------------------------|-------------|
| 1     | Package refactor + uniform external field            | **Complete** |
| 2     | Magnetic force on SPION-labelled cells; M(H) curve  | Planned     |
| 3     | Cell trajectory integration; capture efficiency     | Planned     |
| 4     | Geometry optimisation (strut shape, n_struts, R)    | Planned     |
