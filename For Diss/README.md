# Dissertation Appendix — Magnetic Stent Cell-Capture Code

Standalone code package accompanying the dissertation. This folder contains
everything needed to regenerate the eight figures cited in the main text:

| Dissertation | Script (package path)                                            |
|--------------|------------------------------------------------------------------|
| Fig 7        | `stent_capture/figures/fig07_comsol_gradient_validation.py`      |
| Fig 7b       | `stent_capture/figures/fig07_comsol_gradient_validation.py`      |
| Fig 7 extra  | `stent_capture/figures/fig07extra_comsol_multigeometry.py`       |
| Fig 8        | `stent_capture/figures/fig08_force_parameter.py`                 |
| Fig 9        | `stent_capture/figures/fig09_spion_loading_sweep.py`             |
| Fig 10       | `stent_capture/figures/fig10_single_trajectory.py`               |
| Fig 11       | `stent_capture/figures/fig11_static_vs_trajectory.py`            |
| Fig 12       | `stent_capture/figures/fig12_capture_efficiency.py`              |

Fig 7 and Fig 7b are produced by a single script (bar chart of threshold
crossings + log-log gradient profile overlay).

---

## Requirements

- Python ≥ 3.10
- See `requirements.txt`:
  - numpy, scipy, matplotlib, openpyxl (for the COMSOL Excel reader), pytest.

Install into a fresh virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

---

## Regenerating the figures

From inside this folder:

```bash
# All figures at once → For Diss/results/
python run_all_figures.py

# Or a single figure:
python -m stent_capture.figures.fig07_comsol_gradient_validation
python -m stent_capture.figures.fig07extra_comsol_multigeometry
python -m stent_capture.figures.fig08_force_parameter
python -m stent_capture.figures.fig09_spion_loading_sweep
python -m stent_capture.figures.fig10_single_trajectory
python -m stent_capture.figures.fig11_static_vs_trajectory
python -m stent_capture.figures.fig12_capture_efficiency
```

Every script saves PNG (`dpi=150–300`) and PDF copies to `results/`.
Fig 7 (extra) additionally writes two CSV tables
(`fig7extra_powerlawfit_summary.csv`, `fig7extra_linearity_check.csv`).

Runtimes on a laptop-class CPU:

| Figure | Approximate runtime |
|--------|---------------------|
| Fig 7 / 7b      | < 10 s |
| Fig 7 extra     | < 10 s |
| Fig 8           | ~ 20 s |
| Fig 9           | < 10 s |
| Fig 10          | ~ 30 s |
| Fig 11          | ~ 4 min (binary search × two sweeps) |
| Fig 12          | ~ 6 min (velocity × loading sweeps, multiprocessing) |

Fig 11 and Fig 12 integrate many individual cell trajectories with
`scipy.solve_ivp` and use `multiprocessing.Pool` when available.

---

## Running the unit tests

```bash
python -m pytest stent_capture/tests/ -v
```

Tests cover every physics / simulation module used by the figures:

| Test file                      | Module under test                 | Count |
|--------------------------------|-----------------------------------|-------|
| `test_external_field.py`       | core.field_model, external_field  | 19    |
| `test_magnetic_force.py`       | physics.magnetic_force            | 19    |
| `test_hydrodynamics.py`        | physics.hydrodynamics             | 12    |
| `test_capture_criterion.py`    | physics.capture_criterion         | 11    |
| `test_trajectories.py`         | simulation.trajectories           |  5    |
| `test_capture_efficiency.py`   | simulation.capture_efficiency     |  5    |

All tests should pass on Python ≥ 3.10 / numpy ≥ 1.24 / scipy ≥ 1.10.

---

## Folder layout

```
For Diss/
├── README.md
├── requirements.txt
├── run_all_figures.py               # regenerates every figure into results/
├── data/
│   ├── comsol_gradient_data_cleaned.xlsx   # V1–V4 + 2-D linearity datasets
│   └── V1-Mgrad-Plot-New2.csv              # COMSOL V1 cut-line export
├── results/                         # created on first run
└── stent_capture/                   # Python package
    ├── __init__.py
    ├── core/
    │   ├── field_model.py           # StentRing (Akoun & Yonnet 1984 kernel)
    │   └── gradient.py              # finite-difference |∇|B||
    ├── physics/
    │   ├── external_field.py        # UniformExternalField + TotalField
    │   ├── magnetic_force.py        # SPIONLabelledCell, Langevin χ_eff
    │   ├── hydrodynamics.py         # BloodFlow (Poiseuille) + Stokes drag
    │   └── capture_criterion.py     # static |F_mag| > |F_drag|
    ├── simulation/
    │   ├── trajectories.py          # integrate_trajectory (RK45)
    │   └── capture_efficiency.py    # injection-line & parameter sweeps
    ├── data/
    │   └── comsol_loader.py         # openpyxl-based COMSOL sheet reader
    ├── figures/
    │   ├── common.py                # DEFAULTS, COMSOL calibration, make_ring()
    │   ├── style.py                 # dissertation matplotlib rcParams
    │   ├── fig07_comsol_gradient_validation.py
    │   ├── fig07extra_comsol_multigeometry.py
    │   ├── fig08_force_parameter.py
    │   ├── fig09_spion_loading_sweep.py
    │   ├── fig10_single_trajectory.py
    │   ├── fig11_static_vs_trajectory.py
    │   └── fig12_capture_efficiency.py
    └── tests/                       # pytest suite (71 tests)
```

---

## Physics outline

- **Field.** Each stent strut is a uniformly magnetised rectangular prism;
  the B-field is the closed-form Akoun & Yonnet (1984) kernel. The ring
  field is summed over `n_struts` rotated copies. An optional
  `UniformExternalField` adds a spatially uniform `B₀`; `TotalField`
  composes them.
- **Gradient.** `|∇|B||` is obtained via 3-D central finite differences
  (`core.gradient`). Because `|B_total| ≠ |B_stent| + |B₀|` in general,
  the gradient of the magnitude must be computed on the combined field.
- **Force.** Magnetic force on a SPION-labelled cell uses the
  Furlani & Ng (2006) scalar-gradient form
  `F = (V_p χ_eff / μ₀) · |B| · ∇|B|`. A Langevin saturation model sets
  `χ_eff(B) = χ₀ · 3 L(ξ)/ξ` with `ξ = 3 χ₀ B / (M_sat μ₀)`.
- **Drag.** Poiseuille velocity profile in a 1.54 mm cerebral-artery
  lumen; Stokes drag `6π η R_cell (v_blood − v_cell)`.
- **Capture (static).** `|F_mag(r)| > |F_drag(r)|` evaluated along the
  through-strut axis.
- **Capture (trajectory).** Terminal-velocity ODE
  `dr/dt = v_blood + F_mag / (6π η R_cell)` integrated with
  `scipy.integrate.solve_ivp` (RK45). Capture = strut-proximity or
  lumen-wall contact event.

---

## Default parameters (`stent_capture/figures/common.py`)

| Parameter   | Value      | Description                               |
|-------------|------------|-------------------------------------------|
| R           | 1.5 mm     | Ring radius                               |
| w           | 100 µm     | Strut circumferential width               |
| t           | 80 µm      | Strut radial thickness                    |
| L           | 500 µm     | Strut axial length                        |
| M           | 1.0 MA/m   | Magnetisation (baseline, cold-worked SS)  |
| n_struts    | 12         | V2-2C COMSOL-matched geometry             |
| mag_mode    | radial     | Strut magnetisation direction             |
| M_COMSOL_EFF_B15 | 2.20 MA/m | Calibrated M at B₀ = 1.5 T (MRI)     |

At B₀ ≈ 1.5 T, `make_ring(B0_magnitude=1.5)` automatically substitutes
the COMSOL-calibrated magnetisation. Cell / flow defaults: 10 µm cell
radius, 4 mPa·s blood viscosity, 1.54 mm vessel radius, 0.2 m/s mean
velocity (distal MCA).

---

## Key references

- Akoun, G. & Yonnet, J.-P. (1984). 3D analytical calculation of the
  forces exerted between two cuboidal magnets. *IEEE Trans. Magn.*
  20, 1962.
- Furlani, E.P. & Ng, K.C. (2006). Analytical model of magnetic
  nanoparticle transport and capture in the microvasculature.
  *Phys. Rev. E* 73, 061919.
- Polyak, B. et al. (2008). High field gradient targeting of magnetic
  nanoparticle-loaded endothelial cells to the surfaces of steel
  stents. *PNAS* 105(2), 698–703.
- Tefft, B.J. et al. (2014). Magnetizable stent-grafts enable
  endothelial cell capture. *IEEE Trans. Magn.* 50(11), 1–4.
- Aaslid, R. et al. (1982). Noninvasive transcranial Doppler
  ultrasound recording of flow velocity. *J. Neurosurg.* 57, 769.
