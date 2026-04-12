# Resume note — Part 3 (paracrine module)

Session ran out of usage mid-Part 3. Parts 1 and 2 are committed on `main`.

## Already committed

- **`bc40c15`** — `AUDIT.md` (Part 1). Physics audit of all modules, fig21 reproduces bit-identical, 70/70 tests.
- **`81b08f4`** — `Overview.pdf` + `scripts/build_overview.py` (Part 2). 25 pages, 684 KB, built via matplotlib PdfPages, embeds figs 12/13/17/18/20/21.

## Audit findings (carry into Part 3 — do not re-derive)

1. **Polyak 2008 mislabel** — code comments claim "10 pg default" but the actual Polyak paper uses ~200 pg. Affects labels in `physics/magnetic_force.py`, `figures/fig17_*.py`, `figures/fig20_*.py`. Not fixed yet — leave for a follow-up unless the user asks.
2. **Aaslid 1982 mislabel** — `hydrodynamics.py` default `v_mean=0.2 m/s` labelled "MCA mean (Aaslid 1982)". Aaslid's actual MCA mean is 0.62 m/s. 0.2 m/s is only defensible for distal/diseased vessels.
3. **`stent_capture/physics/shear_stress.py` is orphaned & broken** — no imports anywhere, `for r in range(0, int(self.length))` where `length` is metres → empty lists. Added in commit `ec066b8`. Recommend deletion but do not delete without asking.

## Part 3 — remaining work

### Step 1 — Literature search (DO FIRST, must cite real papers)
Need these VEGF parameters with real sources:
- **Diffusion coefficient** D in tissue — target Mac Gabhann & Popel papers (~10⁻¹¹ m²/s range expected)
- **Secretion rate** per endothelial cell (molecules/cell/s or pg/cell/day)
- **Half-life / clearance rate** k — plasma vs tissue
- **Therapeutic threshold** concentration for angiogenesis (ng/mL)

Use `WebSearch` for: "VEGF diffusion coefficient tissue Mac Gabhann", "VEGF secretion rate endothelial cell", "VEGF half-life tissue clearance", "VEGF therapeutic threshold angiogenesis concentration".

### Step 2 — Implement module
Create `stent_capture/paracrine/`:
- `__init__.py`
- `transport.py` — `ParacrineField` class, 2D finite-difference solver for
  `∂C/∂t + u·∇C = D∇²C + S(r) − kC`
  on a cross-section grid. Use upwind advection, explicit Euler or Crank-Nicolson for stability, CFL check.
- `secretion.py` — `VEGFSource` maps captured-cell positions to a spatial source term S(r). Gaussian kernel around each captured cell is fine.
- `therapeutic.py` — threshold analysis, time-to-threshold at a given radial distance.
  **Important**: ChatGPT previously proposed `t_res = distance / v_slip` — this is dimensionally OK (s = m / (m/s)) but **physically wrong** here because it ignores diffusion. Document why this shortcut was rejected and solve the PDE instead.

### Step 3 — Figures
- `figures/fig22_concentration_field.py` — 2D heatmap of C(x,y) at steady state around a captured cell cluster
- `figures/fig23_concentration_vs_distance.py` — C vs radial distance, overlay therapeutic threshold
- `figures/fig24_time_to_threshold.py` — time to reach threshold at r = 100, 250, 500 μm

### Step 4 — Tests
`stent_capture/tests/test_paracrine.py`, 4 tests:
1. Zero-source → C stays zero everywhere
2. Point source, no decay, no advection → 1/r decay at steady state
3. Pure-decay limit (D = 0, u = 0) → exponential decay matches `exp(-k t)`
4. Literature agreement — steady-state C at 100–500 μm is within the therapeutic zone reported in the VEGF literature

### Step 5 — Commit
```
git -c user.name="Stanz1919" -c user.email="samlukaitis5@gmail.com" \
  commit -m "Paracrine module: VEGF signalling analysis, figs 22-24"
```

### Step 6 — Report
To user: citations found (with real page/value), therapeutic-zone radius, fig22–24 headline numbers.

## Environment reminders

- Python with full scientific stack: `/c/Users/Lenovo/anaconda3/python.exe` (NOT default `python`)
- Tests live in `stent_capture/tests/`, not `tests/`
- Windows `multiprocessing` uses spawn — any new worker fn must be top-level
- Git identity workaround: `git -c user.name="Stanz1919" -c user.email="samlukaitis5@gmail.com" commit …` (do NOT modify git config)

## How to resume

Say "continue Part 3" and I will start from Step 1 (WebSearch for VEGF parameters).
