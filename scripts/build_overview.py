"""
build_overview.py
=================
Generate Overview.pdf — a standalone walkthrough of the stent_capture
project (Stages 1-3).  Uses matplotlib's PdfPages backend so that no
additional dependencies are required beyond matplotlib itself.

Run::

    python scripts/build_overview.py

Output: Overview.pdf in the project root.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg

ROOT    = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT_PDF = ROOT / "Overview.pdf"

# ---------------------------------------------------------------------------
# Typesetting constants
# ---------------------------------------------------------------------------

PAGE_W, PAGE_H = 8.5, 11.0     # US letter, inches

TITLE_FS    = 18
H1_FS       = 15
H2_FS       = 12
BODY_FS     = 9.5
CODE_FS     = 8
CAPTION_FS  = 9
FOOT_FS     = 7.5

MARGIN_L = 0.09    # fraction of page width
MARGIN_R = 0.91
TOP_Y    = 0.93
BOT_Y    = 0.06

BODY_WIDTH = 95    # textwrap char width for serif body font
CODE_WIDTH = 88

SERIF     = {"family": "serif"}
MONO      = {"family": "monospace"}
SERIF_BF  = {"family": "serif", "weight": "bold"}

_page_no = [0]


def _new_page(pdf: PdfPages):
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    _page_no[0] += 1
    fig.text(0.5, 0.025, f"— {_page_no[0]} —",
             ha="center", fontsize=FOOT_FS, color="gray", **SERIF)
    return fig


def _finish(pdf: PdfPages, fig):
    pdf.savefig(fig)
    plt.close(fig)


def _wrap(txt: str, width: int = BODY_WIDTH) -> str:
    """Wrap preserving paragraph breaks (blank lines)."""
    paras = txt.split("\n\n")
    out = []
    for p in paras:
        # collapse internal newlines then fill
        flat = " ".join(line.strip() for line in p.splitlines() if line.strip())
        if flat:
            out.append(textwrap.fill(flat, width=width))
    return "\n\n".join(out)


def _body(fig, y: float, text: str, fs: float = BODY_FS,
          mono: bool = False, bold: bool = False) -> float:
    """Place wrapped body text starting at y and return the y after it."""
    kwargs = MONO if mono else (SERIF_BF if bold else SERIF)
    wrapped = _wrap(text, CODE_WIDTH if mono else BODY_WIDTH)
    n_lines = wrapped.count("\n") + 1
    fig.text(MARGIN_L, y, wrapped, ha="left", va="top",
             fontsize=fs, **kwargs)
    return y - 0.015 * n_lines - 0.01


def _heading(fig, y: float, text: str, level: int = 1) -> float:
    fs = H1_FS if level == 1 else H2_FS
    fig.text(MARGIN_L, y, text, ha="left", va="top",
             fontsize=fs, **SERIF_BF)
    return y - (0.038 if level == 1 else 0.028)


def _eqn(fig, y: float, latex: str, fs: float = 11) -> float:
    """Display-style equation, centred."""
    fig.text(0.5, y, latex, ha="center", va="top", fontsize=fs, **SERIF)
    return y - 0.040


# ---------------------------------------------------------------------------
# Page builders
# ---------------------------------------------------------------------------

def page_title(pdf: PdfPages):
    fig = _new_page(pdf)
    fig.text(0.5, 0.72, "Stent Capture Model",
             ha="center", fontsize=32, **SERIF_BF)
    fig.text(0.5, 0.66, "Magnetic capture of SPION-labelled cells\n"
             "on a cerebral ferromagnetic stent",
             ha="center", fontsize=15, **SERIF)
    fig.text(0.5, 0.56, "Project Overview",
             ha="center", fontsize=18, **SERIF)
    fig.text(0.5, 0.30,
             "A multi-stage Python framework for high-resolution\n"
             "trajectory-level exploration of magnetic cell targeting,\n"
             "built to complement the COMSOL finite-element model.",
             ha="center", fontsize=11, color="#333", **SERIF)
    fig.text(0.5, 0.14, "v0.5.0 — 2026-04-11",
             ha="center", fontsize=10, color="#888", **SERIF)
    _finish(pdf, fig)


def page_toc(pdf: PdfPages):
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "Contents", level=1)
    items = [
        ("1  Introduction",                               "  3"),
        ("2  Physical model",                              "  5"),
        ("    2.1  Geometry and coordinate frames",        "  5"),
        ("    2.2  The Akoun-Yonnet analytical field",     "  6"),
        ("    2.3  Magnetic force on SPION cells",         "  7"),
        ("    2.4  Hydrodynamics: Poiseuille and Stokes",  "  8"),
        ("    2.5  Trajectory ODE derivation",             "  9"),
        ("3  Implementation",                              " 10"),
        ("    3.1  Module structure",                      " 10"),
        ("    3.2  Module-by-module summary",              " 11"),
        ("    3.3  Default parameters",                    " 13"),
        ("    3.4  Test suite",                            " 14"),
        ("4  Results",                                     " 15"),
        ("    4.1  Fig 12 – external field comparison",    " 15"),
        ("    4.2  Fig 13 – force parameter",              " 16"),
        ("    4.3  Fig 17 – SPION loading sweep",          " 17"),
        ("    4.4  Fig 18 – single-cell trajectory",       " 18"),
        ("    4.5  Fig 20 – capture efficiency curves",    " 19"),
        ("    4.6  Fig 21 – static vs trajectory",         " 20"),
        ("5  Limitations and assumptions",                 " 22"),
        ("6  How to use the code",                         " 24"),
        ("7  References",                                  " 26"),
    ]
    for line, pg in items:
        fig.text(MARGIN_L, y, line, fontsize=BODY_FS, **SERIF)
        fig.text(MARGIN_R, y, pg, fontsize=BODY_FS, ha="right", **SERIF)
        y -= 0.027
    _finish(pdf, fig)


# --- Section 1 --------------------------------------------------------------

def section_1(pdf: PdfPages):
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "1  Introduction", level=1)
    body = (
        "This document is a walkthrough and overview of the 'stent_capture_project': "
        "a Python simulation framework that models the magnetic capture of "
        "SPION-labelled endothelial cells on a ferromagnetic stent deployed "
        "under cerebral flow conditions. It is intended as a standalone"
        "analaytical tool to aid the dissertation. This overview is designed to provide"
        "understanding of the physics, code, and the numerical "
        "results without having to read the source tree.\n\n"
        "The code simulates a ring of eight uniformly magnetised "
        "rectangular struts that sit inside a cerebral artery "
        "modelled with Poiseuille flow. A uniform external field B0 is "
        "superposed on the stent's magnetostatic field. SPION-"
        "labelled cells are introduced upstream and their trajectories are "
        "integrated under a terminal-velocity ODE that balances Stokes drag "
        "against the magnetic dipole force. The framework evaluates both "
        "the static capture criterion |F_mag| > |F_drag| and the full "
        "trajectory integration, and compares them at MCA-representative "
        "flow and loading conditions.\n\n"
        "The analysis of magnetic capture uses a "
        "finite-element COMSOL model to solve the full magnetic field gradient"
        "around the stent. However, the model is expensive: each configuration "
        "costs minutes to hours, making parametric sweeps infeasible. "
        "The 'stent_capture_project' framework is a tool built on fast "
        "analytical expressions (Akoun & Yonnet 1984) so that parameter "
        "sweeps across loadings, velocities, and stent geometries run in "
        "seconds to minutes rather than hours; trajectory mapping of "
        "cells can be integrated without meshing; and the "
        "results serve as an independent numerical cross-check on the "
        "COMSOL simulations."
    )
    y = _body(fig, y, body)
    _finish(pdf, fig)

    # Page 2 of Section 1
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "1  Introduction (cont.)", level=1)
    body2 = (
        "The project is  organised into 3"
        "stages, each building on the previous "
        "one, leaving behind its own tests and figures so that "
        "regressions surface immediately.\n\n"
        "Stage 1 — Field model. Implements the Akoun & Yonnet (1984) closed-"
        "form expressions for the B-field of a uniformly magnetised "
        "rectangular prism, rotates and superposes 8 struts into a ring, "
        "and layers an optional uniform external field B0 on top. Stage 1 "
        "establishes the field primitives (field_at, B_magnitude, grad_B) "
        "that every later stage consumes through a shared TotalField "
        "interface.\n\n"
        "Stage 2 — Static capture criterion. Adds the SPION cell model, "
        "the Furlani & Ng (2006) scalar-gradient magnetic force, the "
        "Poiseuille-flow Stokes drag, and the scalar capture criterion "
        "|F_mag| > |F_drag|. Produces 2D capture maps and SPION loading "
        "sweeps. Shows that the static criterion predicts no capture in "
        "the lumen at MCA flow for 10 pg loading — motivating Stage 3.\n\n"
        "Stage 3 — Trajectory integration. Adds a low-Re terminal-velocity "
        "ODE solver built on scipy.integrate.solve_ivp, event-based capture "
        "detection, a multiprocessing capture-efficiency sweep, and the "
        "headline static-vs-trajectory comparison figure. Demonstrates that "
        "trajectory drift extends the effective capture range by ~14x over "
        "the static criterion at MCA conditions.\n\n"
        "Stage 4 — Geometry optimisation (planned). Would sweep "
        "strut count, thickness, length, and magnetisation to identify "
        "optimal stent designs. Not in scope for this overview.\n\n"
        "The rest of this document is organised to mirror that structure: "
        "Section 2 derives the physics; Section 3 describes the code that "
        "implements it; Section 4 walks through the key figures produced; "
        "Section 5 lists the known limitations; Section 6 explains how to "
        "run everything; and Section 7 collects the literature references."
    )
    _body(fig, y, body2)
    _finish(pdf, fig)


# --- Section 2 --------------------------------------------------------------

def section_2(pdf: PdfPages):
    # Page: 2.1 Geometry
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "2  Physical model", level=1)
    y = _heading(fig, y, "2.1  Geometry and coordinate frames", level=2)
    body = (
        "The stent is modelled as a ring of n = 8 identical rectangular "
        "struts equally spaced around a circle of radius R = 1.5 mm in the "
        "z = 0 plane. Each strut has radial thickness t = 80 um, "
        "circumferential width w = 100 um, and axial length L = 500 um. "
        "The ring is placed inside a straight cylindrical vessel of radius "
        "R_ves = 1.54 mm, so the stent outer surface is flush with the "
        "vessel wall and the inner lumen has radius r_inner = R - t/2 "
        "~ 1.46 mm.\n\n"
        "Two coordinate frames are used. The global frame (x, y, z) has "
        "the stent axis along +z and the ring centred at the origin. "
        "Each strut i has a local frame (x', y', z') aligned with its "
        "own radial, circumferential, and axial directions — this is the "
        "frame in which the Akoun-Yonnet kernel is evaluated. The global "
        "field is recovered by rotating the local field by the strut's "
        "azimuth theta_i = 2*pi*i / n and summing the 8 contributions.\n\n"
        "Magnetisation direction is a user parameter (mag_mode = radial / "
        "circumferential / axial). All results in this document use the "
        "default mag_mode = 'radial', i.e. each strut is magnetised "
        "outward in its own local x direction. This matches the geometry "
        "of a cold-worked ferromagnetic stent wire drawn radially across "
        "the ring."
    )
    y = _body(fig, y, body)
    y = _heading(fig, y, "Boundary conditions", level=2)
    bcs = (
        "(i) Lumen: no-slip at r = R_ves, parabolic Poiseuille profile "
        "inside. (ii) Field kernel: free-space exterior; the Akoun-Yonnet "
        "expressions are exact for points outside the prism. (iii) Cell "
        "trajectories: start at an injection point upstream of the ring "
        "(typically z = -2 mm) and terminate on one of three events — "
        "reaching z_end = +2 mm (escape), coming within the capture "
        "tolerance of any strut (proximity capture), or reaching the "
        "lumen inner radius (wall capture). The stent axial extent is "
        "L = 500 um, so the 4 mm integration window brackets 8x the stent "
        "length on each side."
    )
    _body(fig, y, bcs)
    _finish(pdf, fig)

    # Page: 2.2 Akoun-Yonnet
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "2.2  The Akoun-Yonnet analytical field", level=1)
    body = (
        "For a uniformly magnetised rectangular prism with half-dimensions "
        "(a, b, c) and magnetisation M, Akoun & Yonnet (1984) derived a "
        "closed-form expression for the B-field at any exterior point by "
        "integrating the scalar potential of the magnetic surface charges "
        "sigma_m = M . n_hat over each of the six rectangular faces. The "
        "result is a sum over the 8 corners of the prism with alternating "
        "signs:"
    )
    y = _body(fig, y, body)
    y = _eqn(fig, y,
             r"$B_i(\mathbf{r}) = \dfrac{\mu_0}{4\pi}\;"
             r"\sum_{k,l,m=0}^{1} (-1)^{k+l+m}\,"
             r"\mathcal{K}_i\!\left(U_{kl m},V_{kl m},W_{kl m}\right)$")
    body2 = (
        "where U, V, W are the signed corner-to-observation coordinates "
        "and the kernel K_i is a combination of arctangent and logarithmic "
        "terms. For a prism magnetised along local +z, the x-component is "
        "for example"
    )
    y -= 0.005
    y = _body(fig, y, body2)
    y = _eqn(fig, y,
             r"$\mathcal{K}_x(U,V,W) = -\ln\!\left(V + R\right),\quad "
             r"R = \sqrt{U^2 + V^2 + W^2}$")
    body3 = (
        "with analogous cyclic permutations for the Mx and My contributions "
        "and for By, Bz. The full kernel is implemented in "
        "core/field_model.py::_akoun_yonnet_local, with a global sign flip "
        "applied at the return statement so that B = mu_0*H follows the "
        "H = -grad(phi_m) convention (this was a latent bug in the original "
        "legacy script, fixed in Stage 1 and now guarded by the "
        "TestSuperposition regression test to relative tolerance 1e-10).\n\n"
        "Why analytical, not FEM? The alternative would be to mesh the "
        "stent in COMSOL, solve the magnetostatic problem via finite "
        "elements, and interpolate the result onto the trajectory "
        "integration grid. Analytical expressions avoid all of this: there "
        "is no mesh to converge, no linear system to solve, and evaluation "
        "cost is O(n_struts * n_points) with a small constant. Each field "
        "evaluation is numerically exact (to float64 precision) regardless "
        "of the distance to the strut — there is no 'near-strut refinement' "
        "problem. For a framework whose primary purpose is dense parameter "
        "sweeps and trajectory ensembles, the analytical approach is "
        "strictly better."
    )
    _body(fig, y, body3)
    _finish(pdf, fig)

    # Page: 2.3 Magnetic force
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "2.3  Magnetic force on SPION-labelled cells", level=1)
    body = (
        "The force on a small superparamagnetic particle of volume V_p and "
        "volume susceptibility chi_p in a static external field follows "
        "from the gradient of the dipole energy (Furlani & Ng 2006, Eq. 2):"
    )
    y = _body(fig, y, body)
    y = _eqn(fig, y,
             r"$\mathbf{F} = \mu_0\, V_p\, (\mathbf{M}_p \cdot \nabla)\,"
             r"\mathbf{H}_{\mathrm{tot}}$")
    body2 = (
        "In the linear-susceptibility regime (M_p = chi_p H), and for a "
        "free-space curl-free field (grad x H = 0), this simplifies via "
        "the vector identity (H . grad) H = grad(|H|^2 / 2) to"
    )
    y = _body(fig, y, body2)
    y = _eqn(fig, y,
             r"$\mathbf{F} = \dfrac{V_p\,\chi_p}{2\mu_0}\,"
             r"\nabla |\mathbf{B}_{\mathrm{tot}}|^{2} = "
             r"\dfrac{V_p\,\chi_p}{\mu_0}\,"
             r"|\mathbf{B}_{\mathrm{tot}}|\,\nabla|\mathbf{B}_{\mathrm{tot}}|$")
    body3 = (
        "This is the scalar-gradient approximation used throughout the "
        "magnetic drug and cell targeting literature (Furlani & Ng 2006; "
        "Polyak et al. 2008; Tefft et al. 2014). The force direction is "
        "along grad|B_tot|, which points from low |B| toward high |B| — "
        "i.e. from the lumen interior toward the strut surface — and the "
        "magnitude scales with the product |B|*|grad|B||, called the "
        "force parameter.\n\n"
        "Two subtleties matter when an external uniform field B0 is added. "
        "First, V_p is the SPION volume, not the cell volume; the cell "
        "cytoplasm contributes no magnetic force. Second, grad|B_tot| is "
        "NOT equal to grad|B_stent| even though grad B0 = 0, because the "
        "gradient of a magnitude is nonlinear:\n"
    )
    y = _body(fig, y, body3)
    y = _eqn(fig, y,
             r"$\dfrac{\partial\,|B_{\mathrm{tot}}|}{\partial x_i} = "
             r"\dfrac{\mathbf{B}_{\mathrm{tot}}\cdot "
             r"(\partial\mathbf{B}_{\mathrm{stent}}/\partial x_i)}"
             r"{|\mathbf{B}_{\mathrm{tot}}|}$")
    body4 = (
        "When B0 is axial and the stent is magnetised radially, B_tot "
        "rotates toward +z and the projection of d B_stent / d x_i onto "
        "B_tot becomes small — so |grad|B_tot|| is SUPPRESSED relative to "
        "|grad|B_stent||. Counter-intuitively this is physically correct: "
        "adding a uniform field does reduce the scalar gradient, but it "
        "ALSO increases |B_tot|, and the product |B|*|grad|B|| (which is "
        "what enters the force) actually grows. This is why the code "
        "exposes the force parameter |B|*grad|B| as the relevant capture "
        "metric rather than grad|B| alone."
    )
    _body(fig, y, body4)
    _finish(pdf, fig)

    # Page: 2.4 Hydrodynamics
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "2.4  Hydrodynamics: Poiseuille and Stokes", level=1)
    body = (
        "The blood flow through the stented vessel is modelled as steady "
        "axisymmetric Poiseuille flow. For a Newtonian fluid of viscosity "
        "eta in a cylinder of radius R_ves, the velocity profile is "
        "parabolic:"
    )
    y = _body(fig, y, body)
    y = _eqn(fig, y,
             r"$v_z(r) = 2\,\bar{v}\left(1 - (r/R_{\mathrm{ves}})^{2}\right)$")
    body2 = (
        "where v_bar is the cross-sectional mean velocity. The peak "
        "centreline velocity is 2*v_bar and the wall shear stress is "
        "tau_w = 4 eta v_bar / R_ves.\n\n"
        "Why steady Poiseuille, not Womersley? Cerebral arteries pulsate "
        "with a Womersley number alpha = R sqrt(omega rho / eta) ~ 2 at "
        "typical heart rates, so the flow is in the quasi-steady regime "
        "where the instantaneous profile is close to the steady Poiseuille "
        "profile at each instant. Time-averaging the trajectory over a "
        "cardiac cycle would smooth the force balance but not change the "
        "population-level capture efficiency dramatically. Steady "
        "Poiseuille is therefore a reasonable first approximation for a "
        "time-averaged capture analysis.\n\n"
        "The cell experiences a viscous drag from the surrounding blood. "
        "For a rigid sphere of radius R_cell moving with velocity v_cell "
        "through fluid with local velocity v_blood, the drag force is "
        "given by Stokes' law:"
    )
    y = _body(fig, y, body2)
    y = _eqn(fig, y,
             r"$\mathbf{F}_{\mathrm{drag}} = "
             r"6\pi\eta R_{\mathrm{cell}} "
             r"\left(\mathbf{v}_{\mathrm{blood}} - \mathbf{v}_{\mathrm{cell}}\right)$")
    body3 = (
        "The Stokes formula is valid in the low particle Reynolds number "
        "limit Re_p = rho v (2 R_cell) / eta << 1. For a 10 um cell in "
        "whole blood (rho ~ 1060, eta ~ 4 mPa s) at v = 0.2 m/s, "
        "Re_p ~ 1, which is at the upper edge of Stokes validity but still "
        "close enough that the error is a few percent at worst. For the "
        "same cell with a radial magnetic-drift velocity component, which "
        "is typically 1-3 orders of magnitude slower than the axial blood "
        "velocity, Re_p is well below unity and Stokes is comfortably "
        "valid."
    )
    _body(fig, y, body3)
    _finish(pdf, fig)

    # Page: 2.5 Trajectory ODE
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "2.5  Trajectory ODE derivation", level=1)
    body = (
        "Newton's second law for a cell of mass m_cell in the combined "
        "magnetic and hydrodynamic force field reads"
    )
    y = _body(fig, y, body)
    y = _eqn(fig, y,
             r"$m_{\mathrm{cell}}\,\dfrac{d\mathbf{v}_{\mathrm{cell}}}{dt} "
             r"= \mathbf{F}_{\mathrm{mag}} + "
             r"\mathbf{F}_{\mathrm{drag}} + "
             r"\mathbf{F}_{\mathrm{buoyancy}} + \ldots$")
    body2 = (
        "In the low-Reynolds-number limit the inertial term m dv/dt is "
        "negligible compared to the viscous drag, and the cell "
        "instantaneously adopts the terminal velocity at which the drag "
        "exactly balances all other forces. Setting the left-hand side to "
        "zero and solving for v_cell:"
    )
    y = _body(fig, y, body2)
    y = _eqn(fig, y,
             r"$\mathbf{v}_{\mathrm{cell}} = "
             r"\mathbf{v}_{\mathrm{blood}}(\mathbf{r}) + "
             r"\dfrac{\mathbf{F}_{\mathrm{mag}}(\mathbf{r})}"
             r"{6\pi\eta R_{\mathrm{cell}}}$")
    body3 = (
        "Since the position r is updated by integrating v_cell in time, "
        "this gives a first-order ODE for the cell trajectory:"
    )
    y = _body(fig, y, body3)
    y = _eqn(fig, y,
             r"$\dfrac{d\mathbf{r}}{dt} = "
             r"\mathbf{v}_{\mathrm{blood}}(\mathbf{r}) + "
             r"\dfrac{\mathbf{F}_{\mathrm{mag}}(\mathbf{r})}"
             r"{6\pi\eta R_{\mathrm{cell}}}$")
    body4 = (
        "This is the equation integrated by "
        "simulation/trajectories.py::integrate_trajectory using "
        "scipy.integrate.solve_ivp with method 'RK45' (embedded 4/5-order "
        "Runge-Kutta), relative tolerance 1e-6, absolute tolerance 1e-9, "
        "and event-based termination.\n\n"
        "Why the inertial term is negligible. The viscous relaxation time "
        "of a 10 um cell in blood is tau_v = m_cell / (6 pi eta R_cell) "
        "~ 1 us. The magnetic force changes over a length scale ~ 100 um, "
        "so the time it takes the cell to traverse that scale at drift "
        "velocities of ~1 mm/s is ~100 ms — five orders of magnitude "
        "longer than tau_v. The cell is therefore always in terminal-"
        "velocity equilibrium, and neglecting m dv/dt introduces an error "
        "of order (tau_v / 100 ms) ~ 1e-5, far below every other "
        "approximation in the model.\n\n"
        "Three terminal events are defined in the solver. The escape event "
        "fires when z reaches z_end = +2 mm downstream — these cells are "
        "classified as 'escaped'. The wall event fires when the cell's "
        "transverse radius sqrt(x^2 + y^2) reaches R - t/2, the stent "
        "inner surface — these are classified as 'captured'. The strut-"
        "proximity event catches cells that approach a strut from outside "
        "the ring (a cylindrical approximation of the rectangular strut is "
        "used). A fourth 'timeout' outcome is raised if the integrator "
        "hits max_time = 1 s without any event firing."
    )
    _body(fig, y, body4)
    _finish(pdf, fig)


# --- Section 3 --------------------------------------------------------------

MODULE_TREE = """\
stent_capture_project/
├── scripts/
│   └── build_overview.py       this document
├── stent_capture/
│   ├── core/
│   │   ├── field_model.py      StentRing + Akoun-Yonnet kernel
│   │   └── gradient.py         compute_gradient_{magnitude,vector}
│   ├── physics/
│   │   ├── external_field.py   UniformExternalField, TotalField
│   │   ├── magnetic_force.py   SPIONLabelledCell, magnetic_force()
│   │   ├── hydrodynamics.py    BloodFlow (Poiseuille), stokes_drag()
│   │   ├── capture_criterion.py capture_map(), capture_distance()
│   │   └── shear_stress.py     [orphan; see AUDIT.md]
│   ├── simulation/
│   │   ├── trajectories.py     integrate_trajectory, CellTrajectory
│   │   └── capture_efficiency.py sweep_injection_line + sweeps
│   ├── figures/
│   │   ├── common.py           DEFAULTS, make_ring(), save_fig()
│   │   ├── style.py            rcParams (dissertation style)
│   │   ├── fig01-11_*.py       Stage 1 field & geometry
│   │   ├── fig12-16_*.py       Stage 2 force, drag, capture maps
│   │   ├── fig17_*.py          Stage 2.5 loading sweep
│   │   ├── fig18-19_*.py       Stage 3a trajectories
│   │   ├── fig20_*.py          Stage 3b efficiency curves
│   │   └── fig21_*.py          Stage 3c headline comparison
│   └── tests/                  70 pytest tests, 6 files
├── results/                    PNG + PDF output
├── AUDIT.md                    final physics review (this session)
├── CHANGELOG.md
├── README.md
└── Overview.pdf                this document
"""


def section_3(pdf: PdfPages):
    # Page: 3.1 structure
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "3  Implementation", level=1)
    y = _heading(fig, y, "3.1  Module structure", level=2)
    intro = (
        "The codebase is a pure-Python package with a strict "
        "core -> physics -> simulation -> figures dependency direction. "
        "Every module has a single responsibility, every class accepts "
        "dependencies by injection, and every physics quantity is carried "
        "in SI units throughout. No module in core/ imports from physics/, "
        "no module in physics/ imports from simulation/, and figures/ "
        "imports from all lower levels but from nothing else. The tree "
        "below shows the directory layout as of v0.5.0:"
    )
    y = _body(fig, y, intro)
    fig.text(MARGIN_L, y, MODULE_TREE, ha="left", va="top",
             fontsize=CODE_FS, **MONO)
    _finish(pdf, fig)

    # Page: 3.2 module summaries (multi-page)
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "3.2  Module-by-module summary", level=1)
    summaries = [
        ("core/field_model.py",
         "StentRing(n, R, w, t, L, M, mag_mode, assume_saturation). "
         "Public method field_at(points) -> (N, 3) B-vectors. Internal "
         "kernel _akoun_yonnet_local implements the 8-corner closed-form "
         "evaluation. Rotations between each strut's local frame and the "
         "global frame are cached at construction. assume_saturation = True "
         "signals that the supplied M is the saturation value."),
        ("core/gradient.py",
         "compute_gradient_vector(field_func, points, dx) and "
         "compute_gradient_magnitude(...) — 6-evaluation central finite "
         "differences with default dx = 500 nm. Field-model agnostic: "
         "accepts any callable mapping (N, 3) points to (N, 3) B-vectors."),
        ("physics/external_field.py",
         "UniformExternalField(B0_vector) for spatially uniform B0. "
         "TotalField(stent_ring, external_field) for B_stent + B0. Both "
         "expose the same (field_at, B_field, B_magnitude, grad_B) API "
         "as StentRing, so they are drop-in compatible with downstream "
         "force and trajectory code."),
        ("physics/magnetic_force.py",
         "SPIONLabelledCell(radius, spion_mass, spion_susceptibility, "
         "spion_density). magnetic_force(cell, total_field, points) "
         "returns the (N, 3) force vector using the curl-free scalar-"
         "gradient formula F = (V_p chi / mu0) |B| grad|B|. Uses the "
         "SPION volume, not the cell volume."),
        ("physics/hydrodynamics.py",
         "BloodFlow(vessel_radius, mean_velocity, viscosity, density) "
         "with velocity_at, shear_rate_at, wall_shear_stress. "
         "stokes_drag(cell, blood_flow, points, cell_velocities) returns "
         "the (N, 3) drag force. For stationary cells the default "
         "cell_velocities = 0 gives the worst-case drag."),
        ("physics/capture_criterion.py",
         "capture_map(cell, tf, flow, points) returns a dict of "
         "F_mag_vec / F_drag_vec / scalar magnitudes / margin / captured. "
         "capture_distance(cell, tf, flow, direction='inward') sweeps "
         "from the stent inner surface toward the vessel centre and "
         "returns the outermost capture radius."),
    ]
    for name, desc in summaries:
        y = _body(fig, y, name, bold=True)
        y = _body(fig, y, desc)
        y -= 0.005
    _finish(pdf, fig)

    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "3.2  Module summaries (cont.)", level=1)
    more = [
        ("simulation/trajectories.py",
         "CellTrajectory result object (positions, times, status, cell, "
         "field, flow references). integrate_trajectory uses "
         "scipy.integrate.solve_ivp(method='RK45', rtol=1e-6, atol=1e-9) "
         "with three event functions (escape, strut proximity, wall "
         "contact). Returns status in {'captured', 'escaped', 'error'}. "
         "The RHS closure recomputes the 6-evaluation FD magnetic force "
         "at every solver call — no pre-computed field grid."),
        ("simulation/capture_efficiency.py",
         "sweep_injection_line dispatches N trajectories to a "
         "multiprocessing.Pool (Windows 'spawn'-compatible worker at "
         "module level). Two sweep helpers: "
         "capture_efficiency_vs_velocity and "
         "capture_efficiency_vs_loading, each accepting a cell_factory "
         "or fresh BloodFlow per sweep step."),
        ("figures/common.py",
         "Central DEFAULTS dict with R, w, t, L, M, n_struts. "
         "make_ring(**overrides) factory, threshold_lines(ax) helper, "
         "save_fig(fig, stem) that writes both PNG and PDF to results/."),
        ("figures/style.py",
         "Side-effect import that configures matplotlib rcParams for "
         "serif fonts, 300 dpi, grid on, tight layout, consistent font "
         "sizes. Every figure module imports common (which in turn "
         "imports style) before plotting."),
        ("physics/shear_stress.py",
         "Orphan file. Not imported by any module or test. Contains "
         "broken loops (int(length_in_metres) == 0). Recommended for "
         "deletion in the AUDIT report. BloodFlow.wall_shear_stress "
         "already provides the correct tau_w."),
    ]
    for name, desc in more:
        y = _body(fig, y, name, bold=True)
        y = _body(fig, y, desc)
        y -= 0.005
    _finish(pdf, fig)

    # Page: 3.3 defaults table
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "3.3  Default parameters", level=1)
    intro = (
        "The table below lists every numerical parameter carried by the "
        "production code with the value it takes and the literature "
        "source. See AUDIT.md for the full per-parameter verification "
        "and caveats."
    )
    y = _body(fig, y, intro)
    table = (
        "  Parameter            Value        Source / notes\n"
        "  ------------------   ----------   -------------------------------------\n"
        "  n_struts             8            dissertation stent geometry\n"
        "  R (ring radius)      1.5 mm       cerebral stent\n"
        "  w (strut width)      100 um       fabrication-feasible\n"
        "  t (strut thickness)  80 um        cerebral stent\n"
        "  L (strut length)     500 um       cerebral stent\n"
        "  M (magnetisation)    1 MA/m       cold-worked 304 / 430 SS\n"
        "  mag_mode             radial       each strut magnetised outward\n"
        "\n"
        "  Cell radius          10 um        standard endothelial cell\n"
        "  SPION mass / cell    10 pg        Chorny uptake curve lower end [*]\n"
        "  chi_spion            2.0          Furlani & Ng 2006 Table 1\n"
        "  SPION density        5170 kg/m^3  magnetite Fe3O4 handbook\n"
        "\n"
        "  Vessel radius        1.54 mm      M1 MCA upper range\n"
        "  Mean velocity        0.2 m/s      see note [**]\n"
        "  Viscosity            4 mPa s      whole blood, ~40% Hct\n"
        "  Density              1060 kg/m^3  whole blood\n"
        "\n"
        "  B0 (external)        0.5 T        dissertation sweep centre\n"
        "  mag_mode B0          axial (+z)   external coil geometry\n"
        "\n"
        "  RK45 rtol            1e-6         solve_ivp default tightened\n"
        "  RK45 atol            1e-9         to ~10 nm position accuracy\n"
        "  FD dx (gradient)     500 nm       validated dx sweep 1e-8..5e-6\n"
        "  Capture tolerance    5 um         strut-proximity event margin\n"
    )
    fig.text(MARGIN_L, y, table, ha="left", va="top",
             fontsize=CODE_FS, **MONO)
    y -= 0.60
    notes = (
        "[*] 'Polyak 2008 default = 10 pg' is mislabelled throughout "
        "figs 17 / 20 / 21. Polyak et al. (2008) actually selected "
        "0.2 ng = 200 pg per cell. See AUDIT.md issue #2.\n\n"
        "[**] The 'MCA mean = 0.2 m/s' attribution to Aaslid et al. 1982 "
        "is wrong — Aaslid's 1982 TCD value is 0.62 +/- 0.12 m/s. The "
        "0.2 m/s figure is defensible only for distal M2/M3 or diseased "
        "vessels; see AUDIT.md issue #3."
    )
    _body(fig, y, notes)
    _finish(pdf, fig)

    # Page: 3.4 tests
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "3.4  Test suite", level=1)
    intro = (
        "The test suite has 70 pytest tests distributed across six files "
        "and runs in ~2 minutes on a laptop. Every test is a regression "
        "guard against a specific physics claim made in the module "
        "docstrings. The tests are the primary evidence that the code "
        "does what the physics says it does."
    )
    y = _body(fig, y, intro)
    tests = [
        ("test_external_field.py  (19 tests)",
         "TotalField(ring, None) reproduces StentRing.grad_B at 20 points "
         "to rtol 1e-6. Single-strut superposition |B_total| = |B_stent| "
         "+ B0 when B0 parallel M at rtol 1e-10. 45 degree rotation "
         "symmetry for 8-strut ring at rtol 1e-3. Far-field limit "
         "|B_total| -> |B0| and |grad| -> 0 at r = 50 mm. "
         "UniformExternalField broadcast, magnitude, input validation. "
         "TestFDStability sweeps dx = 1e-8..5e-6 and verifies ratio ~ 1."),
        ("test_magnetic_force.py  (15 tests)",
         "Zero-gradient gives zero force. Force points toward the strut. "
         "Linear scaling in chi and V_spion. Cell property correctness. "
         "Order-of-magnitude check at 100 um: |F| ~ 340 pN for 10 pg / "
         "0.5 T."),
        ("test_hydrodynamics.py  (12 tests)",
         "Poiseuille profile at centreline, wall, and mean point. Shear "
         "rate magnitude. Stokes drag scaling in R_cell and velocity. "
         "Wall shear stress formula. Vector direction along flow."),
        ("test_capture_criterion.py  (10 tests)",
         "Dict structure / keys. Capture at the stent surface. No capture "
         "at the vessel centreline. Range decreases with velocity. Range "
         "increases with B0. Inward-sweep loading monotonicity. Maximum "
         "|F_mag|/|F_drag| < 1.0 in the lumen at v = 0.05 m/s "
         "(regression guard for the static no-capture finding)."),
        ("test_trajectories.py  (5 tests)",
         "Straight-line drift with zero magnetisation (|x, y| < 100 nm; "
         "z_end within 1 um). Poiseuille velocity matches prediction at "
         "10 intermediate points. 200 pg near-wall cell captured within "
         "20 um. 10 pg centreline cell escapes. Timeout safety."),
        ("test_capture_efficiency.py  (5 tests)",
         "M = 0 -> efficiency 0. Summary dict key integrity. "
         "Non-zero capture at 200 pg / v = 0.05 m/s near-wall. "
         "Monotonicity in velocity (slow > fast). "
         "Monotonicity in loading (high > low)."),
    ]
    for name, desc in tests:
        y = _body(fig, y, name, bold=True)
        y = _body(fig, y, desc)
        y -= 0.004
    _finish(pdf, fig)


# --- Section 4: figures -----------------------------------------------------

FIGURE_PAGES = [
    ("fig12_external_field_comparison",
     "4.1  Fig 12 — External field comparison",
     "2x2 panel showing the effect of adding a uniform external field B0 "
     "on the stent's gradient and force parameter profiles. Panel (a) is "
     "the B0 = 0 reference. Panel (b) shows B0 = 0.5 T applied axially "
     "(perpendicular to the stent's radial magnetisation) — the scalar "
     "|grad|B|| is suppressed by ~5x, but |B| grows ~17x, and the product "
     "|B|*|grad|B|| (the force parameter) increases net. Panel (c) shows "
     "B0 transverse (parallel to M) where |B| and |grad|B|| simply add. "
     "Panel (d) sweeps capture distance vs B0 magnitude from 0 to 1 T. "
     "Physical interpretation: the correct capture metric is the force "
     "parameter |B|*grad|B|, not grad|B| alone, and an axial external "
     "coil is a valid enhancement strategy despite counter-intuitively "
     "reducing the gradient magnitude."),
    ("fig13_force_parameter",
     "4.2  Fig 13 — Force parameter profiles",
     "Force parameter |B|*|grad|B|| vs distance from the stent outer "
     "surface along the through-strut radial axis, for five values of "
     "the external field B0 from 0 to 1 T axial. Y-axis log scale. The "
     "B0 = 0 curve falls from ~1e5 T^2/m at 10 um to ~1e2 T^2/m at "
     "1 mm — a 1000x drop over a millimetre. The B0 = 0.5 T curve sits "
     "~3x above the B0 = 0 curve at 200 um and ~55x above at 500 um, "
     "quantifying the long-range enhancement from the external field. "
     "This is the figure that justifies running the rest of the project "
     "at B0 = 0.5 T rather than B0 = 0."),
    ("fig17_spion_loading_sweep",
     "4.3  Fig 17 — SPION loading sweep",
     "Stage 2.5. Panel (a): capture distance (um from stent inner "
     "surface) vs SPION loading per cell (1-300 pg log-spaced), for "
     "three velocities (0.05, 0.2, 0.5 m/s) at B0 = 0.5 T axial. Panel "
     "(b): force ratio |F_mag|/|F_drag| at 5 um from the stent inner "
     "surface vs loading, log-log scale, with the ratio = 1 capture "
     "threshold marked. Key result: the STATIC criterion predicts NO "
     "capture anywhere in the lumen at 10 pg / MCA flow; static capture "
     "distance only exceeds ~20 um at loadings above ~50 pg and slow "
     "flow. This is the motivation for the trajectory integration in "
     "Stage 3 — the static picture is pessimistic."),
    ("fig18_single_trajectory",
     "4.4  Fig 18 — Single-cell trajectory (200 pg captured vs 10 pg "
     "escaping)",
     "Stage 3a. Two trajectories from identical injection at "
     "(1.3 mm, 0, -2 mm), v_mean = 0.05 m/s, B0 = 0.5 T, differing only "
     "in SPION loading. Panel (a) is a 3D view with the 8 strut prisms "
     "rendered. The 200 pg cell (blue, captured at t = 62.3 ms) curves "
     "radially inward as it approaches the stent and contacts the wall; "
     "the 10 pg cell (red, escaped) drifts axially downstream with "
     "negligible radial displacement. Panel (b) is the x-z side view. "
     "Panel (c) is the time-series of the 200 pg trajectory: radial "
     "position r(t), cell speed, and |F_mag|(t) on log scale. Key "
     "observation: at injection (z = -2 mm) the static criterion "
     "predicts no capture, yet the 200 pg cell accumulates enough "
     "radial drift during its transit to be captured at the stent. This "
     "is the fundamental mechanism the static criterion misses."),
    ("fig20_capture_efficiency",
     "4.5  Fig 20 — Capture efficiency curves",
     "Stage 3b. Two-panel population study, 20 cells injected along the "
     "near-wall band r in [1.20, 1.45] mm. Panel (a): efficiency vs "
     "v_mean (0.01-0.50 m/s, log x) at fixed 200 pg loading — falls "
     "from ~100% at slow flow to ~40% at 0.5 m/s, with a transition "
     "near 0.03-0.05 m/s. Panel (b): efficiency vs loading (5-300 pg "
     "log x) at fixed v_mean = 0.05 m/s — rises from ~25% at 5 pg to "
     "~85% at 300 pg, with a clear threshold near 10-40 pg. Both "
     "panels have overlaid MCA reference velocity and Polyak loading "
     "reference lines for context."),
    ("fig21_static_vs_trajectory",
     "4.6  Fig 21 — Static vs trajectory (headline result)",
     "Stage 3c — the headline figure of the project. Two panels, both "
     "showing the effective capture range (um from the stent inner "
     "surface) as a function of a parameter. Solid blue is the static "
     "criterion |F_mag| > |F_drag|; dashed red is the trajectory "
     "effective range, computed by a 7-step binary search over the "
     "injection radius of a strut-aligned cell. Panel (a) sweeps SPION "
     "loading (10-200 pg) at v_mean = 0.2 m/s. Panel (b) sweeps velocity "
     "(0.02-0.5 m/s) at 50 pg. Headline numbers, verified in this "
     "session: at MCA conditions (50 pg, v = 0.2 m/s, B0 = 0.5 T) the "
     "trajectory predicts 96.9 um of effective capture range vs static "
     "6.8 um — a 14.2x extension factor. At v = 0.5 m/s the static "
     "criterion predicts ZERO capture but trajectory integration still "
     "finds 74.2 um of range. The ratio falls from 14.2x (50 pg) to "
     "3.7x (200 pg) as the static criterion strengthens. These numbers "
     "reproduce bit-identically across reruns because the entire "
     "pipeline is deterministic."),
]


def section_4(pdf: PdfPages):
    for stem, title, caption in FIGURE_PAGES:
        fig = _new_page(pdf)
        y = _heading(fig, TOP_Y, title, level=1)

        img_path = RESULTS / f"{stem}.png"
        if img_path.exists():
            img = mpimg.imread(img_path)
            # Place image in upper half (axes in figure coords)
            ax_img = fig.add_axes([0.08, 0.38, 0.84, 0.48])
            ax_img.imshow(img)
            ax_img.axis("off")
        else:
            fig.text(0.5, 0.60,
                     f"[Figure file not found: {stem}.png]\n"
                     f"Regenerate with: python -m stent_capture.figures.{stem}",
                     ha="center", va="center", fontsize=CAPTION_FS,
                     color="#a00", **MONO)

        _body(fig, 0.33, caption, fs=CAPTION_FS)
        _finish(pdf, fig)


# --- Section 5: Limitations -------------------------------------------------

def section_5(pdf: PdfPages):
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "5  Limitations and assumptions", level=1)
    intro = (
        "Every modelling framework trades fidelity for tractability. The "
        "table of trades made in stent_capture is explicit and exhaustive "
        "so that a dissertation examiner can see exactly what the results "
        "do and do not claim."
    )
    y = _body(fig, y, intro)
    items = [
        ("Single ring of struts",
         "The model has exactly 8 struts equally spaced at one axial "
         "position. Real stents are woven diamond-cell lattices with "
         "axial periodicity and azimuthal asymmetry. Adding axial "
         "periodicity is straightforward (just replicate the ring at "
         "multiple z) but has not been done because the single-ring "
         "analysis already isolates the relevant physics."),
        ("Straight cylindrical vessel",
         "Cerebral arteries are tortuous, bifurcated, and tapering. The "
         "model assumes a perfectly straight, uniform-radius cylinder. "
         "Curvature-induced secondary flows (Dean vortices) and the "
         "altered pressure gradient near bifurcations are not "
         "represented. This is acceptable for a local-region analysis "
         "near one stent but not for whole-vessel perfusion studies."),
        ("Cylindrical strut approximation in capture event",
         "The strut-proximity event uses a circumscribed cylinder of "
         "radius max(w, t)/2 = 50 um. For cells approaching along the "
         "through-strut axis the wall event fires first, so the "
         "cylindrical approximation has no effect on fig21. For off-axis "
         "cells in fig19 / fig20 it can over-classify some grazing "
         "trajectories as captured when a point-to-rectangle metric "
         "would let them escape."),
        ("No Brownian motion",
         "For a 10 um cell the thermal diffusion coefficient "
         "D_thermal = kT/(6 pi eta R_cell) ~ 2e-14 m^2/s. The "
         "displacement over the ~100 ms transit time is sqrt(2 D t) "
         "~ 60 nm, which is negligible compared to both the lumen "
         "dimensions (mm) and the capture range (um-tens of um). "
         "Brownian contributions to trajectory spread are therefore "
         "excluded."),
        ("No cell-cell interactions",
         "Trajectories are integrated independently. In reality, "
         "captured cells at the stent surface could shadow later cells "
         "or alter the local field (ferrofluid-like cooperative "
         "capture). The independent-trajectory approximation is exact "
         "in the dilute limit."),
    ]
    for name, desc in items:
        y = _body(fig, y, name, bold=True)
        y = _body(fig, y, desc)
        y -= 0.004
    _finish(pdf, fig)

    # Page 2 of limitations
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "5  Limitations (cont.)", level=1)
    more = [
        ("No rotational dynamics",
         "Cells are treated as non-rotating spheres. Real cells can "
         "rotate under shear, and SPION distributions within a rotating "
         "cell produce a time-varying magnetic moment. This coupling is "
         "not represented and is unlikely to matter for capture "
         "prediction but could matter for retention-under-shear studies."),
        ("Saturation magnetisation assumption",
         "The stent's M is passed directly as a saturation value when "
         "the assume_saturation flag is set. A full M(H) hysteresis "
         "curve is not implemented. At B0 >= ~0.1 T the assumption is "
         "physically reasonable for 304 / 430 / 2205 steels; at weaker "
         "fields it over-estimates the stent's contribution."),
        ("Linear-chi SPION force at high B0",
         "The Furlani-Ng scalar-gradient force assumes SPIONs are in "
         "the linear-susceptibility regime M_p = chi_p H. Real SPIONs "
         "saturate at |B| ~ 0.3-0.5 T. At B0 = 0.5 T the particles are "
         "near or at saturation, and the force should transition to "
         "F = V_p mu_0 M_s grad|H|, which scales as grad|B| alone. The "
         "current model thus over-estimates absolute force magnitudes "
         "by a factor of ~(B_total/B_sat) ~ 1-2 at the high-B0 end of "
         "the sweeps. Relative comparisons between static and "
         "trajectory predictions are unaffected because both branches "
         "use the same force model. This is the most physically "
         "significant limitation of the framework."),
        ("Steady (non-pulsatile) flow",
         "Cerebral arteries pulsate. The Womersley number alpha ~ 2 "
         "places the flow in the quasi-steady regime, but time-averaged "
         "efficiency could differ from steady-state efficiency by "
         "~10-20% in a more complete treatment."),
        ("Mislabelled literature references",
         "Two figure labels do not match their cited sources: "
         "'Polyak 2008 default = 10 pg' (actual Polyak working dose is "
         "200 pg) and 'MCA mean 0.2 m/s (Aaslid 1982)' (actual Aaslid "
         "value is 0.62 m/s). See AUDIT.md for the full discussion and "
         "the recommended relabellings. These do not change any "
         "numerical result — only the axis labels and docstrings need "
         "correction."),
        ("Orphaned shear_stress.py module",
         "A shear_stress.py file was added recently containing broken "
         "classes (int(length_in_metres) = 0 loops) with no citations "
         "and not imported anywhere. The AUDIT recommends deletion. "
         "BloodFlow.wall_shear_stress already provides the correct "
         "tau_w value used by every figure that needs it."),
    ]
    for name, desc in more:
        y = _body(fig, y, name, bold=True)
        y = _body(fig, y, desc)
        y -= 0.004
    _finish(pdf, fig)


# --- Section 6: How to use --------------------------------------------------

def section_6(pdf: PdfPages):
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "6  How to use the code", level=1)
    intro = (
        "The project is a standard Python package with no installation "
        "step beyond installing its three runtime dependencies "
        "(numpy, scipy, matplotlib). All scripts are run from the "
        "project root."
    )
    y = _body(fig, y, intro)
    y = _heading(fig, y, "6.1  Installation", level=2)
    inst = (
        "    pip install numpy scipy matplotlib pytest\n"
        "    cd stent_capture_project\n"
    )
    fig.text(MARGIN_L, y, inst, ha="left", va="top",
             fontsize=CODE_FS, **MONO)
    y -= 0.08
    y = _heading(fig, y, "6.2  Running the test suite", level=2)
    tests = (
        "    python -m pytest stent_capture/tests/ -v\n"
        "\n"
        "All 70 tests should pass in ~2 minutes. Any regression in the\n"
        "field kernel, force formula, or capture criterion will surface\n"
        "as a failing assertion with a clear message.\n"
    )
    fig.text(MARGIN_L, y, tests, ha="left", va="top",
             fontsize=CODE_FS, **MONO)
    y -= 0.14
    y = _heading(fig, y, "6.3  Regenerating individual figures", level=2)
    figs = (
        "    python -m stent_capture.figures.fig01_single_strut\n"
        "    python -m stent_capture.figures.fig12_external_field_comparison\n"
        "    python -m stent_capture.figures.fig17_spion_loading_sweep\n"
        "    python -m stent_capture.figures.fig21_static_vs_trajectory\n"
        "\n"
        "Each figure writes PNG + PDF output to results/. The Stage 1/2\n"
        "figures run in a few seconds each; fig18 takes ~15 s; fig20 is\n"
        "the expensive one (320 trajectories via multiprocessing.Pool,\n"
        "typically 3-6 minutes); fig21 is the headline result (90\n"
        "trajectories plus 10 static distance evaluations, ~180 s).\n"
    )
    fig.text(MARGIN_L, y, figs, ha="left", va="top",
             fontsize=CODE_FS, **MONO)
    _finish(pdf, fig)

    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "6  How to use the code (cont.)", level=1)
    y = _heading(fig, y, "6.4  Running parameter sweeps", level=2)
    intro = (
        "The capture_efficiency module exposes two ready-made sweep "
        "helpers. To sweep over velocities:"
    )
    y = _body(fig, y, intro)
    code = (
        "    from stent_capture.figures.common import make_ring\n"
        "    from stent_capture.physics.external_field import (\n"
        "        TotalField, UniformExternalField)\n"
        "    from stent_capture.physics.magnetic_force import SPIONLabelledCell\n"
        "    from stent_capture.simulation.capture_efficiency import (\n"
        "        capture_efficiency_vs_velocity)\n"
        "    import numpy as np\n"
        "\n"
        "    ring = make_ring()\n"
        "    ring.assume_saturation = True\n"
        "    tf   = TotalField(ring, UniformExternalField([0, 0, 0.5]))\n"
        "    cell = SPIONLabelledCell(spion_mass_per_cell=200e-15)\n"
        "\n"
        "    result = capture_efficiency_vs_velocity(\n"
        "        cell, tf, ring,\n"
        "        velocities=np.geomspace(0.01, 0.5, 8),\n"
        "        line_start=np.array([1.20e-3, 0, -2e-3]),\n"
        "        line_end  =np.array([1.45e-3, 0, -2e-3]),\n"
        "        n_cells_per_velocity=20,\n"
        "        z_end=2e-3, max_time=1.0,\n"
        "    )\n"
        "    print(result['velocities'], result['efficiencies'])\n"
    )
    fig.text(MARGIN_L, y, code, ha="left", va="top",
             fontsize=CODE_FS, **MONO)
    y -= 0.44
    y = _heading(fig, y, "6.5  Output file layout", level=2)
    out = (
        "results/figXX_name.png   — 300 dpi PNG for presentations\n"
        "results/figXX_name.pdf   — vector PDF for the dissertation\n"
        "\n"
        "Both files are produced together by the figure scripts'\n"
        "save_fig() helper. No other output locations are used.\n"
    )
    fig.text(MARGIN_L, y, out, ha="left", va="top",
             fontsize=CODE_FS, **MONO)
    _finish(pdf, fig)


# --- Section 7: References --------------------------------------------------

def section_7(pdf: PdfPages):
    fig = _new_page(pdf)
    y = _heading(fig, TOP_Y, "7  References", level=1)
    refs = [
        "Aaslid, R., Markwalder, T.-M., & Nornes, H. (1982). Noninvasive "
        "transcranial Doppler ultrasound recording of flow velocity in "
        "basal cerebral arteries. Journal of Neurosurgery, 57(6), 769-774.",

        "Akoun, G., & Yonnet, J.-P. (1984). 3D analytical calculation of "
        "the forces exerted between two cuboidal magnets. IEEE "
        "Transactions on Magnetics, 20(5), 1962-1964.",

        "Avilés, M. O., Ebner, A. D., & Ritter, J. A. (2007). Theoretical "
        "analysis of a transdermal ferromagnetic implant for retention "
        "of magnetic drug carrier particles. Journal of Magnetism and "
        "Magnetic Materials, 310(2), 428-443.",

        "Chorny, M., Fishbein, I., Yellen, B. B., Alferiev, I. S., "
        "Bakay, M., Ganta, S., Adamo, R., Amiji, M., Friedman, G., & "
        "Levy, R. J. (2007). Targeting stents with locally delivered "
        "paclitaxel-loaded magnetic nanoparticles prevents neointima "
        "formation in porcine arteries. FASEB J, 21(13), 3375-3383.",

        "Furlani, E. P. (2001). Permanent Magnet and Electromechanical "
        "Devices: Materials, Analysis, and Applications. Academic Press.",

        "Furlani, E. P., & Ng, K. C. (2006). Analytical model of magnetic "
        "nanoparticle transport and capture in the microvasculature. "
        "Physical Review E, 73(6), 061919.",

        "Pedley, T. J. (1980). The Fluid Mechanics of Large Blood Vessels. "
        "Cambridge University Press.",

        "Polyak, B., Fishbein, I., Chorny, M., Alferiev, I., Williams, "
        "D., Yellen, B., Friedman, G., & Levy, R. J. (2008). High field "
        "gradient targeting of magnetic nanoparticle-loaded endothelial "
        "cells to the surfaces of steel stents. Proceedings of the "
        "National Academy of Sciences, 105(2), 698-703.",

        "Riegler, J., Allain, B., Southern, P., Wells, J. A., Gauberti, "
        "M., Kaempfl, S., Lythgoe, M. F., Pankhurst, Q. A., & Davies, K. "
        "E. (2011). Superparamagnetic iron oxide nanoparticle targeting "
        "of MSCs in vascular injury. Biomaterials, 32(7), 1942-1951.",

        "Tefft, B. J., Uthamaraj, S., Harburn, J. J., Hlinomaz, O., "
        "Lerman, A., Dragomir-Daescu, D., & Sandhu, G. S. (2014). Design "
        "and Validation of a Novel Ferromagnetic Bare Metal Stent "
        "Capable of Capturing and Retaining Endothelial Cells. Annals "
        "of Biomedical Engineering, 42(11), 2269-2280.",

        "Tefft, B. J., Uthamaraj, S., Harburn, J. J., Klabusay, M., "
        "Dragomir-Daescu, D., & Sandhu, G. S. (2017). Magnetizable "
        "stent-grafts enable endothelial cell capture. Journal of "
        "Magnetism and Magnetic Materials, 427, 100-104.",

        "Womersley, J. R. (1955). Method for the calculation of velocity, "
        "rate of flow and viscous drag in arteries when the pressure "
        "gradient is known. Journal of Physiology, 127(3), 553-563.",
    ]
    for r in refs:
        y = _body(fig, y, r, fs=CAPTION_FS)
        y -= 0.002
    _finish(pdf, fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"  Building Overview.pdf ...")
    with PdfPages(OUT_PDF) as pdf:
        page_title(pdf)
        page_toc(pdf)
        section_1(pdf)
        section_2(pdf)
        section_3(pdf)
        section_4(pdf)
        section_5(pdf)
        section_6(pdf)
        section_7(pdf)

    import os
    sz = os.path.getsize(OUT_PDF) / 1024
    print(f"  [OK] {OUT_PDF.name} written ({_page_no[0]} pages, {sz:.0f} KB)")


if __name__ == "__main__":
    main()
