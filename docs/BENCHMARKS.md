# Validation benchmark ladder

This project should not introduce new topological operators before the numerical machinery reproduces problems whose expected behavior is already known.

The benchmark ladder is therefore ordered from exact solutions to increasingly complex reference problems.

Every benchmark used here should satisfy two documentation rules:

1. record the published/original reference whenever one exists;
2. record the exact numerical formulation used by this repository, because a modern benchmark specialization may differ from the historical paper that introduced the flow family.

A consolidated bibliography is maintained in [`REFERENCES.md`](REFERENCES.md).

## B0 — Taylor–Green vortex (2-D periodic exact solution)

Purpose:

- verify spatial discretization;
- verify time integration;
- verify incompressibility;
- verify viscous decay;
- establish grid-convergence/error metrics;
- validate the periodic-domain machinery before any twisted identification is added.

Historical source:

- G. I. Taylor and A. E. Green, *Mechanism of the Production of Small Eddies from Large Ones*, Proc. Royal Society A 158 (1937), 499–521. DOI: https://doi.org/10.1098/rspa.1937.0036

Useful modern verification formulation:

- TensorMesh Taylor–Green convergence example: https://docs.tensor-mesh.com/example_gallery/fluid/taylor_green.html

The original 1937 paper studies the Taylor–Green vortex family; the exact 2-D doubly-periodic field used in this repository is a standard verification specialization.

For the incompressible Navier–Stokes equations on a doubly periodic domain, one standard 2-D Taylor–Green solution is

```math
u_x(x,y,t)=-U_0\cos(kx)\sin(ky)e^{-2\nu k^2 t},
```

```math
u_y(x,y,t)= U_0\sin(kx)\cos(ky)e^{-2\nu k^2 t}.
```

A compatible pressure field is

```math
p(x,y,t)=-\frac{\rho U_0^2}{4}\left[\cos(2kx)+\cos(2ky)\right]e^{-4\nu k^2t}+C.
```

The exact solution gives us a pointwise reference, so no qualitative judgement is needed.

Primary diagnostics:

```math
E_u=\frac{\|\mathbf u_{num}-\mathbf u_{exact}\|_2}
{\|\mathbf u_{exact}\|_2},
```

```math
E_{div}=\|\nabla\cdot\mathbf u\|_2,
```

and kinetic-energy decay.

Acceptance criterion for later stages: demonstrated convergence under mesh/time-step refinement.

---

## B1 — Lid-driven cavity

Purpose:

- verify the nonlinear convective term;
- verify no-slip wall handling;
- compare centerline velocity profiles and vortex locations with standard reference datasets;
- test moderate/high Reynolds-number behavior before topology is modified.

Canonical reference:

- U. Ghia, K. N. Ghia and C. T. Shin, *High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method*, Journal of Computational Physics 48 (1982), 387–411. DOI: https://doi.org/10.1016/0021-9991(82)90058-4

Additional high-accuracy tabulated benchmark data:

- AceNumerics lid-driven cavity benchmark tables: https://www.acenumerics.com/the-benchmarks.html

The cavity is deliberately conventional. If the solver does not reproduce this problem, no interpretation of D2/D3, Gray ordering, or topology is trustworthy.

Reference observables:

- `u(y)` on the vertical centerline;
- `v(x)` on the horizontal centerline;
- primary-vortex center;
- secondary-corner vortices;
- mass conservation and residual convergence.

---

## B2 — Low-Re reversible / journal-bearing style flow

Purpose:

- validate the Stokes-limit forward/reverse protocol;
- separate kinematic reversibility from diffusion and numerical irreversibility;
- introduce material tracers/passive scalar before turbulence;
- study complex Lagrangian trajectories at low Reynolds number.

Reference literature:

- P. Dutta and R. Chevray, *Inertial effects in chaotic mixing with diffusion*, Journal of Fluid Mechanics. The paper studies passive-scalar dispersion in a low-Re periodic journal-bearing flow and compares numerical and experimental observations: https://www.cambridge.org/core/product/identifier/S0022112095000437/type/journal_article
- For historical analytic Stokes solutions in eccentric-cylinder geometry, see the review discussion in *Chaotic advection and heat transfer enhancement in Stokes flows*: https://www.sciencedirect.com/science/article/abs/pii/S0142727X03000225

Protocol:

```text
forcing sequence:  S1, S2, ..., SN
reverse sequence:  SN^-1, ..., S2^-1, S1^-1
```

Measure

```math
E_R=\frac{\|c_{final}-c_{initial}\|_2}
{\|c_{initial}\|_2}.
```

This benchmark is important because complicated stretching/folding can occur without turbulence. That distinction is central to the project.

---

## B3 — Kolmogorov flow on a periodic 2-D domain

Purpose:

- establish a nonlinear periodic-flow benchmark;
- provide a natural baseline before replacing an ordinary periodic identification with an orientation-reversing one;
- study bifurcation/instability while keeping the domain free of solid-wall corner effects.

Reference examples:

- K. Seshasayanan, V. Dallas and S. Fauve, *Bifurcations of a plane parallel flow with Kolmogorov forcing*, Physical Review Fluids 6, 103902 (2021): https://doi.org/10.1103/PhysRevFluids.6.103902
- NASA/ICASE report on 3-D Kolmogorov flow with fully periodic boundaries: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960020961.pdf

This is the first benchmark that is structurally close to the later topological experiments.

Baseline:

```math
(x+L_x,y)\sim(x,y),\qquad (x,y+L_y)\sim(x,y).
```

Later controlled modification:

```math
(x,y+L_y)\sim(L_x-x,y).
```

Only one change should be introduced at a time.

---

## B4 — 3-D periodic reference flow (ABC-type)

Purpose:

- move from 2-D boundary identifications to genuinely 3-D transport;
- study braid-like Lagrangian structure;
- prepare the ground for D2/D3 realizations as spatial operators rather than drawings of specific vessels.

This stage comes only after B0–B3 are quantitatively stable. A specific ABC-flow reference will be fixed before implementation so the forcing convention, nondimensionalization and target diagnostics are unambiguous.

---

# Experimental logic

For each benchmark we maintain two layers:

```text
P0 = published/reference problem
P1 = same problem + exactly one project-specific modification
```

Examples:

```text
P0: ordinary periodic boundary
P1: one orientation-reversing boundary map
```

or

```text
P0: baseline forcing
P1: baseline forcing + one D2 operation
```

A Gray-code schedule can later order the binary feature set so that consecutive runs differ in one feature only.

# Documentation and code-comment rule

For every new solver, boundary map or benchmark implementation:

- comments should explain the physical/numerical reason for non-obvious operations, not merely restate the line of code;
- the module docstring should name the benchmark/reference it implements;
- the relevant DOI or stable source URL should appear either in the module docstring or in this benchmark document;
- assumptions and nondimensionalization must be written down before reference numbers are compared;
- no tabulated number should be copied into tests without recording its source and the precise convention used.

# What must not be inferred from a benchmark

Passing a benchmark only validates the part of the implementation that benchmark exercises. It does not validate a new physical claim.

In particular:

- Taylor–Green validates numerics, not useful separation;
- reversible Stokes tests validate reversibility, not energy extraction;
- a twisted boundary map validates a mathematically consistent identification, not a physical Klein-bottle advantage;
- any future ionic or electrochemical claim requires explicit Poisson–Nernst–Planck coupling and energy/electrochemical-potential balances.

The project should preserve this separation between **validation**, **controlled modification**, and **new hypothesis** throughout.