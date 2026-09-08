# Validation benchmark ladder

This project should not introduce new topological operators before the numerical machinery reproduces problems whose expected behavior is already known.

The benchmark ladder is therefore ordered from exact solutions to increasingly complex reference problems.

## B0 — Taylor–Green vortex (2-D periodic exact solution)

Purpose:

- verify spatial discretization;
- verify time integration;
- verify incompressibility;
- verify viscous decay;
- establish grid-convergence/error metrics;
- validate the periodic-domain machinery before any twisted identification is added.

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

This stage comes only after B0–B3 are quantitatively stable.

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

# What must not be inferred from a benchmark

Passing a benchmark only validates the part of the implementation that benchmark exercises. It does not validate a new physical claim.

In particular:

- Taylor–Green validates numerics, not useful separation;
- reversible Stokes tests validate reversibility, not energy extraction;
- a twisted boundary map validates a mathematically consistent identification, not a physical Klein-bottle advantage;
- any future ionic or electrochemical claim requires explicit Poisson–Nernst–Planck coupling and energy/electrochemical-potential balances.

The project should preserve this separation between **validation**, **controlled modification**, and **new hypothesis** throughout.