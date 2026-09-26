# Taylor–Green quantitative convergence plan

The first numerical gate of the project is not visual agreement. It is quantitative convergence against an exact solution.

## Why Taylor–Green comes first

The 2-D Taylor–Green vortex used here is smooth, doubly periodic, divergence-free, and has a closed-form viscous decay. That makes it suitable for checking the periodic spectral operators and time integrator before introducing walls, twisted boundary identifications, D2/D3 operators, Gray schedules, or electrokinetics.

The reference field is regenerated analytically for every run. No remembered or context-derived numerical value is used as ground truth.

## Stage TG-1 — temporal convergence

The current solver uses classical RK4. For the single-mode Taylor–Green field, the spatial Fourier representation is exact to roundoff on a modest grid, so holding the grid fixed and reducing `dt` isolates the temporal error.

For consecutive time steps `dt_c > dt_f` and corresponding velocity errors `E_c` and `E_f`, the observed order is

```math
p=\frac{\log(E_c/E_f)}{\log(dt_c/dt_f)}.
```

For RK4 we expect the asymptotic temporal order to approach four before floating-point roundoff dominates.

The automated test therefore checks:

1. incompressibility remains at roundoff level;
2. the exact-solution error decreases under refinement;
3. the measured order is consistent with fourth-order time integration;
4. every requested `dt` lies inside a conservative explicit-diffusion RK4 stability bound.

### Why the diffusion stability guard matters

A useful failure was found while setting up the first refinement sequence. A Taylor–Green initial condition occupies only a low Fourier mode, so an overly large `dt` can appear harmless at first. However, floating-point roundoff seeds high-frequency modes. For explicit RK4, viscous diffusion contributes eigenvalues

```math
\lambda_k=-\nu |k|^2.
```

Classical RK4 is stable along the negative real axis only out to approximately

```math
|\lambda_k\,dt|\lesssim 2.785.
```

Therefore a necessary viscous stability bound is

```math
dt \lesssim \frac{2.785}{\nu\,k_{max}^2},
```

with a safety factor applied in code. This is not the full CFL condition: nonlinear advection may require a smaller step. But it prevents a convergence study from accidentally comparing unstable solutions and interpreting their growth as discretization error.

The solver exposes this check as `PeriodicVorticitySolver2D.diffusive_rk4_dt_limit()`, and the convergence runner rejects unsafe time steps by default.

The acceptance band in the regression test is intentionally broad. It is meant to detect a broken convergence rate, not to encode a machine-specific floating-point result.

## Stage TG-2 — spatial convergence

Taylor–Green's single Fourier mode is *not* a useful test of spectral spatial convergence once it is resolved, because refinement then changes almost nothing except roundoff.

For that reason we will not manufacture a fake spatial convergence slope from Taylor–Green. The next spatial-order check will use either:

- a smooth non-band-limited periodic analytic field to verify the spectral differential operators, followed by
- a published/non-trivial flow benchmark whose PDE solution has measurable spatial truncation error.

The first item checks the spatial machinery itself; the second checks the complete flow solver. This distinction is deliberate: a convergence study must refine an error source that actually exists.

## Stage TG-3 — invariant and balance checks

Alongside pointwise/field error, the solver records:

- RMS divergence;
- kinetic energy;
- enstrophy;
- mean vorticity/circulation density.

For unforced inviscid Taylor–Green, the 2-D single mode is stationary under the nonlinear term and should preserve energy up to numerical error. With viscosity, its decay is known analytically.

## Reproducible command

Run

```bash
python scripts/run_taylor_green_convergence.py
```

to print the current refinement table, the conservative diffusion stability limit, relative velocity errors, RMS divergence, and observed pairwise orders. The script recomputes every quantity from the current code; no result table is hard-coded into the repository.

## Source provenance

See [`REFERENCES.md`](REFERENCES.md) for the original Taylor–Green citation and the modern verification source used to document the exact 2-D specialization.

The repository rule remains:

> If a comparison needs a specific published number, that number must be traceable to a source, coordinate convention, nondimensionalization, and transformation used by this code.

Taylor–Green TG-1 avoids this issue entirely by regenerating the exact field from the analytical formula.
