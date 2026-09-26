# Lid-driven cavity validation

This stage introduces the first wall-bounded nonlinear benchmark after the
periodic Taylor--Green checks.

## Reference problem

We use the unit-square lid-driven cavity at `Re = 100` with a unit horizontal
lid velocity.  The numerical solution is compared with the centerline velocity
samples attributed to Ghia, Ghia & Shin (1982).  The exact values used by the
code and their provenance are stored in `src/navier_stokes_gray/cavity_reference.py`.

Primary source:

- U. Ghia, K. N. Ghia, C. T. Shin (1982), *High-Re solutions for incompressible
  flow using the Navier-Stokes equations and a multigrid method*, Journal of
  Computational Physics 48, 387--411.
- DOI: https://doi.org/10.1016/0021-9991(82)90058-4

Independent public transcription used to cross-check the entered values:

- Sandia SIERRA/Fuego Verification Manual, lid-driven-cavity benchmark,
  Table 22.1-1: https://www.sandia.gov/app/uploads/sites/315/2025/11/Fuego_Verification_5_18.pdf

## What is compared

For the computed solution we extract

```math
u(0.5,y)
```

and

```math
v(x,0.5).
```

The numerical centerlines are interpolated to the *published sample
coordinates*.  We do not substitute the nearest grid node, because the Ghia
coordinates are nonuniform.

For each profile the report records:

```math
\mathrm{RMSE}=\sqrt{\frac{1}{N}\sum_i (q_i-q_i^{ref})^2},
```

```math
E_\infty=\max_i |q_i-q_i^{ref}|,
```

and

```math
E_2=\frac{\|q-q^{ref}\|_2}{\|q^{ref}\|_2}.
```

The cavity comparison additionally records the pseudo-time convergence update
and an RMS divergence diagnostic.

## Important validation rule

The current baseline cavity solver is intentionally simple: central finite
differences, a streamfunction-vorticity formulation, Jacobi Poisson iteration,
and explicit pseudo-time marching.  Therefore we do **not** freeze a Ghia error
threshold into CI before observing grid/refinement behavior.

The validation sequence is:

1. verify the reference arrays, conventions and interpolation path;
2. require a cheap finite/stability smoke test in CI;
3. run the full `Re=100` comparison and record the profile errors;
4. repeat at increasing odd grid resolutions;
5. only after an asymptotic trend is visible, choose a documented regression
   threshold.

This prevents a coarse-grid result from being declared "validated" merely
because an arbitrary tolerance happened to pass.

## Separation from later topology experiments

No Gray ordering, D2/D3 operator, twisted identification or braid operation is
allowed into this benchmark.  The purpose of this stage is to establish how
the baseline nonlinear wall-bounded solver behaves against published data.
Only after that baseline is understood do we modify one feature at a time.
