# References and provenance

This file records the published sources and stable benchmark links used by the project.

The intent is practical: every benchmark number, forcing convention, exact solution, or comparison profile used in code or tests should be traceable to a source and to a clearly stated nondimensionalization/convention.

## Taylor–Green vortex

### Historical source

G. I. Taylor and A. E. Green, **Mechanism of the Production of Small Eddies from Large Ones**, *Proceedings of the Royal Society A*, 158 (1937), 499–521.

DOI: https://doi.org/10.1098/rspa.1937.0036

The repository uses a standard 2-D doubly-periodic Taylor–Green specialization as an exact Navier–Stokes verification field. The historical paper establishes the Taylor–Green vortex family; the exact modern 2-D verification formula used in code is documented separately.

### Modern verification implementation

TensorMesh, **Taylor-Green Vortex (Convergence Study)**:

https://docs.tensor-mesh.com/example_gallery/fluid/taylor_green.html

This is useful as an independently implemented modern benchmark with an explicit closed-form 2-D solution and convergence study.

## Lid-driven cavity

U. Ghia, K. N. Ghia and C. T. Shin, **High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method**, *Journal of Computational Physics* 48 (1982), 387–411.

DOI: https://doi.org/10.1016/0021-9991(82)90058-4

This is the canonical comparison source for centerline velocity profiles and high-Re cavity structure.

Additional high-accuracy benchmark tables:

https://www.acenumerics.com/the-benchmarks.html

Any value copied into repository tests must identify whether it comes from Ghia et al., AceNumerics, or another source.

## Low-Reynolds-number journal-bearing / chaotic-advection references

P. Dutta and R. Chevray, **Inertial effects in chaotic mixing with diffusion**, *Journal of Fluid Mechanics*.

Stable Cambridge article page:

https://www.cambridge.org/core/product/identifier/S0022112095000437/type/journal_article

This work is relevant because it compares numerical and experimental passive-scalar dispersion in a low-Re periodic journal-bearing flow and discusses transient inertial and diffusive effects.

Review containing the journal-bearing Stokes-flow background and references to analytic eccentric-cylinder solutions:

**Chaotic advection and heat transfer enhancement in Stokes flows**:

https://www.sciencedirect.com/science/article/abs/pii/S0142727X03000225

The project will select one precise eccentric-cylinder protocol before implementing this benchmark, rather than mixing data from different rotation schedules or geometries.

## Kolmogorov flow

K. Seshasayanan, V. Dallas and S. Fauve, **Bifurcations of a plane parallel flow with Kolmogorov forcing**, *Physical Review Fluids* 6, 103902 (2021).

DOI: https://doi.org/10.1103/PhysRevFluids.6.103902

For a fully periodic three-dimensional Kolmogorov-flow reference:

J. V. Shebalin and S. L. Woodruff, **Kolmogorov Flow in Three Dimensions**, NASA/ICASE report (1996):

https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960020961.pdf

## Source-handling rule

When implementing a benchmark:

1. write down the governing equations and nondimensionalization used by this repository;
2. cite the source of the reference solution or tabulated data;
3. state any coordinate/sign/pressure-gauge conversion explicitly;
4. put the source near the relevant code when the connection is not obvious;
5. do not treat a secondary implementation as the original source when the original paper is available;
6. do not compare numbers from two references until their conventions are confirmed compatible.

The goal is reproducibility, not merely obtaining visually similar plots.
