"""Run the quantitative Taylor--Green temporal-convergence study.

This script is intentionally small and reproducible. It does not embed
published comparison numbers. For every run, the numerical velocity is
compared with the closed-form Taylor--Green solution evaluated at the same
final time.

Usage
-----

    python scripts/run_taylor_green_convergence.py

The default sequence halves ``dt`` three times. The chosen steps lie inside
the explicit RK4 viscous-stability bound for the configured grid and viscosity.
The resulting table can be copied directly into a report or issue without
hand-calculating errors or observed orders.
"""

from __future__ import annotations

import numpy as np

from navier_stokes_gray.benchmarks import taylor_green_vorticity
from navier_stokes_gray.convergence import observed_orders, temporal_convergence_study
from navier_stokes_gray.solver2d import PeriodicVorticitySolver2D
from navier_stokes_gray.spectral import make_grid


def main() -> None:
    """Execute the standard temporal-refinement study and print a compact table."""

    # These are numerical experiment settings, not literature benchmark values.
    # The exact reference field is regenerated analytically by the package.
    dts = [0.008, 0.004, 0.002, 0.001]
    resolution = 32
    nu = 0.5
    final_time = 0.4

    # Print the explicit diffusion bound alongside the requested study so a
    # future parameter change cannot silently turn the refinement sequence into
    # an unstable run. This bound is necessary, though not sufficient, for the
    # full nonlinear solver because convective CFL may impose a smaller step.
    grid = make_grid(resolution)
    x, y = np.meshgrid(grid.x, grid.y, indexing="ij")
    omega0 = taylor_green_vorticity(x, y, 0.0, nu=nu)
    stability_solver = PeriodicVorticitySolver2D(grid=grid, nu=nu, omega=omega0)
    diffusive_limit = stability_solver.diffusive_rk4_dt_limit()

    points = temporal_convergence_study(
        dts,
        resolution=resolution,
        nu=nu,
        final_time=final_time,
    )
    rates = observed_orders(points)

    print("Taylor--Green temporal convergence")
    print(f"resolution={resolution}, nu={nu}, final_time={final_time}")
    print(f"conservative RK4 diffusion dt limit={diffusive_limit:.8g}")
    print()
    print(f"{'dt':>10} {'relative L2 velocity error':>28} {'rms(div u)':>16} {'order':>10}")
    print("-" * 70)

    for index, point in enumerate(points):
        order_text = "-" if index == 0 else f"{rates[index - 1].order:.6f}"
        print(
            f"{point.dt:10.5g} "
            f"{point.velocity_relative_l2:28.12e} "
            f"{point.rms_divergence:16.8e} "
            f"{order_text:>10}"
        )


if __name__ == "__main__":
    main()
