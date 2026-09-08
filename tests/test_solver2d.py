import numpy as np

from navier_stokes_gray.benchmarks import (
    taylor_green_velocity,
    taylor_green_vorticity,
    vector_relative_l2,
)
from navier_stokes_gray.solver2d import PeriodicVorticitySolver2D
from navier_stokes_gray.spectral import make_grid


def _mesh(grid):
    return np.meshgrid(grid.x, grid.y, indexing="ij")


def test_taylor_green_velocity_is_recovered_from_vorticity():
    grid = make_grid(32)
    x, y = _mesh(grid)
    omega0 = taylor_green_vorticity(x, y, 0.0, nu=0.02)
    solver = PeriodicVorticitySolver2D(grid=grid, nu=0.02, omega=omega0)

    numerical_velocity = solver.velocity()
    exact_velocity = taylor_green_velocity(x, y, 0.0, nu=0.02)

    assert vector_relative_l2(exact_velocity, numerical_velocity) < 1e-12


def test_taylor_green_remains_divergence_free():
    grid = make_grid(32)
    x, y = _mesh(grid)
    omega0 = taylor_green_vorticity(x, y, 0.0, nu=0.01)
    solver = PeriodicVorticitySolver2D(grid=grid, nu=0.01, omega=omega0)

    assert solver.diagnostics().rms_divergence < 1e-12


def test_taylor_green_viscous_decay_matches_exact_solution():
    """Numerical evolution is compared with the exact field, not tabulated memory.

    Taylor--Green is especially useful here because the exact time dependence is
    known analytically.  This test therefore does not depend on copied benchmark
    numbers: it recomputes the reference field from the published solution.
    """

    grid = make_grid(32)
    x, y = _mesh(grid)
    nu = 0.05
    final_time = 0.2
    omega0 = taylor_green_vorticity(x, y, 0.0, nu=nu)
    solver = PeriodicVorticitySolver2D(grid=grid, nu=nu, omega=omega0)

    # The solution is smooth and single-mode, so this moderate step is enough
    # for a regression test while keeping the test suite fast.
    solver.run(final_time=final_time, dt=0.005)

    numerical_velocity = solver.velocity()
    exact_velocity = taylor_green_velocity(x, y, final_time, nu=nu)
    error = vector_relative_l2(exact_velocity, numerical_velocity)

    assert error < 1e-9


def test_unforced_inviscid_taylor_green_preserves_energy_over_short_run():
    grid = make_grid(32)
    x, y = _mesh(grid)
    omega0 = taylor_green_vorticity(x, y, 0.0, nu=0.0)
    solver = PeriodicVorticitySolver2D(grid=grid, nu=0.0, omega=omega0)

    initial_energy = solver.diagnostics().kinetic_energy
    solver.run(final_time=0.1, dt=0.005)
    final_energy = solver.diagnostics().kinetic_energy

    # Taylor--Green's nonlinear Jacobian vanishes for this 2-D mode, so the
    # inviscid field should remain stationary apart from roundoff/RK error.
    assert abs(final_energy - initial_energy) < 1e-12
