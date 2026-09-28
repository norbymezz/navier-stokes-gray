import numpy as np

from navier_stokes_gray.benchmarks import (
    kinetic_energy_density,
    relative_l2,
    reversibility_error,
    scalar_variance,
    taylor_green_pressure,
    taylor_green_velocity,
    taylor_green_vorticity,
    vector_relative_l2,
)


def test_taylor_green_is_divergence_free_analytically_on_grid():
    # Evaluate the analytic derivatives explicitly on a regular grid.
    n = 64
    x = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    y = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    xx, yy = np.meshgrid(x, y, indexing="ij")

    nu = 0.02
    t = 0.7
    u0 = 1.3
    k = 1.0
    decay = np.exp(-2.0 * nu * k * k * t)

    du_dx = u0 * k * np.sin(k * xx) * np.sin(k * yy) * decay
    dv_dy = -u0 * k * np.sin(k * xx) * np.sin(k * yy) * decay

    assert np.max(np.abs(du_dx + dv_dy)) < 1e-14


def test_taylor_green_velocity_decay_factor():
    x = np.array([[0.4]])
    y = np.array([[1.1]])
    nu = 0.1
    dt = 0.6

    u0_field = taylor_green_velocity(x, y, 0.0, nu=nu)
    u1_field = taylor_green_velocity(x, y, dt, nu=nu)

    expected = np.exp(-2.0 * nu * dt)
    assert np.allclose(u1_field[0], u0_field[0] * expected)
    assert np.allclose(u1_field[1], u0_field[1] * expected)


def test_taylor_green_pressure_has_expected_decay():
    x = np.array([[0.2]])
    y = np.array([[0.8]])
    nu = 0.05
    dt = 0.4

    p0 = taylor_green_pressure(x, y, 0.0, nu=nu)
    p1 = taylor_green_pressure(x, y, dt, nu=nu)

    expected = np.exp(-4.0 * nu * dt)
    assert np.allclose(p1, p0 * expected)


def test_vorticity_matches_definition_at_simple_point():
    x = np.array([[0.0]])
    y = np.array([[0.0]])
    omega = taylor_green_vorticity(x, y, 0.0, u0=2.0, k=3.0)
    assert np.allclose(omega, 12.0)


def test_error_metrics_are_zero_for_identical_fields():
    a = np.arange(9, dtype=float).reshape(3, 3)
    assert relative_l2(a, a) == 0.0
    assert vector_relative_l2((a, a), (a, a)) == 0.0
    assert reversibility_error(a, a) == 0.0


def test_scalar_variance_and_energy_are_nonnegative():
    c = np.array([0.0, 1.0, 2.0, 3.0])
    assert scalar_variance(c) >= 0.0

    u = np.array([1.0, -2.0])
    v = np.array([3.0, 4.0])
    e = kinetic_energy_density(u, v)
    assert np.all(e >= 0.0)
