import numpy as np
import pytest

from navier_stokes_gray.convergence import (
    ConvergencePoint,
    observed_orders,
    run_taylor_green_case,
    temporal_convergence_study,
)


def test_single_taylor_green_case_remains_divergence_free():
    """The convergence runner must preserve the solver's incompressibility."""

    point = run_taylor_green_case(
        resolution=32,
        nu=0.05,
        final_time=0.2,
        dt=0.01,
    )

    assert point.velocity_relative_l2 < 1e-7
    assert point.rms_divergence < 1e-12


def test_temporal_refinement_reduces_taylor_green_error():
    """Halving dt should monotonically reduce error before roundoff dominates.

    These are not remembered benchmark values. Every error is computed against
    the exact Taylor--Green velocity regenerated at the requested final time.
    The refinement sequence is chosen to remain inside the solver's conservative
    RK4 diffusion-stability limit for this grid and viscosity.
    """

    points = temporal_convergence_study(
        [0.008, 0.004, 0.002],
        resolution=32,
        nu=0.5,
        final_time=0.4,
    )
    errors = np.array([p.velocity_relative_l2 for p in points])

    assert np.all(np.diff(errors) < 0.0)


def test_observed_temporal_order_is_consistent_with_rk4():
    """Verify fourth-order temporal convergence without using external tables.

    The Taylor--Green mode is represented exactly in space on this grid, so the
    dominant refinement error is the RK4 integration of exponential viscous
    decay. The observed order should therefore approach four until floating-
    point roundoff becomes important.
    """

    points = temporal_convergence_study(
        [0.008, 0.004, 0.002],
        resolution=32,
        nu=0.5,
        final_time=0.4,
    )
    rates = observed_orders(points)

    # Use a deliberately broad numerical acceptance band. The purpose is to
    # detect loss of RK4 convergence, not to hard-code a particular machine's
    # floating-point result.
    for rate in rates:
        assert 3.5 < rate.order < 4.5


def test_unstable_explicit_diffusion_step_is_rejected():
    """A refinement study must fail loudly before entering RK4 instability.

    Even a low-mode Taylor--Green initial condition can acquire high Fourier
    components at roundoff level. Explicit viscous integration outside RK4's
    negative-real-axis stability interval can amplify those modes and produce
    meaningless 'convergence' data. The runner therefore checks a conservative
    diffusion bound before evolving the solution.
    """

    with pytest.raises(ValueError, match="diffusion"):
        run_taylor_green_case(
            resolution=32,
            nu=0.5,
            final_time=0.4,
            dt=0.02,
        )


def test_observed_order_formula_with_synthetic_exact_ratio():
    """Unit-test the order estimator independently of the CFD solver."""

    points = [
        ConvergencePoint(32, 0.1, 1.0, 1.6e-4, 0.0),
        ConvergencePoint(32, 0.05, 1.0, 1.0e-5, 0.0),
    ]
    rate = observed_orders(points)[0]

    # Error ratio 16 under a factor-two refinement corresponds exactly to p=4.
    assert np.isclose(rate.order, 4.0)
