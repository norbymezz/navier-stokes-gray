import numpy as np
import pytest

from navier_stokes_gray.cavity import LidDrivenCavity2D
from navier_stokes_gray.cavity_reference import ghia_re100
from navier_stokes_gray.cavity_validation import _profile_error, compare_centerlines


def test_profile_error_matches_hand_calculation():
    reference = np.array([1.0, -1.0, 2.0])
    numerical = np.array([1.5, -1.5, 1.0])
    result = _profile_error(reference, numerical)

    delta = numerical - reference
    assert np.isclose(result.rmse, np.sqrt(np.mean(delta * delta)))
    assert np.isclose(result.max_abs, 1.0)
    assert np.isclose(
        result.relative_l2,
        np.linalg.norm(delta) / np.linalg.norm(reference),
    )


def test_profile_error_rejects_zero_reference_norm():
    with pytest.raises(ValueError, match="non-zero"):
        _profile_error(np.zeros(3), np.ones(3))


def test_ghia_reference_is_monotone_in_coordinates_and_complete():
    """Catch transcription/order mistakes before comparing a CFD result."""

    ref = ghia_re100()
    assert ref.reynolds == 100.0
    assert len(ref.y_for_u) == len(ref.u_centerline) == 17
    assert len(ref.x_for_v) == len(ref.v_centerline) == 17
    assert np.all(np.diff(ref.y_for_u) > 0.0)
    assert np.all(np.diff(ref.x_for_v) > 0.0)
    assert ref.y_for_u[0] == 0.0 and ref.y_for_u[-1] == 1.0
    assert ref.x_for_v[0] == 0.0 and ref.x_for_v[-1] == 1.0
    assert ref.u_centerline[0] == 0.0 and ref.u_centerline[-1] == 1.0
    assert ref.v_centerline[0] == 0.0 and ref.v_centerline[-1] == 0.0


def test_compare_centerlines_rejects_wrong_reynolds_number():
    solver = LidDrivenCavity2D(n=17, reynolds=50.0, dt=1.0e-4)
    with pytest.raises(ValueError, match="Reynolds"):
        compare_centerlines(solver, ghia_re100())


def test_short_cavity_march_stays_finite():
    """Cheap CI smoke test: catch immediate instability without claiming validation.

    This deliberately does not compare to Ghia.  A few pseudo-time steps are
    only a software/stability check; quantitative validation requires marching
    to steady state and reporting the published centerline errors.
    """

    solver = LidDrivenCavity2D(
        n=17,
        reynolds=100.0,
        dt=2.0e-4,
        poisson_iterations=20,
    )
    for _ in range(25):
        solver.step()

    assert np.all(np.isfinite(solver.psi))
    assert np.all(np.isfinite(solver.omega))
    assert np.all(np.isfinite(solver.u))
    assert np.all(np.isfinite(solver.v))
    assert solver.diagnostics().max_speed >= 1.0
