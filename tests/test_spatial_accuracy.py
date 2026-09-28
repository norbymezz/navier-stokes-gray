import numpy as np

from navier_stokes_gray.spatial_accuracy import spatial_accuracy_study


def test_spectral_derivative_errors_drop_under_refinement():
    """The non-band-limited analytic field must show genuine spatial refinement.

    Unlike the single-mode Taylor--Green vortex, this manufactured periodic
    field contains infinitely many Fourier modes.  Increasing resolution should
    therefore reduce truncation error rapidly until roundoff dominates.
    """

    points = spatial_accuracy_study([8, 16, 24, 32])

    ex = np.array([p.derivative_x_relative_l2 for p in points])
    ey = np.array([p.derivative_y_relative_l2 for p in points])
    elap = np.array([p.laplacian_relative_l2 for p in points])

    assert np.all(np.diff(ex) < 0.0)
    assert np.all(np.diff(ey) < 0.0)
    assert np.all(np.diff(elap) < 0.0)

    # These are deliberately loose regression gates, not reference-table
    # values.  Their role is to catch a damaged spectral implementation while
    # avoiding dependence on one machine's roundoff floor.
    assert points[-1].derivative_x_relative_l2 < 1e-10
    assert points[-1].derivative_y_relative_l2 < 1e-8
    assert points[-1].laplacian_relative_l2 < 1e-8
