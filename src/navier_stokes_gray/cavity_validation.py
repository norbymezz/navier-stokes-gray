"""Validation helpers for the lid-driven-cavity benchmark.

The purpose of this module is to keep *comparison logic* separate from the
cavity solver itself.  The solver produces a field; this module samples that
field at the published Ghia centerline coordinates and reports errors against
reference values whose provenance is documented in ``cavity_reference.py``.

No acceptance threshold is hard-coded here.  Thresholds belong in benchmark
reports/tests and should be justified by grid resolution and solver accuracy.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cavity import LidDrivenCavity2D
from .cavity_reference import CavityCenterlineReference, ghia_re100


@dataclass(frozen=True)
class ProfileError:
    """Error metrics for one sampled velocity profile."""

    rmse: float
    max_abs: float
    relative_l2: float


@dataclass(frozen=True)
class CavityValidationResult:
    """Comparison of numerical centerline profiles with one reference case."""

    reynolds: float
    resolution: int
    iterations: int
    final_update_inf: float
    u_error: ProfileError
    v_error: ProfileError
    rms_divergence: float


def _profile_error(reference: np.ndarray, numerical: np.ndarray) -> ProfileError:
    """Return RMSE, max absolute error and relative L2 error.

    The relative norm is useful across grids, while RMSE and max error retain
    an intuitive velocity scale.  A zero reference norm is rejected because a
    relative error would then be undefined.
    """

    reference = np.asarray(reference, dtype=float)
    numerical = np.asarray(numerical, dtype=float)
    if reference.shape != numerical.shape:
        raise ValueError("reference and numerical profiles must have equal shape")

    delta = numerical - reference
    reference_norm = float(np.linalg.norm(reference))
    if reference_norm == 0.0:
        raise ValueError("reference profile must have non-zero L2 norm")

    return ProfileError(
        rmse=float(np.sqrt(np.mean(delta * delta))),
        max_abs=float(np.max(np.abs(delta))),
        relative_l2=float(np.linalg.norm(delta) / reference_norm),
    )


def compare_centerlines(
    solver: LidDrivenCavity2D,
    reference: CavityCenterlineReference,
) -> CavityValidationResult:
    """Sample ``solver`` at the published coordinates and compare profiles.

    The numerical centerline is interpolated from the solver grid to the exact
    reference coordinates.  This is intentional: the published Ghia samples
    are not uniformly spaced and should not be replaced by nearby grid nodes.
    """

    if not np.isclose(solver.reynolds, reference.reynolds):
        raise ValueError("solver and reference Reynolds numbers must match")

    y, u, x, v = solver.centerline_profiles()
    u_at_reference = np.interp(reference.y_for_u, y, u)
    v_at_reference = np.interp(reference.x_for_v, x, v)
    diagnostics = solver.diagnostics()

    return CavityValidationResult(
        reynolds=float(solver.reynolds),
        resolution=int(solver.n),
        iterations=int(solver.iterations),
        final_update_inf=float(solver.final_update_inf),
        u_error=_profile_error(reference.u_centerline, u_at_reference),
        v_error=_profile_error(reference.v_centerline, v_at_reference),
        rms_divergence=float(diagnostics.rms_divergence),
    )


def run_ghia_re100_validation(
    *,
    resolution: int = 33,
    dt: float = 1.0e-3,
    max_iterations: int = 30000,
    tolerance: float = 1.0e-7,
) -> CavityValidationResult:
    """Solve the project's baseline ``Re=100`` cavity and compare with Ghia.

    This routine deliberately returns metrics rather than asserting that the
    result is "good".  The first task is to observe grid/refinement behavior;
    only then should a regression threshold be frozen into CI.
    """

    solver = LidDrivenCavity2D(n=resolution, reynolds=100.0, dt=dt)
    solver.solve(max_iterations=max_iterations, tolerance=tolerance)
    return compare_centerlines(solver, ghia_re100())
