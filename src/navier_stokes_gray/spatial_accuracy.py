"""Spatial-accuracy checks for the periodic spectral operators.

Taylor--Green is a poor spatial-convergence test once its single Fourier mode
is resolved: further grid refinement changes almost nothing except roundoff.
This module therefore uses a smooth, periodic, non-band-limited analytic field
whose Fourier coefficients decay rapidly but do not terminate.

The field is deliberately elementary so its derivatives are available in
closed form and no external table or remembered number is required:

    f(x,y) = exp(sin(x) + 0.5 cos(2y)).

Because ``exp(sin(x))`` and ``exp(cos(2y))`` contain infinitely many Fourier
modes, increasing the grid resolution genuinely reduces spectral truncation
error until floating-point roundoff dominates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .benchmarks import relative_l2
from .spectral import derivative, laplacian, make_grid


@dataclass(frozen=True)
class SpatialAccuracyPoint:
    """Errors for one periodic-grid resolution."""

    resolution: int
    dx: float
    derivative_x_relative_l2: float
    derivative_y_relative_l2: float
    laplacian_relative_l2: float


def analytic_field_and_derivatives(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``f``, ``f_x``, ``f_y`` and ``Laplacian(f)`` analytically.

    The chosen function is smooth and exactly periodic on ``[0, 2*pi)^2``.
    It is not band-limited, which makes it suitable for measuring spectral
    spatial truncation error under grid refinement.
    """

    f = np.exp(np.sin(x) + 0.5 * np.cos(2.0 * y))
    fx = np.cos(x) * f
    fy = -np.sin(2.0 * y) * f

    # If f = exp(g), then f_xx = (g_xx + g_x^2) f and similarly in y.
    # Here g = sin(x) + 0.5 cos(2y).
    lap = (
        -np.sin(x)
        + np.cos(x) ** 2
        - 2.0 * np.cos(2.0 * y)
        + np.sin(2.0 * y) ** 2
    ) * f
    return f, fx, fy, lap


def run_spatial_accuracy_case(resolution: int) -> SpatialAccuracyPoint:
    """Compare spectral derivatives with the analytic periodic field."""

    if resolution < 4:
        raise ValueError("resolution must be at least 4")

    grid = make_grid(resolution)
    x, y = np.meshgrid(grid.x, grid.y, indexing="ij")
    field, fx_exact, fy_exact, lap_exact = analytic_field_and_derivatives(x, y)

    fx_num = derivative(field, grid, axis=0)
    fy_num = derivative(field, grid, axis=1)
    lap_num = laplacian(field, grid)

    return SpatialAccuracyPoint(
        resolution=resolution,
        dx=grid.length / grid.n,
        derivative_x_relative_l2=relative_l2(fx_exact, fx_num),
        derivative_y_relative_l2=relative_l2(fy_exact, fy_num),
        laplacian_relative_l2=relative_l2(lap_exact, lap_num),
    )


def spatial_accuracy_study(resolutions: Iterable[int]) -> list[SpatialAccuracyPoint]:
    """Run the analytic spatial-accuracy check from coarse to fine grids."""

    values = sorted({int(n) for n in resolutions})
    if len(values) < 2:
        raise ValueError("at least two distinct resolutions are required")
    return [run_spatial_accuracy_case(n) for n in values]
