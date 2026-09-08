"""Quantitative convergence utilities for the Taylor--Green validation stage.

This module keeps convergence analysis separate from the solver so that the
same machinery can later be reused for other benchmarks.  No published table
values are embedded here: Taylor--Green references are regenerated from the
closed-form solution in :mod:`navier_stokes_gray.benchmarks`.

The first study deliberately varies the time step while keeping the spatial
problem exactly resolved.  The 2-D Taylor--Green field is a single Fourier
mode, so a modest periodic grid represents it to roundoff.  This isolates the
order of the RK4 time integrator instead of mixing temporal and spatial error.

A later spatial-convergence study should use a smooth multi-mode manufactured
solution or a benchmark whose spatial truncation error is non-zero; refining a
single exactly represented Fourier mode would not measure spectral spatial
convergence in a meaningful way.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log
from typing import Iterable

import numpy as np

from .benchmarks import (
    taylor_green_velocity,
    taylor_green_vorticity,
    vector_relative_l2,
)
from .solver2d import PeriodicVorticitySolver2D
from .spectral import make_grid


@dataclass(frozen=True)
class ConvergencePoint:
    """One run in a refinement study."""

    resolution: int
    dt: float
    final_time: float
    velocity_relative_l2: float
    rms_divergence: float


@dataclass(frozen=True)
class ConvergenceRate:
    """Observed order between two consecutive refinement points."""

    coarse_dt: float
    fine_dt: float
    coarse_error: float
    fine_error: float
    order: float


def _mesh(resolution: int) -> tuple[np.ndarray, np.ndarray]:
    """Return an ``ij`` mesh on the project's standard 2*pi periodic domain."""

    grid = make_grid(resolution)
    return np.meshgrid(grid.x, grid.y, indexing="ij")


def run_taylor_green_case(
    *,
    resolution: int = 32,
    nu: float = 0.05,
    final_time: float = 0.2,
    dt: float = 0.01,
) -> ConvergencePoint:
    """Run one Taylor--Green case and compare it with the exact solution.

    The exact velocity at ``final_time`` is evaluated directly from the
    analytical expression used by the benchmark module.  Therefore the error
    reported here measures the numerical evolution relative to ground truth,
    rather than relative to another numerical run.
    """

    if resolution < 4:
        raise ValueError("resolution must be at least 4")
    if nu < 0.0:
        raise ValueError("nu must be non-negative")
    if final_time <= 0.0:
        raise ValueError("final_time must be positive")
    if dt <= 0.0:
        raise ValueError("dt must be positive")

    grid = make_grid(resolution)
    x, y = np.meshgrid(grid.x, grid.y, indexing="ij")

    omega0 = taylor_green_vorticity(x, y, 0.0, nu=nu)
    solver = PeriodicVorticitySolver2D(grid=grid, nu=nu, omega=omega0)
    solver.run(final_time=final_time, dt=dt)

    numerical_velocity = solver.velocity()
    exact_velocity = taylor_green_velocity(x, y, final_time, nu=nu)
    error = vector_relative_l2(exact_velocity, numerical_velocity)
    diagnostics = solver.diagnostics()

    return ConvergencePoint(
        resolution=resolution,
        dt=float(dt),
        final_time=float(final_time),
        velocity_relative_l2=float(error),
        rms_divergence=float(diagnostics.rms_divergence),
    )


def temporal_convergence_study(
    dts: Iterable[float],
    *,
    resolution: int = 32,
    nu: float = 0.05,
    final_time: float = 0.2,
) -> list[ConvergencePoint]:
    """Run a Taylor--Green time-step refinement study.

    ``dts`` are sorted from coarse to fine.  Using one fixed, sufficiently
    resolved grid intentionally isolates time-integration error.
    """

    values = sorted((float(dt) for dt in dts), reverse=True)
    if len(values) < 2:
        raise ValueError("at least two time steps are required")
    if any(dt <= 0.0 for dt in values):
        raise ValueError("all time steps must be positive")

    return [
        run_taylor_green_case(
            resolution=resolution,
            nu=nu,
            final_time=final_time,
            dt=dt,
        )
        for dt in values
    ]


def observed_orders(points: Iterable[ConvergencePoint]) -> list[ConvergenceRate]:
    """Compute pairwise observed order ``p`` from ``error ~ C * dt**p``.

    For two consecutive refinements,

        p = log(E_coarse / E_fine) / log(dt_coarse / dt_fine).

    Zero or non-finite errors are rejected because they do not define a useful
    logarithmic convergence rate.  This can occur once roundoff dominates.
    """

    pts = list(points)
    rates: list[ConvergenceRate] = []

    for coarse, fine in zip(pts, pts[1:]):
        if not (coarse.dt > fine.dt > 0.0):
            raise ValueError("points must be ordered from coarse to fine dt")

        ec = coarse.velocity_relative_l2
        ef = fine.velocity_relative_l2
        if not (np.isfinite(ec) and np.isfinite(ef) and ec > 0.0 and ef > 0.0):
            raise ValueError("positive finite errors are required")

        order = log(ec / ef) / log(coarse.dt / fine.dt)
        rates.append(
            ConvergenceRate(
                coarse_dt=coarse.dt,
                fine_dt=fine.dt,
                coarse_error=ec,
                fine_error=ef,
                order=float(order),
            )
        )

    return rates
