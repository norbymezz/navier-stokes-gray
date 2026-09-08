"""Minimal 2-D incompressible periodic Navier--Stokes solver.

This module is deliberately narrow: it evolves scalar vorticity on a square
periodic domain using a pseudo-spectral spatial discretization and classical
RK4 time stepping.  It exists to validate the numerical core against standard
reference problems before D2/D3, Gray ordering, twisted identifications or
ionic transport are introduced.

Governing equation
------------------

With the streamfunction convention used in :mod:`navier_stokes_gray.spectral`,

    u = d psi / d y
    v = - d psi / d x
    omega = - Laplacian(psi)

the 2-D incompressible vorticity equation is

    d omega / dt = J(psi, omega) + nu Laplacian(omega) + f_omega,

where J(a,b) = a_x b_y - a_y b_x.

The nonlinear Jacobian is evaluated pseudo-spectrally and de-aliased with the
2/3 rule.  No pressure solve is required in this vorticity formulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .spectral import (
    SpectralGrid2D,
    jacobian,
    laplacian,
    streamfunction_from_vorticity,
    velocity_from_vorticity,
)


VorticityForcing = Callable[[float, np.ndarray, SpectralGrid2D], np.ndarray]


@dataclass(frozen=True)
class SolverDiagnostics:
    """Compact diagnostics returned with a solver state."""

    time: float
    kinetic_energy: float
    enstrophy: float
    rms_divergence: float
    circulation: float


@dataclass
class PeriodicVorticitySolver2D:
    """Pseudo-spectral solver for 2-D periodic incompressible flow."""

    grid: SpectralGrid2D
    nu: float
    omega: np.ndarray
    time: float = 0.0
    forcing: VorticityForcing | None = None

    def __post_init__(self) -> None:
        if self.nu < 0.0:
            raise ValueError("nu must be non-negative")

        self.omega = np.asarray(self.omega, dtype=float).copy()
        if self.omega.shape != (self.grid.n, self.grid.n):
            raise ValueError("omega shape must match grid")

    def rhs(self, t: float, omega: np.ndarray) -> np.ndarray:
        """Evaluate the vorticity RHS at one RK stage."""

        psi = streamfunction_from_vorticity(omega, self.grid)

        # With u = psi_y and v = -psi_x, advection u·grad(omega) equals
        # -J(psi, omega).  Moving advection to the RHS therefore gives +J.
        rhs_value = jacobian(psi, omega, self.grid, dealias=True)
        rhs_value += self.nu * laplacian(omega, self.grid)

        if self.forcing is not None:
            forcing_value = np.asarray(self.forcing(t, omega, self.grid), dtype=float)
            if forcing_value.shape != omega.shape:
                raise ValueError("forcing must return an array matching omega")
            rhs_value += forcing_value

        return rhs_value

    def step(self, dt: float) -> None:
        """Advance one time step with classical fourth-order Runge--Kutta."""

        if dt <= 0.0:
            raise ValueError("dt must be positive")

        t0 = self.time
        w0 = self.omega

        k1 = self.rhs(t0, w0)
        k2 = self.rhs(t0 + 0.5 * dt, w0 + 0.5 * dt * k1)
        k3 = self.rhs(t0 + 0.5 * dt, w0 + 0.5 * dt * k2)
        k4 = self.rhs(t0 + dt, w0 + dt * k3)

        self.omega = w0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        self.time = t0 + dt

    def run(self, final_time: float, dt: float) -> None:
        """Advance until ``final_time`` without overshooting it."""

        if final_time < self.time:
            raise ValueError("final_time cannot be earlier than current time")
        if dt <= 0.0:
            raise ValueError("dt must be positive")

        while self.time < final_time:
            remaining = final_time - self.time
            self.step(min(dt, remaining))

    def velocity(self) -> tuple[np.ndarray, np.ndarray]:
        """Return velocity reconstructed from the current vorticity."""

        return velocity_from_vorticity(self.omega, self.grid)

    def diagnostics(self) -> SolverDiagnostics:
        """Compute integral and incompressibility diagnostics.

        The averages are discrete periodic-domain averages.  Multiplication by
        the physical domain area would convert the energy/enstrophy densities
        into total integrals; for regression tests the normalized averages are
        more convenient and independent of grid resolution.
        """

        u, v = self.velocity()
        speed2 = u * u + v * v
        kinetic_energy = 0.5 * float(np.mean(speed2))
        enstrophy = 0.5 * float(np.mean(self.omega * self.omega))

        # Spectral derivatives are used here rather than finite differences so
        # this diagnostic measures the solver representation itself.
        from .spectral import derivative

        divergence = derivative(u, self.grid, axis=0) + derivative(v, self.grid, axis=1)
        rms_divergence = float(np.sqrt(np.mean(divergence * divergence)))

        # In a periodic domain, the mean vorticity equals circulation density.
        circulation = float(np.mean(self.omega))

        return SolverDiagnostics(
            time=self.time,
            kinetic_energy=kinetic_energy,
            enstrophy=enstrophy,
            rms_divergence=rms_divergence,
            circulation=circulation,
        )
