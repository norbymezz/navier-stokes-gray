"""Steady lid-driven-cavity solver used for the first wall-bounded benchmark.

This implementation is intentionally transparent rather than optimized.  It
uses the same streamfunction-vorticity variables as the Ghia et al. benchmark,
which makes sign conventions and centerline comparisons easier to audit.

Equations (nondimensional)
--------------------------

    Laplacian(psi) = -omega

    omega_t + u omega_x + v omega_y = (1/Re) Laplacian(omega)

with

    u = psi_y
    v = -psi_x.

The square cavity spans ``[0,1] x [0,1]``.  All walls are no-slip and the top
wall moves in ``+x`` with unit speed.

Numerical method
----------------

- uniform Cartesian grid;
- central differences in the interior;
- Jacobi iteration for the streamfunction Poisson equation;
- explicit Euler pseudo-time marching for vorticity;
- second-order wall-vorticity formulas (Thom-style boundary treatment).

This is not intended to compete with Ghia's multigrid solver.  It is a small,
auditable baseline whose centerline profiles can be compared with published
reference values before more sophisticated numerics are introduced.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CavityDiagnostics:
    """Convergence and field diagnostics for one cavity solution."""

    reynolds: float
    iterations: int
    final_update_inf: float
    max_speed: float
    rms_divergence: float


@dataclass
class LidDrivenCavity2D:
    """Uniform-grid streamfunction-vorticity lid-driven cavity solver."""

    n: int = 33
    reynolds: float = 100.0
    dt: float = 1.0e-3
    poisson_iterations: int = 50
    poisson_tolerance: float = 1.0e-8

    def __post_init__(self) -> None:
        if self.n < 9:
            raise ValueError("n must be at least 9")
        if self.reynolds <= 0.0:
            raise ValueError("reynolds must be positive")
        if self.dt <= 0.0:
            raise ValueError("dt must be positive")
        if self.poisson_iterations < 1:
            raise ValueError("poisson_iterations must be positive")

        self.h = 1.0 / (self.n - 1)
        self.x = np.linspace(0.0, 1.0, self.n)
        self.y = np.linspace(0.0, 1.0, self.n)
        self.psi = np.zeros((self.n, self.n), dtype=float)
        self.omega = np.zeros_like(self.psi)
        self.u = np.zeros_like(self.psi)
        self.v = np.zeros_like(self.psi)
        self.iterations = 0
        self.final_update_inf = float("inf")

        self._apply_vorticity_boundary_conditions()
        self._update_velocity()

    def _apply_vorticity_boundary_conditions(self) -> None:
        """Apply second-order wall-vorticity formulas.

        With ``psi = 0`` on every solid wall and the convention
        ``omega = -Laplacian(psi)``, a stationary wall has approximately

            omega_wall = -2 psi_adjacent / h^2.

        At the moving lid, ``u = psi_y = 1`` adds ``-2 U/h``.
        Corner values are averaged from their neighboring wall estimates to
        avoid choosing one wall formula arbitrarily.
        """

        h2 = self.h * self.h

        # Bottom and top walls (varying x, fixed y).
        self.omega[:, 0] = -2.0 * self.psi[:, 1] / h2
        self.omega[:, -1] = -2.0 * self.psi[:, -2] / h2 - 2.0 / self.h

        # Left and right stationary walls.
        self.omega[0, :] = -2.0 * self.psi[1, :] / h2
        self.omega[-1, :] = -2.0 * self.psi[-2, :] / h2

        # Corner values are not benchmark observables and the continuum
        # velocity is discontinuous there because the moving lid meets a
        # stationary side wall.  Averaging keeps the discrete treatment
        # symmetric and avoids injecting a one-sided arbitrary corner choice.
        self.omega[0, 0] = 0.5 * (self.omega[1, 0] + self.omega[0, 1])
        self.omega[-1, 0] = 0.5 * (self.omega[-2, 0] + self.omega[-1, 1])
        self.omega[0, -1] = 0.5 * (self.omega[1, -1] + self.omega[0, -2])
        self.omega[-1, -1] = 0.5 * (self.omega[-2, -1] + self.omega[-1, -2])

    def _solve_streamfunction_poisson(self) -> None:
        """Approximately solve ``Laplacian(psi) = -omega`` by Jacobi sweeps."""

        h2 = self.h * self.h

        for _ in range(self.poisson_iterations):
            previous = self.psi.copy()
            self.psi[1:-1, 1:-1] = 0.25 * (
                previous[2:, 1:-1]
                + previous[:-2, 1:-1]
                + previous[1:-1, 2:]
                + previous[1:-1, :-2]
                + h2 * self.omega[1:-1, 1:-1]
            )

            # psi remains zero on the walls, representing an impermeable closed
            # cavity.  Stop the inner Poisson iteration once the update is tiny.
            if np.max(np.abs(self.psi - previous)) < self.poisson_tolerance:
                break

    def _update_velocity(self) -> None:
        """Recover velocity from the current streamfunction."""

        self.u.fill(0.0)
        self.v.fill(0.0)

        self.u[1:-1, 1:-1] = (
            self.psi[1:-1, 2:] - self.psi[1:-1, :-2]
        ) / (2.0 * self.h)
        self.v[1:-1, 1:-1] = -(
            self.psi[2:, 1:-1] - self.psi[:-2, 1:-1]
        ) / (2.0 * self.h)

        # No-slip on the stationary walls, unit +x speed on the lid.
        self.u[:, -1] = 1.0

    def _vorticity_rhs(self) -> np.ndarray:
        """Return the interior vorticity transport RHS."""

        rhs = np.zeros_like(self.omega)
        h = self.h
        h2 = h * h

        omega_x = (self.omega[2:, 1:-1] - self.omega[:-2, 1:-1]) / (2.0 * h)
        omega_y = (self.omega[1:-1, 2:] - self.omega[1:-1, :-2]) / (2.0 * h)
        lap_omega = (
            self.omega[2:, 1:-1]
            + self.omega[:-2, 1:-1]
            + self.omega[1:-1, 2:]
            + self.omega[1:-1, :-2]
            - 4.0 * self.omega[1:-1, 1:-1]
        ) / h2

        rhs[1:-1, 1:-1] = -(
            self.u[1:-1, 1:-1] * omega_x
            + self.v[1:-1, 1:-1] * omega_y
        ) + lap_omega / self.reynolds
        return rhs

    def step(self) -> float:
        """Advance one pseudo-time iteration and return ``||delta omega||_inf``."""

        self._solve_streamfunction_poisson()
        self._apply_vorticity_boundary_conditions()
        self._update_velocity()

        rhs = self._vorticity_rhs()
        previous = self.omega.copy()
        self.omega[1:-1, 1:-1] += self.dt * rhs[1:-1, 1:-1]
        self._apply_vorticity_boundary_conditions()

        update = float(np.max(np.abs(self.omega - previous)))
        self.iterations += 1
        self.final_update_inf = update
        return update

    def solve(self, *, max_iterations: int = 30000, tolerance: float = 1.0e-7) -> None:
        """March to a steady state according to the vorticity update norm."""

        if max_iterations < 1:
            raise ValueError("max_iterations must be positive")
        if tolerance <= 0.0:
            raise ValueError("tolerance must be positive")

        for _ in range(max_iterations):
            update = self.step()
            if self.iterations > 100 and update < tolerance:
                break

        # Finish with a tighter streamfunction reconstruction so exported
        # centerline velocities are not limited by the inner Poisson tolerance.
        old_iterations = self.poisson_iterations
        old_tolerance = self.poisson_tolerance
        self.poisson_iterations = max(500, old_iterations)
        self.poisson_tolerance = min(1.0e-10, old_tolerance)
        self._solve_streamfunction_poisson()
        self.poisson_iterations = old_iterations
        self.poisson_tolerance = old_tolerance
        self._update_velocity()

    def centerline_profiles(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return centerline coordinates and velocities for benchmark comparison.

        The grid is required to have odd ``n`` so ``x=0.5`` and ``y=0.5`` are
        represented exactly.  This avoids adding interpolation error before the
        first comparison with the published Ghia centerline samples.
        """

        if self.n % 2 == 0:
            raise ValueError("centerline extraction requires odd n")

        mid = self.n // 2
        return self.y.copy(), self.u[mid, :].copy(), self.x.copy(), self.v[:, mid].copy()

    def diagnostics(self) -> CavityDiagnostics:
        """Return compact convergence and incompressibility diagnostics."""

        speed = np.sqrt(self.u * self.u + self.v * self.v)

        # Divergence is evaluated only where centered differences are available.
        du_dx = (self.u[2:, 1:-1] - self.u[:-2, 1:-1]) / (2.0 * self.h)
        dv_dy = (self.v[1:-1, 2:] - self.v[1:-1, :-2]) / (2.0 * self.h)
        divergence = du_dx + dv_dy

        return CavityDiagnostics(
            reynolds=float(self.reynolds),
            iterations=int(self.iterations),
            final_update_inf=float(self.final_update_inf),
            max_speed=float(np.max(speed)),
            rms_divergence=float(np.sqrt(np.mean(divergence * divergence))),
        )
