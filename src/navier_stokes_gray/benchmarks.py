"""Reference benchmark fields and diagnostics.

The functions here intentionally begin with exact/reference quantities rather
than with a custom CFD solver. They provide ground truth that later numerical
implementations must reproduce before project-specific topology is added.

Primary reference used for the 2-D Taylor–Green verification field:
- Taylor, G. I. & Green, A. E. (1937), *Mechanism of the Production of Small
  Eddies from Large Ones*, Proc. Royal Society A, 158, 499–521.
  DOI: https://doi.org/10.1098/rspa.1937.0036

The exact 2-D periodic form used below is a standard verification specialization
of the Taylor–Green family. The project documentation records both the original
paper and modern benchmark formulations so that the numerical target is always
traceable.
"""

from __future__ import annotations

import numpy as np


def taylor_green_velocity(
    x: np.ndarray,
    y: np.ndarray,
    t: float,
    *,
    nu: float = 0.01,
    u0: float = 1.0,
    k: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the 2-D Taylor–Green exact velocity field.

    Convention used by this repository::

        u = -U0 cos(kx) sin(ky) exp(-2 nu k^2 t)
        v =  U0 sin(kx) cos(ky) exp(-2 nu k^2 t)

    Why this function exists:
        Before trusting any custom boundary identification, D2/D3 operator, or
        Gray-code experiment, the numerical solver must reproduce this known
        periodic solution with convergent error under grid/time-step refinement.
    """

    # Both velocity components decay with the same viscous factor. Keeping this
    # factor explicit makes it easy to test temporal convergence independently
    # from the spatial trigonometric structure.
    decay = np.exp(-2.0 * nu * k * k * t)

    # The signs and sine/cosine pairing are chosen so that div(u) = 0 exactly in
    # the continuum problem. This makes divergence error a direct numerical test.
    u = -u0 * np.cos(k * x) * np.sin(k * y) * decay
    v = u0 * np.sin(k * x) * np.cos(k * y) * decay
    return u, v


def taylor_green_pressure(
    x: np.ndarray,
    y: np.ndarray,
    t: float,
    *,
    nu: float = 0.01,
    rho: float = 1.0,
    u0: float = 1.0,
    k: float = 1.0,
    constant: float = 0.0,
) -> np.ndarray:
    """Return a pressure field compatible with the 2-D Taylor–Green vortex.

    Pressure is defined only up to an additive constant for incompressible flow,
    so ``constant`` is intentionally exposed instead of silently fixing a gauge.
    """

    # Pressure decays twice as fast in the exponent because it scales with U^2.
    decay = np.exp(-4.0 * nu * k * k * t)
    dynamic = -0.25 * rho * u0 * u0 * (
        np.cos(2.0 * k * x) + np.cos(2.0 * k * y)
    ) * decay
    return dynamic + constant


def taylor_green_vorticity(
    x: np.ndarray,
    y: np.ndarray,
    t: float,
    *,
    nu: float = 0.01,
    u0: float = 1.0,
    k: float = 1.0,
) -> np.ndarray:
    """Return scalar vorticity ``omega_z = dv/dx - du/dy``.

    Keeping an analytic vorticity field gives us a second, derivative-sensitive
    validation target in addition to velocity and pressure.
    """

    decay = np.exp(-2.0 * nu * k * k * t)
    return 2.0 * u0 * k * np.cos(k * x) * np.cos(k * y) * decay


def relative_l2(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Relative L2 error with a strict zero-reference guard."""

    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    denom = np.linalg.norm(reference.ravel())
    num = np.linalg.norm((candidate - reference).ravel())

    # Do not hide a failure behind division by zero. If the reference is exactly
    # zero, only an exactly zero candidate is considered zero relative error.
    if denom == 0.0:
        return 0.0 if num == 0.0 else float("inf")
    return float(num / denom)


def vector_relative_l2(
    reference: tuple[np.ndarray, np.ndarray],
    candidate: tuple[np.ndarray, np.ndarray],
) -> float:
    """Relative L2 error for a 2-component vector field."""

    ru, rv = reference
    cu, cv = candidate
    r = np.concatenate((np.ravel(ru), np.ravel(rv)))
    c = np.concatenate((np.ravel(cu), np.ravel(cv)))
    return relative_l2(r, c)


def kinetic_energy_density(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Pointwise kinetic energy per unit mass: ``1/2 |u|^2``."""

    return 0.5 * (np.asarray(u) ** 2 + np.asarray(v) ** 2)


def scalar_variance(c: np.ndarray) -> float:
    """Spatial variance-like mixing diagnostic using the discrete mean.

    This is deliberately simple. Later mixing studies may add entropy-like or
    multiscale diagnostics, but the first benchmark should remain transparent.
    """

    c = np.asarray(c, dtype=float)
    return float(np.mean((c - np.mean(c)) ** 2))


def reversibility_error(initial: np.ndarray, final: np.ndarray) -> float:
    """Normalized error after a forward/reverse transport protocol.

    At very low Reynolds number, with an exactly reversed boundary protocol and
    negligible diffusion/numerical error, this quantity should approach zero.
    Departures from zero are therefore a useful diagnostic, not something to be
    interpreted automatically as turbulence.
    """

    return relative_l2(np.asarray(initial, dtype=float), np.asarray(final, dtype=float))
