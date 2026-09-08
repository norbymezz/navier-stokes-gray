"""Reference benchmark fields and diagnostics.

The functions here intentionally begin with exact/reference quantities rather
than with a custom CFD solver. They provide ground truth that later numerical
implementations must reproduce before project-specific topology is added.
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

    Convention:
        u = -U0 cos(kx) sin(ky) exp(-2 nu k^2 t)
        v =  U0 sin(kx) cos(ky) exp(-2 nu k^2 t)
    """

    decay = np.exp(-2.0 * nu * k * k * t)
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
    """Return a pressure field compatible with the 2-D Taylor–Green vortex."""

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
    """Return scalar vorticity omega_z = dv/dx - du/dy."""

    decay = np.exp(-2.0 * nu * k * k * t)
    return 2.0 * u0 * k * np.cos(k * x) * np.cos(k * y) * decay


def relative_l2(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Relative L2 error with a strict zero-reference guard."""

    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    denom = np.linalg.norm(reference.ravel())
    num = np.linalg.norm((candidate - reference).ravel())
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
    """Pointwise kinetic energy per unit mass: 1/2 |u|^2."""

    return 0.5 * (np.asarray(u) ** 2 + np.asarray(v) ** 2)


def scalar_variance(c: np.ndarray) -> float:
    """Spatial variance-like mixing diagnostic using the discrete mean."""

    c = np.asarray(c, dtype=float)
    return float(np.mean((c - np.mean(c)) ** 2))


def reversibility_error(initial: np.ndarray, final: np.ndarray) -> float:
    """Normalized error after a forward/reverse transport protocol."""

    return relative_l2(np.asarray(initial, dtype=float), np.asarray(final, dtype=float))
