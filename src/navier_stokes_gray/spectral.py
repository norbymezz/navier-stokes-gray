"""Small periodic 2-D spectral operators used by the validation solver.

The first numerical stage of the project is intentionally restricted to a
smooth doubly-periodic domain.  That removes solid-wall corner singularities
and gives us a clean setting in which Taylor--Green can be compared against an
exact solution before any project-specific topology is introduced.

Conventions
-----------
We use a streamfunction ``psi`` such that

    u = d psi / d y
    v = - d psi / d x

and scalar vorticity

    omega = d v / d x - d u / d y = - Laplacian(psi).

With this convention the Fourier-space inversion is

    psi_hat = omega_hat / |k|^2

for non-zero modes.  The zero mode of ``psi`` is a gauge and is set to zero.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SpectralGrid2D:
    """Uniform square periodic grid and its Fourier wave numbers."""

    n: int
    length: float
    x: np.ndarray
    y: np.ndarray
    kx: np.ndarray
    ky: np.ndarray
    k2: np.ndarray
    dealias_mask: np.ndarray


def make_grid(n: int, length: float = 2.0 * np.pi) -> SpectralGrid2D:
    """Build a square periodic Fourier grid with a 2/3 de-aliasing mask.

    Parameters
    ----------
    n:
        Number of grid points in each direction.  An even value is preferred
        because the Nyquist mode is then represented in the usual FFT layout.
    length:
        Physical side length of the square domain.

    Notes
    -----
    The de-aliasing mask removes modes above two thirds of the largest
    resolvable one-dimensional wave number.  This is the classical 2/3 rule
    used for pseudo-spectral evaluation of quadratic nonlinearities.
    """

    if n < 4:
        raise ValueError("n must be at least 4")
    if length <= 0.0:
        raise ValueError("length must be positive")

    x = np.linspace(0.0, length, n, endpoint=False)
    y = x.copy()
    dx = length / n
    k_1d = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    kx, ky = np.meshgrid(k_1d, k_1d, indexing="ij")
    k2 = kx * kx + ky * ky

    # 2/3 rule: keep only modes whose individual Cartesian components lie
    # inside the central two-thirds band.  We compute the cutoff from the
    # largest absolute FFT wave number rather than assuming a particular L.
    cutoff = (2.0 / 3.0) * np.max(np.abs(k_1d))
    dealias_mask = (np.abs(kx) <= cutoff) & (np.abs(ky) <= cutoff)

    return SpectralGrid2D(
        n=n,
        length=length,
        x=x,
        y=y,
        kx=kx,
        ky=ky,
        k2=k2,
        dealias_mask=dealias_mask,
    )


def derivative(field: np.ndarray, grid: SpectralGrid2D, axis: int) -> np.ndarray:
    """Return a periodic spectral derivative along x (axis=0) or y (axis=1)."""

    field = np.asarray(field, dtype=float)
    if field.shape != (grid.n, grid.n):
        raise ValueError("field shape must match grid")
    if axis not in (0, 1):
        raise ValueError("axis must be 0 (x) or 1 (y)")

    field_hat = np.fft.fftn(field)
    multiplier = 1j * (grid.kx if axis == 0 else grid.ky)
    return np.fft.ifftn(multiplier * field_hat).real


def laplacian(field: np.ndarray, grid: SpectralGrid2D) -> np.ndarray:
    """Return the periodic spectral Laplacian of a scalar field."""

    field = np.asarray(field, dtype=float)
    if field.shape != (grid.n, grid.n):
        raise ValueError("field shape must match grid")
    return np.fft.ifftn(-grid.k2 * np.fft.fftn(field)).real


def velocity_from_vorticity(
    omega: np.ndarray,
    grid: SpectralGrid2D,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover an incompressible velocity field from scalar vorticity.

    The inversion solves ``-Laplacian(psi) = omega`` in Fourier space, with
    the mean streamfunction fixed to zero.  Velocity is then differentiated
    spectrally from ``psi``.  The resulting field is divergence free up to
    floating-point roundoff for resolved periodic modes.
    """

    omega = np.asarray(omega, dtype=float)
    if omega.shape != (grid.n, grid.n):
        raise ValueError("omega shape must match grid")

    omega_hat = np.fft.fftn(omega)
    psi_hat = np.zeros_like(omega_hat, dtype=complex)
    nonzero = grid.k2 != 0.0
    psi_hat[nonzero] = omega_hat[nonzero] / grid.k2[nonzero]

    u_hat = 1j * grid.ky * psi_hat
    v_hat = -1j * grid.kx * psi_hat
    u = np.fft.ifftn(u_hat).real
    v = np.fft.ifftn(v_hat).real
    return u, v


def jacobian(
    streamfunction: np.ndarray,
    scalar: np.ndarray,
    grid: SpectralGrid2D,
    *,
    dealias: bool = True,
) -> np.ndarray:
    """Return J(psi, scalar) = psi_x scalar_y - psi_y scalar_x.

    For the velocity convention used in this module, ``u grad(scalar)`` is
    ``-J(psi, scalar)``.  The vorticity equation can therefore be written

        omega_t = J(psi, omega) + nu Laplacian(omega).

    The product is formed in physical space and optionally de-aliased after
    transforming back to Fourier space.
    """

    psi_x = derivative(streamfunction, grid, axis=0)
    psi_y = derivative(streamfunction, grid, axis=1)
    s_x = derivative(scalar, grid, axis=0)
    s_y = derivative(scalar, grid, axis=1)
    value = psi_x * s_y - psi_y * s_x

    if not dealias:
        return value

    value_hat = np.fft.fftn(value)
    value_hat *= grid.dealias_mask
    return np.fft.ifftn(value_hat).real


def streamfunction_from_vorticity(
    omega: np.ndarray,
    grid: SpectralGrid2D,
) -> np.ndarray:
    """Return zero-mean streamfunction satisfying -Laplacian(psi)=omega."""

    omega = np.asarray(omega, dtype=float)
    omega_hat = np.fft.fftn(omega)
    psi_hat = np.zeros_like(omega_hat, dtype=complex)
    nonzero = grid.k2 != 0.0
    psi_hat[nonzero] = omega_hat[nonzero] / grid.k2[nonzero]
    return np.fft.ifftn(psi_hat).real
