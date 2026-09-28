"""Boundary-map primitives used by the first numerical experiments."""

from __future__ import annotations

import numpy as np


def periodic_scalar_top_to_bottom(field: np.ndarray) -> np.ndarray:
    """Map the top edge of a 2-D scalar field to the bottom edge periodically."""
    if field.ndim != 2:
        raise ValueError("field must be a 2-D scalar array")
    return field[-1, :].copy()


def twisted_scalar_top_to_bottom(field: np.ndarray) -> np.ndarray:
    """Orientation-reversing scalar map: x -> Lx - x."""
    if field.ndim != 2:
        raise ValueError("field must be a 2-D scalar array")
    return field[-1, ::-1].copy()


def twisted_vector_top_to_bottom(
    u_x: np.ndarray,
    u_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply a reflected top-to-bottom boundary map to a 2-D vector field.

    Reflection reverses the x component while preserving the y component.
    This function only implements the algebraic boundary map; a CFD solver
    must still enforce discrete flux consistency and incompressibility.
    """
    if u_x.shape != u_y.shape or u_x.ndim != 2:
        raise ValueError("u_x and u_y must be 2-D arrays with identical shape")

    mapped_x = -u_x[-1, ::-1].copy()
    mapped_y = u_y[-1, ::-1].copy()
    return mapped_x, mapped_y
