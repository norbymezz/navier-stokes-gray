"""Published lid-driven-cavity centerline reference data.

The first cavity comparison uses only the Reynolds-number-100 centerline
profiles attributed to Ghia, Ghia & Shin (1982).  The values below were not
entered from memory or conversation context.

Primary paper
-------------
U. Ghia, K. N. Ghia, C. T. Shin (1982),
"High-Re solutions for incompressible flow using the Navier-Stokes equations
and a multigrid method", Journal of Computational Physics 48, 387--411.
DOI: https://doi.org/10.1016/0021-9991(82)90058-4

Independent public transcription used for data entry/cross-check
---------------------------------------------------------------
Sandia SIERRA/Fuego Verification Manual, lid-driven cavity benchmark,
Table 22.1-1 (Re=100):
https://www.sandia.gov/app/uploads/sites/315/2025/11/Fuego_Verification_5_18.pdf

The Sandia table reproduces both centerline profiles and explicitly states the
unit-square configuration with lid speed 1 and viscosity 0.01, hence Re=100.
The repository keeps this provenance beside the numbers so future edits can be
audited against the source instead of relying on copied values with unknown
coordinate/sign conventions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CavityCenterlineReference:
    """Centerline samples for one nondimensional lid-driven-cavity case."""

    reynolds: float
    y_for_u: np.ndarray
    u_centerline: np.ndarray
    x_for_v: np.ndarray
    v_centerline: np.ndarray


def ghia_re100() -> CavityCenterlineReference:
    """Return Ghia-style centerline reference values for ``Re = 100``.

    Coordinate convention used by this project:

    - square cavity: ``0 <= x <= 1``, ``0 <= y <= 1``;
    - lid at ``y = 1`` moves in ``+x`` with nondimensional speed ``U = 1``;
    - ``u_centerline`` is sampled on ``x = 0.5`` as a function of ``y``;
    - ``v_centerline`` is sampled on ``y = 0.5`` as a function of ``x``.
    """

    # Table 22.1-1 in the Sandia verification manual reproduces the Ghia
    # Re=100 benchmark values.  The order here is ascending coordinate to make
    # interpolation with numpy straightforward; many publications print the
    # u-profile from y=1 downward instead.
    y_for_u = np.array(
        [
            0.0000,
            0.0547,
            0.0625,
            0.0703,
            0.1016,
            0.1719,
            0.2813,
            0.4531,
            0.5000,
            0.6172,
            0.7344,
            0.8516,
            0.9531,
            0.9609,
            0.9688,
            0.9766,
            1.0000,
        ],
        dtype=float,
    )
    u_centerline = np.array(
        [
            0.00000,
            -0.03717,
            -0.04192,
            -0.04775,
            -0.06434,
            -0.10150,
            -0.15662,
            -0.21090,
            -0.20581,
            -0.13641,
            0.00332,
            0.23151,
            0.68717,
            0.73722,
            0.78871,
            0.84123,
            1.00000,
        ],
        dtype=float,
    )

    x_for_v = np.array(
        [
            0.0000,
            0.0625,
            0.0703,
            0.0781,
            0.0938,
            0.1563,
            0.2266,
            0.2344,
            0.5000,
            0.8047,
            0.8594,
            0.9063,
            0.9453,
            0.9531,
            0.9609,
            0.9688,
            1.0000,
        ],
        dtype=float,
    )
    v_centerline = np.array(
        [
            0.00000,
            0.09233,
            0.10091,
            0.10890,
            0.12317,
            0.16077,
            0.17507,
            0.17527,
            0.05454,
            -0.24533,
            -0.22445,
            -0.16914,
            -0.10313,
            -0.08864,
            -0.07391,
            -0.05906,
            0.00000,
        ],
        dtype=float,
    )

    return CavityCenterlineReference(
        reynolds=100.0,
        y_for_u=y_for_u,
        u_centerline=u_centerline,
        x_for_v=x_for_v,
        v_centerline=v_centerline,
    )
