"""Run the first wall-bounded benchmark: lid-driven cavity at Re=100.

The numerical solution is compared with the provenance-tracked centerline
profiles in :mod:`navier_stokes_gray.cavity_reference`.  The script reports
errors at the published sample coordinates rather than relying on a visual
plot.
"""

from __future__ import annotations

import numpy as np

from navier_stokes_gray.cavity import LidDrivenCavity2D
from navier_stokes_gray.cavity_reference import ghia_re100


def _rmse(reference: np.ndarray, candidate: np.ndarray) -> float:
    return float(np.sqrt(np.mean((candidate - reference) ** 2)))


def main() -> None:
    reference = ghia_re100()

    solver = LidDrivenCavity2D(
        n=33,
        reynolds=reference.reynolds,
        dt=1.0e-3,
        poisson_iterations=50,
        poisson_tolerance=1.0e-8,
    )
    solver.solve(max_iterations=30000, tolerance=1.0e-7)

    y, u_center, x, v_center = solver.centerline_profiles()
    u_at_reference = np.interp(reference.y_for_u, y, u_center)
    v_at_reference = np.interp(reference.x_for_v, x, v_center)

    u_rmse = _rmse(reference.u_centerline, u_at_reference)
    v_rmse = _rmse(reference.v_centerline, v_at_reference)
    u_max = float(np.max(np.abs(u_at_reference - reference.u_centerline)))
    v_max = float(np.max(np.abs(v_at_reference - reference.v_centerline)))

    diagnostics = solver.diagnostics()

    print("Lid-driven cavity benchmark: Re=100")
    print("reference: Ghia et al. (1982), centerline profiles")
    print(f"grid={solver.n}x{solver.n}, dt={solver.dt:g}")
    print(f"iterations={diagnostics.iterations}")
    print(f"final ||delta omega||_inf={diagnostics.final_update_inf:.8e}")
    print(f"rms(div u)={diagnostics.rms_divergence:.8e}")
    print()
    print(f"u(x=0.5) RMSE at 17 reference points = {u_rmse:.8e}")
    print(f"u(x=0.5) max abs error              = {u_max:.8e}")
    print(f"v(y=0.5) RMSE at 17 reference points = {v_rmse:.8e}")
    print(f"v(y=0.5) max abs error              = {v_max:.8e}")


if __name__ == "__main__":
    main()
