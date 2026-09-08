"""Run the quantitative Taylor--Green temporal-convergence study.

This script is intentionally small and reproducible.  It does not embed
published comparison numbers.  For every run, the numerical velocity is
compared with the closed-form Taylor--Green solution evaluated at the same
final time.

Usage
-----

    python scripts/run_taylor_green_convergence.py

The default sequence halves ``dt`` three times.  The resulting table can be
copied directly into a report or issue without hand-calculating errors or
observed orders.
"""

from __future__ import annotations

from navier_stokes_gray.convergence import observed_orders, temporal_convergence_study


def main() -> None:
    """Execute the standard temporal-refinement study and print a compact table."""

    # These are numerical experiment settings, not literature benchmark values.
    # The exact reference field is regenerated analytically by the package.
    dts = [0.08, 0.04, 0.02, 0.01]
    resolution = 32
    nu = 0.5
    final_time = 0.4

    points = temporal_convergence_study(
        dts,
        resolution=resolution,
        nu=nu,
        final_time=final_time,
    )
    rates = observed_orders(points)

    print("Taylor--Green temporal convergence")
    print(f"resolution={resolution}, nu={nu}, final_time={final_time}")
    print()
    print(f"{'dt':>10} {'relative L2 velocity error':>28} {'rms(div u)':>16} {'order':>10}")
    print("-" * 70)

    for index, point in enumerate(points):
        order_text = "-" if index == 0 else f"{rates[index - 1].order:.6f}"
        print(
            f"{point.dt:10.5g} "
            f"{point.velocity_relative_l2:28.12e} "
            f"{point.rms_divergence:16.8e} "
            f"{order_text:>10}"
        )


if __name__ == "__main__":
    main()
