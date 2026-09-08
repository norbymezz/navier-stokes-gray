"""Print the periodic spectral-operator spatial-accuracy study."""

from __future__ import annotations

from navier_stokes_gray.spatial_accuracy import spatial_accuracy_study


def main() -> None:
    """Run a coarse-to-fine analytic derivative/Laplacian comparison."""

    resolutions = [8, 12, 16, 24, 32]
    points = spatial_accuracy_study(resolutions)

    print("Periodic spectral spatial accuracy")
    print("analytic field: exp(sin(x) + 0.5*cos(2y))")
    print()
    print(
        f"{'N':>6} {'dx':>12} {'rel L2 fx':>16} "
        f"{'rel L2 fy':>16} {'rel L2 lap':>16}"
    )
    print("-" * 72)

    for point in points:
        print(
            f"{point.resolution:6d} "
            f"{point.dx:12.5e} "
            f"{point.derivative_x_relative_l2:16.8e} "
            f"{point.derivative_y_relative_l2:16.8e} "
            f"{point.laplacian_relative_l2:16.8e}"
        )


if __name__ == "__main__":
    main()
