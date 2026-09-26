"""Discrete structure helpers for the exploratory Gray/scale track."""

from __future__ import annotations

from collections.abc import Iterable


def power_sum(n: int, power: int) -> int:
    """Return sum(k**power for k=1..n) exactly."""
    if n < 0:
        raise ValueError("n must be >= 0")
    if power < 0:
        raise ValueError("power must be >= 0")
    return sum(k**power for k in range(1, n + 1))


def finite_difference(values: Iterable[int | float]) -> list[int | float]:
    """Return first forward differences."""
    seq = list(values)
    return [b - a for a, b in zip(seq, seq[1:])]


def binary_to_gray(value: int) -> int:
    """Convert a non-negative integer to reflected binary Gray code."""
    if value < 0:
        raise ValueError("value must be non-negative")
    return value ^ (value >> 1)


def gray_to_binary(gray: int) -> int:
    """Invert reflected binary Gray code by cumulative XOR."""
    if gray < 0:
        raise ValueError("gray must be non-negative")
    value = 0
    while gray:
        value ^= gray
        gray >>= 1
    return value


def two_adic_valuation(value: int) -> int:
    """Return nu_2(value), the exponent of 2 dividing a positive integer."""
    if value <= 0:
        raise ValueError("value must be positive")
    return (value & -value).bit_length() - 1


def ruler_value(transition_index: int) -> int:
    """Return the 1-based bit changed at this reflected-Gray transition."""
    return two_adic_valuation(transition_index) + 1


def gray_change_sequence(n_bits: int) -> list[int]:
    """Return 1-based changed-bit positions for a full n-bit Gray traversal."""
    if n_bits < 1:
        raise ValueError("n_bits must be >= 1")
    return [ruler_value(k) for k in range(1, 1 << n_bits)]


def relative_turn_path(
    lengths: Iterable[float],
    relative_turns_degrees: Iterable[float],
    *,
    initial_angle_degrees: float = 0.0,
) -> list[tuple[float, float]]:
    """Return vertices of a planar path with relative turns."""
    import math

    lengths_list = list(lengths)
    turns = list(relative_turns_degrees)

    if not lengths_list:
        return [(0.0, 0.0)]
    if len(turns) != len(lengths_list) - 1:
        raise ValueError("need exactly len(lengths)-1 relative turns")

    angle = math.radians(initial_angle_degrees)
    x = y = 0.0
    points = [(x, y)]

    for index, length in enumerate(lengths_list):
        if length < 0:
            raise ValueError("lengths must be non-negative")
        if index > 0:
            angle += math.radians(turns[index - 1])
        x += length * math.cos(angle)
        y += length * math.sin(angle)
        points.append((x, y))

    return points
