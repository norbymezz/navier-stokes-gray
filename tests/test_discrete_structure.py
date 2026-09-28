from __future__ import annotations

import math

import pytest

from navier_stokes_gray.discrete_structure import (
    binary_to_gray,
    finite_difference,
    gray_change_sequence,
    gray_to_binary,
    power_sum,
    relative_turn_path,
    ruler_value,
)


def test_power_sum_relations() -> None:
    for n in range(1, 20):
        s1 = power_sum(n, 1)
        s2 = power_sum(n, 2)
        s3 = power_sum(n, 3)
        assert 6 * s2 == n * (n + 1) * (2 * n + 1)
        assert s3 == s1 * s1
        assert 3 * s2 == s1 * (2 * n + 1)


def test_integer_scale_subsequence_is_n_congruent_one_mod_three() -> None:
    integer_cases = [n for n in range(1, 40) if (2 * n + 1) % 3 == 0]
    assert integer_cases == list(range(1, 40, 3))


def test_finite_difference_recovers_terms() -> None:
    cumulative = [power_sum(n, 3) for n in range(0, 9)]
    assert finite_difference(cumulative) == [n**3 for n in range(1, 9)]


def test_gray_round_trip() -> None:
    for value in range(256):
        assert gray_to_binary(binary_to_gray(value)) == value


def test_ruler_sequence_matches_reflected_gray_changes() -> None:
    assert gray_change_sequence(4) == [
        1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1
    ]
    for k in range(1, 16):
        a = binary_to_gray(k - 1)
        b = binary_to_gray(k)
        changed_mask = a ^ b
        assert changed_mask == 1 << (ruler_value(k) - 1)


def test_relative_turn_path_uses_relative_angles() -> None:
    points = relative_turn_path([2.0, 1.0], [90.0])
    assert points[1] == pytest.approx((2.0, 0.0))
    assert points[2] == pytest.approx((2.0, 1.0))

    diagonal = relative_turn_path([math.sqrt(2.0)], [], initial_angle_degrees=45.0)
    assert diagonal[-1] == pytest.approx((1.0, 1.0))
