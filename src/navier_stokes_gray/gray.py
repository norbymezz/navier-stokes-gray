"""Gray-code utilities for one-change-at-a-time simulation schedules."""

from __future__ import annotations


def gray_code(n_bits: int) -> list[int]:
    """Return the reflected binary Gray sequence as integers.

    Consecutive entries differ by exactly one bit.
    """
    if n_bits < 1:
        raise ValueError("n_bits must be >= 1")
    return [i ^ (i >> 1) for i in range(1 << n_bits)]


def hamming_distance(a: int, b: int) -> int:
    """Return the Hamming distance between two non-negative integers."""
    if a < 0 or b < 0:
        raise ValueError("Gray-state integers must be non-negative")
    return (a ^ b).bit_count()


def changed_bit(a: int, b: int) -> int:
    """Return the index of the unique changed bit between adjacent Gray states."""
    diff = a ^ b
    if diff == 0 or diff & (diff - 1):
        raise ValueError("states must differ by exactly one bit")
    return diff.bit_length() - 1


def gray_schedule(n_bits: int) -> list[dict[str, int | None]]:
    """Return a compact schedule with the changed bit at each transition."""
    states = gray_code(n_bits)
    schedule: list[dict[str, int | None]] = []
    previous: int | None = None
    for state in states:
        schedule.append(
            {
                "state": state,
                "changed_bit": None if previous is None else changed_bit(previous, state),
            }
        )
        previous = state
    return schedule
