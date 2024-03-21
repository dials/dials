from __future__ import annotations

from dials.util.system import CPU_COUNT


def test_available_cpus():
    # We can't make any assumptions about the testing environment,
    # but we know there will be at least one available core, and
    # the function must return a positive integer in any case.
    assert CPU_COUNT >= 1
