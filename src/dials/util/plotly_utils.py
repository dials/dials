from __future__ import annotations

from typing import Sequence


def round_for_json(seq: Sequence, ndigits: int = 3) -> list:
    """
    Round floating-point numbers in a sequence to a specified number of digits
    for JSON serialization (https://github.com/dials/dials/issues/3190)

    Parameters
    ----------
    seq : Sequence
        The sequence of numbers to round.
    ndigits : int, optional
        The number of digits to round to (default is 3).

    Returns
    -------
    list
        The rounded sequence.
    """
    return [round(x, ndigits) if isinstance(x, float) else x for x in seq]
