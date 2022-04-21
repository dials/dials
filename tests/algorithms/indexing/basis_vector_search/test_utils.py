from __future__ import annotations

import pytest

from scitbx import matrix

from dials.algorithms.indexing.basis_vector_search.utils import (
    group_vectors,
    is_approximate_integer_multiple,
    vector_group,
)


def test_is_approximate_integer_multiple():
    v1 = matrix.col((34, 50, 20))
    v2 = 2.1 * v1.rotate_around_origin(matrix.col((1, 0, 0)), 3, deg=True)
    assert is_approximate_integer_multiple(v1, v2)
    assert not is_approximate_integer_multiple(v1, v2, relative_length_tolerance=0.01)
    assert not is_approximate_integer_multiple(v1, v2, angular_tolerance=2)


def test_vector_group():
    group = vector_group()
    group.append(matrix.col((12, 13, 14)), weight=10)
    group.append(matrix.col((12.1, 13.1, 14.1)), weight=20)
    group.append(matrix.col((12.2, 13.2, 14.2)), weight=30)
    assert group.mean.elems == (12.1, 13.1, 14.1)
    assert group.weights == [10, 20, 30]


def test_group_vectors():
    vectors = [
        matrix.col(v)
        for v in (
            (12, 13, 14),
            (12.1, 13.1, 14.1),
            (12.2, 13.2, 14.2),
            (14, 15, 16),
            (14.1, 15.1, 16.1),
            (14.2, 15.2, 16.2),
        )
    ]
    weights = [1, 2, 3, 4, 5, 6]
    groups = group_vectors(vectors, weights)
    assert len(groups) == 2
    groups[0].mean.elems == pytest.approx((12.1, 13.1, 14.1))
    groups[0].weights == [1, 2, 3]
    groups[1].mean.elems == pytest.approx((14.1, 15.1, 16.1))
    groups[1].weights == [4, 5, 6]

    # test max_groups parameter
    groups = group_vectors(vectors, weights, max_groups=1)
    assert len(groups) == 1
    groups[0].mean.elems == pytest.approx((12.1, 13.1, 14.1))
    groups[0].weights == [1, 2, 3]
