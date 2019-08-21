from __future__ import absolute_import, division, print_function

import copy
import pytest

from scitbx import matrix

from . import find_unique_vectors, is_approximate_integer_multiple
from . import vector_group, group_vectors


def test_is_approximate_integer_multiple():
    v1 = matrix.col((34, 50, 20))
    v2 = 2.1 * v1.rotate_around_origin(matrix.col((1, 0, 0)), 3, deg=True)
    assert is_approximate_integer_multiple(v1, v2)
    assert not is_approximate_integer_multiple(v1, v2, relative_length_tolerance=0.01)
    assert not is_approximate_integer_multiple(v1, v2, angular_tolerance=2)


def test_find_unique_vectors():
    expected_unique_vectors = [
        matrix.col((34, 50, 20)),
        matrix.col((45, 67, 12)),
        matrix.col((12, 23, 34)),
    ]
    vectors = copy.deepcopy(expected_unique_vectors)
    vectors.append(
        2.1 * vectors[0].rotate_around_origin(matrix.col((1, 0, 0)), 2, deg=True)
    )
    vectors.append(
        3.9 * vectors[1].rotate_around_origin(matrix.col((0, 1, 0)), -1.2, deg=True)
    )
    vectors.append(
        3.1 * vectors[2].rotate_around_origin(matrix.col((1, 1, 1)), 0.5, deg=True)
    )
    sort_func = lambda v: v.length()
    vectors.sort(key=sort_func)
    unique_vectors = find_unique_vectors(vectors, max_vectors=30)
    assert len(unique_vectors) == len(expected_unique_vectors)
    assert sorted(unique_vectors, key=sort_func) == sorted(
        expected_unique_vectors, key=sort_func
    )


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
