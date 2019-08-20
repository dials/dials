from __future__ import absolute_import, division, print_function

import copy

from scitbx import matrix

from . import find_unique_vectors, is_approximate_integer_multiple


def test_is_approximate_integer_multiple():
    v1 = matrix.col((34, 50, 20))
    v2 = 2.1 * v1.rotate(matrix.col((1, 0, 0)), 3, deg=True)
    assert is_approximate_integer_multiple(v1, v2)
    assert not is_approximate_integer_multiple(v1, v2, relative_tolerance=0.01)
    assert not is_approximate_integer_multiple(v1, v2, angular_tolerance=2)


def test_find_unique_vectors():
    expected_unique_vectors = [
        matrix.col((34, 50, 20)),
        matrix.col((45, 67, 12)),
        matrix.col((12, 23, 34)),
    ]
    vectors = copy.deepcopy(expected_unique_vectors)
    vectors.append(2.1 * vectors[0].rotate(matrix.col((1, 0, 0)), 2, deg=True))
    vectors.append(3.9 * vectors[1].rotate(matrix.col((0, 1, 0)), -1.2, deg=True))
    vectors.append(3.1 * vectors[2].rotate(matrix.col((1, 1, 1)), 0.5, deg=True))
    sort_func = lambda v: v.length()
    vectors.sort(key=sort_func)
    unique_vectors = find_unique_vectors(vectors, max_vectors=30)
    assert len(unique_vectors) == len(expected_unique_vectors)
    assert sorted(unique_vectors, key=sort_func) == sorted(
        expected_unique_vectors, key=sort_func
    )
