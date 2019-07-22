from __future__ import absolute_import, division, print_function

import math


def test_manhattan():
    from dials.algorithms.image.filter import manhattan_distance
    from scitbx.array_family import flex

    data = flex.bool(flex.grid(100, 100), True)

    for j in range(100):
        for i in range(100):
            if math.sqrt((j - 50) ** 2 + (i - 50) ** 2) <= 10.5:
                data[j, i] = False

    distance = manhattan_distance(data, True)

    M = distance[1:-1, 1:-1]
    D = data[1:-1, 1:-1]

    selection = D.as_1d().select(M.as_1d() == 0)
    assert selection.all_eq(True)

    while True:
        N = data[0:-2, 1:-1]
        S = data[2:, 1:-1]
        E = data[1:-1, 0:-2]
        W = data[1:-1, 2:]
        selection = M.as_1d() == 1
        neighbours = (
            N.select(selection)
            | S.select(selection)
            | E.select(selection)
            | W.select(selection)
        )
        assert neighbours.all_eq(True)
        indices = flex.size_t(range(len(D))).select(selection)
        for i in indices:
            D[i] = True
        data[1:-1, 1:-1] = D
        M -= 1
        if (M > 0).count(True) == 0:
            break


def test_chebyshev():
    from dials.algorithms.image.filter import chebyshev_distance
    from scitbx.array_family import flex

    data = flex.bool(flex.grid(100, 100), True)

    for j in range(100):
        for i in range(100):
            if math.sqrt((j - 50) ** 2 + (i - 50) ** 2) <= 10.5:
                data[j, i] = False

    distance = chebyshev_distance(data, True)

    M = distance[1:-1, 1:-1]
    D = data[1:-1, 1:-1]
    selection = D.as_1d().select(M.as_1d() == 0)
    assert selection.all_eq(True)

    while True:
        N = data[0:-2, 1:-1]
        S = data[2:, 1:-1]
        E = data[1:-1, 0:-2]
        W = data[1:-1, 2:]
        NE = data[0:-2, 0:-2]
        NW = data[0:-2, 2:]
        SE = data[2:, 0:-2]
        SW = data[2:, 2:]
        selection = M.as_1d() == 1
        neighbours = (
            N.select(selection)
            | S.select(selection)
            | E.select(selection)
            | W.select(selection)
            | NE.select(selection)
            | NW.select(selection)
            | SE.select(selection)
            | SW.select(selection)
        )
        assert neighbours.all_eq(True)
        indices = flex.size_t(range(len(D))).select(selection)
        for i in indices:
            D[i] = True
        data[1:-1, 1:-1] = D
        M -= 1
        if (M > 0).count(True) == 0:
            break
