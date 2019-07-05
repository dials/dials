from __future__ import absolute_import, division, print_function


def test_polygon():
    from dials.algorithms.polygon import polygon

    x = 1
    y = 1
    vertices = [(0, 0), (2, 0), (2, 2), (0, 2)]
    poly = polygon(vertices)
    assert poly.is_inside(x, y)

    poly = polygon([(3, 5), (40, 90), (80, 70), (50, 50), (70, 20)])
    for p in [(42, 40), (16, 25), (64, 67), (16, 30), (48, 45), (21, 30)]:
        assert poly.is_inside(p[0], p[1])
    for p in [(59, 4), (70, 15), (57, 14), (21, 78), (37, 100), (88, 89)]:
        assert not poly.is_inside(p[0], p[1])
