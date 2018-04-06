from __future__ import absolute_import, division, print_function

def test_polygon():
  from dials.algorithms.polygon import polygon

  x = 1
  y = 1
  vertices = [(0,0), (2,0), (2,2), (0,2)]
  poly = polygon(vertices)
  assert poly.is_inside(x,y)

  poly = polygon([(3,5), (40,90), (80,70), (50,50), (70,20)])
  for p in [(42,40), (16,25), (64,67), (16,30), (48, 45), (21,30)]:
    assert poly.is_inside(p[0], p[1])
  for p in [(59,4), (70,15), (57,14), (21,78), (37,100), (88,89)]:
    assert not poly.is_inside(p[0], p[1])

  if 0:
    # for visual confirmation of algorithm
    from scitbx.array_family import flex
    inside_points = flex.vec2_double()
    outside_points = flex.vec2_double()
    import random
    x_max = 100
    y_max = 100
    for i in range(1000):
      x = random.randint(0,x_max)
      y = random.randint(0,y_max)
      is_inside = poly.is_inside(x, y)
      if is_inside:
        inside_points.append((x,y))
      else:
        outside_points.append((x,y))
    from matplotlib import pyplot
    from matplotlib.patches import Polygon
    v = poly.vertices + poly.vertices[:1]
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.add_patch(Polygon(poly.vertices, closed=True, fill=False))
    inside_x, inside_y = inside_points.parts()
    outside_x, outside_y = outside_points.parts()
    ax.scatter(inside_x, inside_y, marker='+', c='r')
    ax.scatter(outside_x, outside_y, marker='+', c='b')
    pyplot.show()
