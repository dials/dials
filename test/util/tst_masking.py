from __future__ import division

def exercise_polygon():
  from dials.util import is_inside_polygon
  from scitbx.array_family import flex

  poly = flex.vec2_double(((0,0), (1,0), (1,1), (0,1)))
  assert is_inside_polygon(poly, 0, 0)
  assert is_inside_polygon(poly, 0.5, 0.5)
  assert not is_inside_polygon(poly, 1, 1.01)

  points = flex.vec2_double(((0.3, 0.8), (0.3, 1.5), (-8,9), (0.00001, 0.9999)))
  assert list(is_inside_polygon(poly, points)) == [True, False, False, True]

def run():
  exercise_polygon()
  print 'OK'

if __name__ == '__main__':
  run()
