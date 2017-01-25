from __future__ import absolute_import, division
from scitbx.array_family import flex # import dependency
from dials_algorithms_polygon_ext import *

class polygon(object):
  def __init__(self, vertices):
    assert len(vertices) > 2
    self.vertices = vertices

  def is_inside(self, x, y):
    # http://en.wikipedia.org/wiki/Point_in_polygon
    # http://en.wikipedia.org/wiki/Even-odd_rule
    poly = self.vertices
    num = len(poly)
    i = 0
    j = num - 1
    inside = False
    for i in range(num):
      if  ((poly[i][1] > y) != (poly[j][1] > y)) and \
          (x < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / (poly[j][1] - poly[i][1]) + poly[i][0]):
        inside = not inside
      j = i
    return inside
