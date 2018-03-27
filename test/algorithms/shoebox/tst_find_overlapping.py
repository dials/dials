
from __future__ import absolute_import, division, print_function

from dials.algorithms.shoebox import find_overlapping

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_single_panel()
    self.tst_multiple_panels()

  def tst_single_panel(self):
    from dials.array_family import flex
    from random import randint

    nrefl = 1000

    # Generate bboxes
    bbox = flex.int6(nrefl)
    for i in range(nrefl):
      x0 = randint(0, 500)
      y0 = randint(0, 500)
      z0 = randint(0, 10)
      x1 = x0 + randint(2, 10)
      y1 = y0 + randint(2, 10)
      z1 = z0 + randint(2, 10)
      bbox[i] = (x0, x1, y0, y1, z0, z1)

    # Find the overlaps
    overlaps = find_overlapping(bbox)

    assert(overlaps.num_vertices() == nrefl)
    overlaps2 = self.brute_force(bbox)
    assert(overlaps.num_edges() == len(overlaps2))
    edges = {}
    for edge in overlaps2:
      edge = (min(edge), max(edge))
      edges[edge] = None
    for edge in overlaps.edges():
      edge = overlaps.source(edge), overlaps.target(edge)
      edge = (min(edge), max(edge))
      assert(edge in edges)

  def tst_multiple_panels(self):
    from dials.array_family import flex
    from random import randint

    nrefl = 1000

    # Generate bboxes
    bbox = flex.int6(nrefl)
    panel = flex.size_t(nrefl)
    for i in range(nrefl):
      x0 = randint(0, 500)
      y0 = randint(0, 500)
      z0 = randint(0, 10)
      x1 = x0 + randint(2, 10)
      y1 = y0 + randint(2, 10)
      z1 = z0 + randint(2, 10)
      bbox[i] = (x0, x1, y0, y1, z0, z1)
      panel[i] = randint(0,2)

    # Find the overlaps
    overlaps = find_overlapping(bbox, panel)
    assert(overlaps.num_vertices() == nrefl)
    overlaps2 = self.brute_force(bbox, panel)
    assert(overlaps.num_edges() == len(overlaps2))
    edges = {}
    for edge in overlaps2:
      edge = (min(edge), max(edge))
      edges[edge] = None
    for edge in overlaps.edges():
      edge = (overlaps.source(edge), overlaps.target(edge))
      edge = (min(edge), max(edge))
      assert(edge in edges)

  def brute_force(self, bbox, panel = None):
    overlaps = []
    if panel is None:
      for j in range(len(bbox)-1):
        jx0, jx1, jy0, jy1, jz0, jz1 = bbox[j]
        for i in range(j+1, len(bbox)):
          ix0, ix1, iy0, iy1, iz0, iz1 = bbox[i]
          if (not (ix0 >= jx1 or jx0 >= ix1 or
                   iy0 >= jy1 or jy0 >= iy1 or
                   iz0 >= jz1 or jz0 >= iz1)):
            overlaps.append((i, j))
    else:
      for j in range(len(bbox)-1):
        jx0, jx1, jy0, jy1, jz0, jz1 = bbox[j]
        for i in range(j+1, len(bbox)):
          ix0, ix1, iy0, iy1, iz0, iz1 = bbox[i]
          if panel[j] != panel[i]:
            continue
          if (not (ix0 >= jx1 or jx0 >= ix1 or
                   iy0 >= jy1 or jy0 >= iy1 or
                   iz0 >= jz1 or jz0 >= iz1)):
            overlaps.append((i, j))

    return overlaps

if __name__ == '__main__':

  test = Test()
  test.run()
