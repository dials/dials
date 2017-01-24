from __future__ import absolute_import, division

from dials.algorithms.spatial_indexing import make_spatial_index

class Test(object):

  def __init__(self):
    num = 100
    self.vec2_double = self.generate_vec2_double(num)
    self.vec3_double = self.generate_vec3_double(num)

  def run(self):
    num = 10
    self.tst_vec2_double(num)
    self.tst_vec3_double(num)

  def generate_vec2_double(self, num):
    from scitbx.array_family import flex
    from random import uniform
    x0, x1, y0, y1 = 0, 100, 0, 100
    result = flex.vec2_double(num)
    for i in range(num):
      result[i] = (uniform(x0, x1), uniform(y0, y1))
    return result

  def generate_vec3_double(self, num):
    from scitbx.array_family import flex
    from random import uniform
    x0, x1, y0, y1, z0, z1 = 0, 100, 0, 100, 0, 100
    result = flex.vec3_double(num)
    for i in range(num):
      result[i] = (uniform(x0, x1), uniform(y0, y1), uniform(z0, z1))
    return result

  def tst_vec2_double(self, num):
    from random import randint
    index = make_spatial_index(self.vec2_double)
    x, y = zip(*self.vec2_double)
    x0, x1, y0, y1 = (int(min(x)), int(max(x)),
                      int(min(y)), int(max(y)))
    for i in range(num):
      xx0 = randint(x0, x1)
      xx1 = randint(x0, x1)
      yy0 = randint(y0, y1)
      yy1 = randint(y0, y1)
      rg = (min(xx0, xx1), max(xx0, xx1)+1,
            min(yy0, yy1), max(yy0, yy1)+1)
      idx = index.query_range(rg)
      for j in idx:
        assert(self.inside2d(rg, self.vec2_double[j]))
      for j in set(range(len(self.vec2_double))).difference(set(idx)):
        assert(not self.inside2d(rg, self.vec2_double[j]))


    print 'OK'

  def tst_vec3_double(self, num):
    from random import randint
    index = make_spatial_index(self.vec3_double)
    x, y, z = zip(*self.vec3_double)
    x0, x1, y0, y1, z0, z1 = (int(min(x)), int(max(x)),
                              int(min(y)), int(max(y)),
                              int(min(z)), int(max(z)))
    for i in range(num):
      xx0 = randint(x0, x1)
      xx1 = randint(x0, x1)
      yy0 = randint(y0, y1)
      yy1 = randint(y0, y1)
      zz0 = randint(y0, y1)
      zz1 = randint(y0, y1)
      rg = (min(xx0, xx1), max(xx0, xx1)+1,
            min(yy0, yy1), max(yy0, yy1)+1,
            min(zz0, zz1), max(zz0, zz1)+1)
      idx = index.query_range(rg)
      for j in idx:
        assert(self.inside3d(rg, self.vec3_double[j]))
      for j in set(range(len(self.vec3_double))).difference(set(idx)):
        assert(not self.inside3d(rg, self.vec3_double[j]))

    print 'OK'

  def inside2d(self, rg, a):
    inside = (rg[0] <= a[0] and rg[1] >= a[0] and
              rg[2] <= a[1] and rg[3] >= a[1])
    return inside

  def inside3d(self, rg, a):
    inside = (rg[0] <= a[0] and rg[1] >= a[0] and
              rg[2] <= a[1] and rg[3] >= a[1] and
              rg[4] <= a[2] and rg[5] >= a[2])
    return inside

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
