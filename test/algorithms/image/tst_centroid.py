from __future__ import absolute_import, division


class CentroidTest(object):

  EPS = 1e-7

  def __init__(self):
    self.generate_data()
    self.calculate_gold()

  def run(self):
    self.tst_centroid_points()
    self.tst_centroid_image()
    self.tst_centroid_masked_image()
    self.tst_centroid_bias()

  def tst_centroid_points(self):
    self.tst_centroid_points2d()
    self.tst_centroid_points3d()

  def tst_centroid_points2d(self):
    from dials.algorithms.image.centroid import centroid_points
    from scitbx import matrix

    centroid = centroid_points(self.pixels2d.as_1d(), self.points2d.as_1d())

    assert(abs(self.gold2d - matrix.col(centroid.mean())) < self.EPS)
    assert(abs(self.gold2dvar - matrix.col(centroid.variance())) < self.EPS)
    assert(abs(self.gold2dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS)
    print 'OK'

  def tst_centroid_points3d(self):
    from dials.algorithms.image.centroid import centroid_points
    from scitbx import matrix

    centroid = centroid_points(self.pixels3d.as_1d(), self.points3d.as_1d())

    assert(abs(self.gold3d - matrix.col(centroid.mean())) < self.EPS)
    assert(abs(self.gold3dvar - matrix.col(centroid.variance())) < self.EPS)
    assert(abs(self.gold3dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS)
    print 'OK'

  def tst_centroid_image(self):
    self.tst_centroid_image2d()
    self.tst_centroid_image3d()

  def tst_centroid_image2d(self):
    from dials.algorithms.image.centroid import centroid_image
    from scitbx import matrix

    centroid = centroid_image(self.pixels2d)

    assert(abs(self.gold2d - matrix.col(centroid.mean())) < self.EPS)
    assert(abs(self.gold2dvar - matrix.col(centroid.variance())) < self.EPS)
    assert(abs(self.gold2dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS)
    print 'OK'

  def tst_centroid_image3d(self):
    from dials.algorithms.image.centroid import centroid_image
    from scitbx import matrix

    centroid = centroid_image(self.pixels3d)

    assert(abs(self.gold3d - matrix.col(centroid.mean())) < self.EPS)
    assert(abs(self.gold3dvar - matrix.col(centroid.variance())) < self.EPS)
    assert(abs(self.gold3dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS)
    print 'OK'

  def tst_centroid_masked_image(self):
    self.tst_centroid_masked_image2d()
    self.tst_centroid_masked_image3d()

  def tst_centroid_masked_image2d(self):
    from dials.algorithms.image.centroid import centroid_image
    from scitbx import matrix

    centroid = centroid_image(self.pixels2d, self.mask2d)

    assert(abs(self.goldmasked2d - matrix.col(centroid.mean())) < self.EPS)
    assert(abs(self.goldmasked2dvar - matrix.col(centroid.variance())) < self.EPS)
    assert(abs(self.goldmasked2dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS)
    print 'OK'

  def tst_centroid_masked_image3d(self):
    from dials.algorithms.image.centroid import centroid_image
    from scitbx import matrix

    centroid = centroid_image(self.pixels3d, self.mask3d)

    assert(abs(self.goldmasked3d - matrix.col(centroid.mean())) < self.EPS)
    assert(abs(self.goldmasked3dvar - matrix.col(centroid.variance())) < self.EPS)
    assert(abs(self.goldmasked3dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS)
    print 'OK'

  def tst_centroid_bias(self):

    from dials.algorithms.image.centroid import centroid_image
    from scitbx.array_family import flex
    pixels = flex.double(flex.grid(5,5), 0)
    pixels[2,2] = 10
    centroid = centroid_image(pixels)
    #assert abs(centroid.average_bias_estimate()[0] - 1/12.0) < 1e-7
    #assert abs(centroid.average_bias_estimate()[1] - 1/12.0) < 1e-7

    pixels = flex.double(flex.grid(5,5), 0)
    pixels[1,2] = 5
    pixels[2,2] = 10
    pixels[3,2] = 5
    pixels[2,1] = 5
    pixels[2,3] = 5
    centroid = centroid_image(pixels)
    assert centroid.average_bias_estimate()[0] < 1e-7
    assert centroid.average_bias_estimate()[1] < 1e-7

    print 'OK'

  def generate_data(self):
    from scitbx.array_family import flex
    from random import random, randint

    # Generate a 3d array of pixels and points
    self.points3d = flex.vec3_double(flex.grid(5, 5, 5))
    self.pixels3d = flex.double(flex.grid(5, 5, 5))
    self.mask3d = flex.bool(flex.grid(5, 5, 5))
    for k in range(0, 5):
      for j in range(0, 5):
        for i in range(0, 5):
          self.points3d[k,j,i] = (i + 0.5, j + 0.5, k + 0.5)
          self.pixels3d[k,j,i] = random()
          self.mask3d[k,j,i] = bool(randint(0, 1))

    self.points2d = flex.vec2_double(flex.grid(5, 5))
    self.pixels2d = flex.double(flex.grid(5, 5))
    self.mask2d = flex.bool(flex.grid(5, 5))
    for j in range(0, 5):
      for i in range(0, 5):
        self.points2d[j,i] = self.points3d[0, j, i][0:2]
        self.pixels2d[j,i] = self.pixels3d[0, j, i]
        self.mask2d[j,i] = self.mask3d[0, j, i]

  def calculate_gold(self):
    self.calculate_gold2d()
    self.calculate_gold3d()
    self.calculate_gold_masked2d()
    self.calculate_gold_masked3d()

  def calculate_gold2d(self):

    from scitbx.array_family import flex
    from scitbx import matrix

    r_tot = 0.0
    c_tot = 0.0
    d_tot = 0.0

    for (r, c), d in zip(self.points2d, self.pixels2d):
      r_tot += d * r
      c_tot += d * c
      d_tot += d

    self.gold2d = matrix.col((r_tot / d_tot, c_tot / d_tot))
    _r, _c = self.gold2d

    r_tot = 0.0
    c_tot = 0.0

    for (r, c), d in zip(self.points2d, self.pixels2d):
      r_tot += d * (r - _r) ** 2
      c_tot += d * (c - _c) ** 2

    _sr = r_tot / d_tot
    _sc = c_tot / d_tot

    self.gold2dvar = matrix.col((_sr, _sc))

    #r_tot = 0.0
    #c_tot = 0.0

    #for (r, c), d in zip(self.points2d, self.pixels2d):
    #  r_tot += d * (r - _r) ** 2
    #  c_tot += d * (c - _c) ** 2

    #_sr = r_tot / (d_tot-1)
    #_sc = c_tot / (d_tot-1)

    #self.gold2dubvar = matrix.col((_sr, _sc))

    pixel_x, pixel_y = zip(*self.points2d)
    xc = flex.mean_and_variance(flex.double(pixel_x), self.pixels2d.as_1d())
    yc = flex.mean_and_variance(flex.double(pixel_y), self.pixels2d.as_1d())
    self.gold2dubvar = matrix.col((xc.gsl_stats_wvariance(),
                                   yc.gsl_stats_wvariance()))

  def calculate_gold3d(self):

    from scitbx.array_family import flex
    from scitbx import matrix

    f_tot = 0.0
    r_tot = 0.0
    c_tot = 0.0
    d_tot = 0.0

    for (f, r, c), d in zip(self.points3d, self.pixels3d):
      f_tot += d * f
      r_tot += d * r
      c_tot += d * c
      d_tot += d

    self.gold3d = matrix.col((f_tot / d_tot, r_tot / d_tot, c_tot / d_tot))

    _f, _r, _c = self.gold3d

    f_tot = 0.0
    r_tot = 0.0
    c_tot = 0.0

    for (f, r, c), d in zip(self.points3d, self.pixels3d):
      f_tot += d * (f - _f) ** 2
      r_tot += d * (r - _r) ** 2
      c_tot += d * (c - _c) ** 2

    _sf = f_tot / d_tot
    _sr = r_tot / d_tot
    _sc = c_tot / d_tot

    self.gold3dvar = matrix.col((_sf, _sr, _sc))

    #f_tot = 0.0
    #r_tot = 0.0
    #c_tot = 0.0

    #for (f, r, c), d in zip(self.points3d, self.pixels3d):
    #  f_tot += d * (f - _f) ** 2
    #  r_tot += d * (r - _r) ** 2
    #  c_tot += d * (c - _c) ** 2

    #_sf = f_tot / (d_tot-1)
    #_sr = r_tot / (d_tot-1)
    #_sc = c_tot / (d_tot-1)

    #self.gold3dubvar = matrix.col((_sf, _sr, _sc))

    pixel_x, pixel_y, pixel_z = zip(*self.points3d)
    xc = flex.mean_and_variance(flex.double(pixel_x), self.pixels3d.as_1d())
    yc = flex.mean_and_variance(flex.double(pixel_y), self.pixels3d.as_1d())
    zc = flex.mean_and_variance(flex.double(pixel_z), self.pixels3d.as_1d())
    self.gold3dubvar = matrix.col((xc.gsl_stats_wvariance(),
                                   yc.gsl_stats_wvariance(),
                                   zc.gsl_stats_wvariance()))

  def calculate_gold_masked2d(self):

    from scitbx.array_family import flex
    from scitbx import matrix

    r_tot = 0.0
    c_tot = 0.0
    d_tot = 0.0

    for (r, c), d, m in zip(self.points2d, self.pixels2d, self.mask2d):
      if m:
        r_tot += d * r
        c_tot += d * c
        d_tot += d

    self.goldmasked2d = matrix.col((r_tot / d_tot, c_tot / d_tot))

    _r, _c = self.goldmasked2d

    r_tot = 0.0
    c_tot = 0.0

    for (r, c), d, m in zip(self.points2d, self.pixels2d, self.mask2d):
      if m:
        r_tot += d * (r - _r) ** 2
        c_tot += d * (c - _c) ** 2

    _sr = r_tot / d_tot
    _sc = c_tot / d_tot

    self.goldmasked2dvar = matrix.col((_sr, _sc))

    #r_tot = 0.0
    #c_tot = 0.0

    #for (r, c), d, m in zip(self.points2d, self.pixels2d, self.mask2d):
    #  if m:
    #    r_tot += d * (r - _r) ** 2
    #    c_tot += d * (c - _c) ** 2

    #_sr = r_tot / (d_tot-1)
    #_sc = c_tot / (d_tot-1)

    #self.goldmasked2dubvar = matrix.col((_sr, _sc))

    pixel_x = flex.double()
    pixel_y = flex.double()
    pixel_d = flex.double()
    for (x, y), d, m in zip(self.points2d, self.pixels2d, self.mask2d):
      if m:
        pixel_x.append(x)
        pixel_y.append(y)
        pixel_d.append(d)

    xc = flex.mean_and_variance(flex.double(pixel_x), pixel_d)
    yc = flex.mean_and_variance(flex.double(pixel_y), pixel_d)
    self.goldmasked2dubvar = matrix.col((xc.gsl_stats_wvariance(),
                                         yc.gsl_stats_wvariance()))

  def calculate_gold_masked3d(self):

    from scitbx.array_family import flex
    from scitbx import matrix

    f_tot = 0.0
    r_tot = 0.0
    c_tot = 0.0
    d_tot = 0.0

    for (f, r, c), d, m in zip(self.points3d, self.pixels3d, self.mask3d):
      if m:
        f_tot += d * f
        r_tot += d * r
        c_tot += d * c
        d_tot += d

    self.goldmasked3d = matrix.col((f_tot / d_tot,
        r_tot / d_tot, c_tot / d_tot))

    _f, _r, _c = self.goldmasked3d

    f_tot = 0.0
    r_tot = 0.0
    c_tot = 0.0

    for (f, r, c), d, m in zip(self.points3d, self.pixels3d, self.mask3d):
      if m:
        f_tot += d * (f - _f) ** 2
        r_tot += d * (r - _r) ** 2
        c_tot += d * (c - _c) ** 2

    _sf = f_tot / d_tot
    _sr = r_tot / d_tot
    _sc = c_tot / d_tot

    self.goldmasked3dvar = matrix.col((_sf, _sr, _sc))

    #f_tot = 0.0
    #r_tot = 0.0
    #c_tot = 0.0

    #for (f, r, c), d, m in zip(self.points3d, self.pixels3d, self.mask3d):
    #  if m:
    #    f_tot += d * (f - _f) ** 2
    #    r_tot += d * (r - _r) ** 2
    #    c_tot += d * (c - _c) ** 2

    #_sf = f_tot / (d_tot-1)
    #_sr = r_tot / (d_tot-1)
    #_sc = c_tot / (d_tot-1)

    #self.goldmasked3dubvar = matrix.col((_sf, _sr, _sc))

    pixel_x = flex.double()
    pixel_y = flex.double()
    pixel_z = flex.double()
    pixel_d = flex.double()
    for (x, y, z), d, m in zip(self.points3d, self.pixels3d, self.mask3d):
      if m:
        pixel_x.append(x)
        pixel_y.append(y)
        pixel_z.append(z)
        pixel_d.append(d)

    xc = flex.mean_and_variance(flex.double(pixel_x), pixel_d)
    yc = flex.mean_and_variance(flex.double(pixel_y), pixel_d)
    zc = flex.mean_and_variance(flex.double(pixel_z), pixel_d)
    self.goldmasked3dubvar = matrix.col((xc.gsl_stats_wvariance(),
                                         yc.gsl_stats_wvariance(),
                                         zc.gsl_stats_wvariance()))


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = CentroidTest()
    test.run()
