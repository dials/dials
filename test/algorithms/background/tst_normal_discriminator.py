from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_no_mask()

  def is_correct(self, shoebox, mask, n_sigma, min_pixels):
    from scitbx.array_family import flex
    flags = [m & (1 << 0) for m in mask]
    pixels = flex.select(shoebox, flags=flags)
    assert(len(pixels) >= min_pixels)
    meanp = flex.mean_and_variance(flex.double(pixels))
    maxp = flex.max(flex.double(pixels))
    m = meanp.mean()
    s = meanp.unweighted_sample_standard_deviation()
    assert(maxp <= m + n_sigma * s)

  def tst_no_mask(self):
    from scitbx.array_family import flex
    from dials.algorithms.background import NormalDiscriminator
    discriminate = NormalDiscriminator(min_data=10)
    shoebox_d = flex.random_double(5 * 5 * 5) * 100
    shoebox = flex.int([int(s) for s in shoebox_d])
    shoebox.reshape(flex.grid((5, 5, 5)))
    mask = discriminate(shoebox)
    self.is_correct(shoebox, mask, 3.0, 10)

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
