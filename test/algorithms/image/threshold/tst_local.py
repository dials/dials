from __future__ import absolute_import, division


class Test:

  def __init__(self):
    from scitbx.array_family import flex

    # Create an image
    self.image = flex.random_double(2000 * 2000, 10)
    self.image.reshape(flex.grid(2000, 2000))
    self.mask = flex.random_bool(2000 * 2000, 0.99)
    self.mask.reshape(flex.grid(2000, 2000))
    self.gain = flex.random_double(2000 * 2000) + 0.5
    self.gain.reshape(flex.grid(2000, 2000))
    self.size = (3, 3)
    self.min_count = 2

  def run(self):
    self.tst_niblack()
    self.tst_sauvola()
    self.tst_fano()
    self.tst_fano_masked()
    self.tst_gain()
    self.tst_kabsch()
    self.tst_kabsch_w_gain()
    self.tst_kabsch_debug()
    self.tst_dispersion_threshold()

  def tst_niblack(self):
    from dials.algorithms.image.threshold import niblack
    n_sigma = 3
    result = niblack(self.image, self.size, n_sigma)
    print 'OK'

  def tst_sauvola(self):
    from dials.algorithms.image.threshold import sauvola
    k = 3
    r = 3
    result = sauvola(self.image, self.size, k, r)
    print 'OK'

  def tst_fano(self):
    from dials.algorithms.image.threshold import fano
    n_sigma = 3
    result = fano(self.image, self.size, n_sigma)
    print 'OK'

  def tst_fano_masked(self):
    from dials.algorithms.image.threshold import fano_masked
    n_sigma = 3
    result = fano_masked(self.image, self.mask, self.size, self.min_count, n_sigma)
    print 'OK'

  def tst_gain(self):
    from dials.algorithms.image.threshold import gain
    n_sigma = 3
    result = gain(self.image, self.mask, self.gain, self.size, self.min_count, n_sigma)
    print 'OK'

  def tst_kabsch(self):
    from dials.algorithms.image.threshold import kabsch
    nsig_b = 3
    nsig_s = 3
    result = kabsch(self.image, self.mask, self.size, nsig_b, nsig_s, self.min_count)
    print 'OK'

  def tst_kabsch_w_gain(self):
    from dials.algorithms.image.threshold import kabsch_w_gain
    nsig_b = 3
    nsig_s = 3
    result1 = kabsch_w_gain(self.image, self.mask, self.gain, self.size, nsig_b, nsig_s, self.min_count)

    # scaling both the image and the gain should not affect the result
    result2 = kabsch_w_gain(2.*self.image, self.mask, 2.*self.gain, self.size,
        nsig_b, nsig_s, self.min_count)
    assert (result1 == result2)

    # should get the same result as kabsch if the gain is unity
    from dials.algorithms.image.threshold import kabsch
    result3 = kabsch(self.image, self.mask, self.size, nsig_b, nsig_s,
        self.min_count)
    result4 = kabsch_w_gain(self.image, self.mask, (0*self.gain+1), self.size,
        nsig_b, nsig_s, self.min_count)
    assert (result3 == result4)
    print 'OK'

  def tst_kabsch_debug(self):
    from dials.algorithms.image.threshold import kabsch
    from dials.algorithms.image.threshold import kabsch_w_gain
    from dials.algorithms.image.threshold import KabschDebug
    nsig_b = 3
    nsig_s = 3
    result1 = kabsch(self.image, self.mask, self.size, nsig_b, nsig_s, self.min_count)
    debug = KabschDebug(self.image, self.mask, self.size, nsig_b, nsig_s,  0, self.min_count)
    result2 = debug.final_mask()
    assert(result1.all_eq(result2))

    result3 = kabsch_w_gain(self.image, self.mask, self.gain, self.size,
        nsig_b, nsig_s, self.min_count)
    debug = KabschDebug(self.image, self.mask, self.gain, self.size, nsig_b,
        nsig_s,  0, self.min_count)
    result4 = debug.final_mask()
    assert (result3 == result4)
    print 'OK'

  def tst_dispersion_threshold(self):
    from dials.algorithms.image.threshold import kabsch, kabsch_w_gain
    from dials.algorithms.image.threshold import DispersionThreshold
    from dials.array_family import flex
    nsig_b = 3
    nsig_s = 3
    algorithm = DispersionThreshold(
      self.image.all(),
      self.size,
      nsig_b,
      nsig_s,
      0,
      self.min_count)

    result1 = kabsch(
      self.image,
      self.mask,
      self.size,
      nsig_b,
      nsig_s,
      self.min_count)

    result2 = kabsch_w_gain(
      self.image,
      self.mask,
      self.gain,
      self.size,
      nsig_b,
      nsig_s,
      self.min_count)

    result3 = flex.bool(flex.grid(self.image.all()))
    result4 = flex.bool(flex.grid(self.image.all()))
    algorithm(self.image, self.mask, result3)
    algorithm(self.image, self.mask, self.gain, result4)
    assert(result1 == result3)
    assert(result2 == result4)

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
