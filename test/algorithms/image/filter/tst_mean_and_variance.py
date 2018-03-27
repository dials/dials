from __future__ import absolute_import, division, print_function

class Test:

  def __init__(self):
    pass

  def run(self):
    self.tst_mean_filter()
    self.tst_masked_mean_filter()
    self.tst_mean_and_variance_filter()
    self.tst_masked_mean_and_variance_filter()

  def tst_mean_filter(self):
    from dials.algorithms.image.filter import mean_filter
    from scitbx.array_family import flex
    from random import randint

    # Create an image
    image = flex.random_double(2000 * 2000)
    image.reshape(flex.grid(2000, 2000))

    # Calculate the summed area table
    mean = mean_filter(image, (3, 3))

    # For a selection of random points, ensure that the value is the
    # sum of the area under the kernel
    eps = 1e-7
    for i in range(10000):
      i = randint(10, 1990)
      j = randint(10, 1990)
      m1 = mean[j,i]
      p = image[j-3:j+4,i-3:i+4]
      mv = flex.mean_and_variance(p.as_1d())
      m2 = mv.mean()
      assert(abs(m1 - m2) <= eps)

    # Test passed

  def tst_masked_mean_filter(self):
    from dials.algorithms.image.filter import mean_filter
    from scitbx.array_family import flex
    from random import randint

    # Create an image
    image = flex.random_double(2000 * 2000)
    image.reshape(flex.grid(2000, 2000))
    mask = flex.random_bool(2000 * 2000, 0.99).as_int()
    mask.reshape(flex.grid(2000, 2000))

    # Calculate the summed area table
    mask2 = mask.deep_copy()
    mean = mean_filter(image, mask2, (3, 3), 1)

    # For a selection of random points, ensure that the value is the
    # sum of the area under the kernel
    eps = 1e-7
    for i in range(10000):
      i = randint(10, 1990)
      j = randint(10, 1990)
      m1 = mean[j,i]
      p = image[j-3:j+4,i-3:i+4]
      m = mask[j-3:j+4,i-3:i+4]
      if mask[j,i] == 0:
        m2 = 0.0
      else:
        p = flex.select(p, flags=m)
        mv = flex.mean_and_variance(flex.double(p))
        m2 = mv.mean()
        s1 = flex.sum(flex.double(p))
        s2 = flex.sum(m.as_1d())
      assert(abs(m1 - m2) <= eps)

    # Test passed

  def tst_mean_and_variance_filter(self):
    from dials.algorithms.image.filter import mean_and_variance_filter
    from scitbx.array_family import flex
    from random import randint

    # Create an image
    image = flex.random_double(2000 * 2000)
    image.reshape(flex.grid(2000, 2000))

    # Calculate the summed area table
    mean_and_variance = mean_and_variance_filter(image, (3, 3))
    mean = mean_and_variance.mean()
    variance = mean_and_variance.variance()
    sample_variance = mean_and_variance.sample_variance()

    # For a selection of random points, ensure that the value is the
    # sum of the area under the kernel
    eps = 1e-7
    for i in range(10000):
      i = randint(10, 1990)
      j = randint(10, 1990)
      m1 = mean[j,i]
      v1 = variance[j,i]
      sv1 = sample_variance[j,i]
      p = image[j-3:j+4,i-3:i+4]
      mv = flex.mean_and_variance(p.as_1d())
      m2 = mv.mean()
      sv2 = mv.unweighted_sample_variance()
      assert(abs(m1 - m2) <= eps)
      assert(abs(sv1 - sv2) <= eps)

    # Test passed

  def tst_masked_mean_and_variance_filter(self):
    from dials.algorithms.image.filter import mean_and_variance_filter
    from scitbx.array_family import flex
    from random import randint

    # Create an image
    image = flex.random_double(2000 * 2000)
    image.reshape(flex.grid(2000, 2000))
    mask = flex.random_bool(2000 * 2000, 0.99).as_int()
    mask.reshape(flex.grid(2000, 2000))

    # Calculate the summed area table
    mask2 = mask.deep_copy()
    mv = mean_and_variance_filter(image, mask2, (3, 3), 2)
    mean = mv.mean()
    var = mv.sample_variance()

    # For a selection of random points, ensure that the value is the
    # sum of the area under the kernel
    eps = 1e-7
    for i in range(10000):
      i = randint(10, 1990)
      j = randint(10, 1990)
      m1 = mean[j,i]
      v1 = var[j,i]
      p = image[j-3:j+4,i-3:i+4]
      m = mask[j-3:j+4,i-3:i+4]
      if mask[j,i] == 0:
        m2 = 0.0
        v2 = 0.0
      else:
        p = flex.select(p, flags=m)
        mv = flex.mean_and_variance(flex.double(p))
        m2 = mv.mean()
        v2 = mv.unweighted_sample_variance()
      assert(abs(m1 - m2) <= eps)
      assert(abs(v1 - v2) <= eps)

    # Test passed


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
