from __future__ import absolute_import, division

class Test:

  def __init__(self):
    pass

  def run(self):
    from dials.algorithms.image.filter import fano_filter
    from scitbx.array_family import flex
    from random import randint

    # Create an image
    image = flex.random_double(2000 * 2000)
    image.reshape(flex.grid(2000, 2000))
    mask = flex.random_bool(2000 * 2000, 0.99).as_int()
    mask.reshape(flex.grid(2000, 2000))

    # Calculate the summed area table
    mask2 = mask.deep_copy()
    fano_filter = fano_filter(image, mask2, (3, 3), 2)
    mean = fano_filter.mean()
    var = fano_filter.sample_variance()
    fano = fano_filter.fano()

    # For a selection of random points, ensure that the value is the
    # sum of the area under the kernel
    eps = 1e-7
    for i in range(10000):
      i = randint(10, 1990)
      j = randint(10, 1990)
      m1 = mean[j,i]
      v1 = var[j,i]
      f1 = fano[j,i]
      p = image[j-3:j+4,i-3:i+4]
      m = mask[j-3:j+4,i-3:i+4]
      if mask[j,i] == 0:
        m2 = 0.0
        v2 = 0.0
        f2 = 1.0
      else:
        p = flex.select(p, flags=m)
        mv = flex.mean_and_variance(flex.double(p))
        m2 = mv.mean()
        v2 = mv.unweighted_sample_variance()
        f2 = v2 / m2
      assert(abs(m1 - m2) <= eps)
      assert(abs(v1 - v2) <= eps)
      assert(abs(f1 - f2) <= eps)

    # Test passed
    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
