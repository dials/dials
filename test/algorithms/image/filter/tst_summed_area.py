from __future__ import absolute_import, division

class Test:

  def __init__(self):
    pass

  def run(self):
    from dials.algorithms.image.filter import summed_area
    from scitbx.array_family import flex
    from random import randint

    # Create an image
    image = flex.random_double(2000 * 2000)
    image.reshape(flex.grid(2000, 2000))

    # Calculate the summed area table
    sa = summed_area(image, (3, 3))

    # For a selection of random points, ensure that the value is the
    # sum of the area under the kernel
    eps = 1e-7
    for i in range(10000):
      i = randint(10, 1990)
      j = randint(10, 1990)
      v = sa[j,i]
      e = flex.sum(image[j-3:j+4,i-3:i+4])
      assert(abs(e - v) <= eps)

    # Test passed
    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
