from __future__ import absolute_import, division, print_function

class Test:

  def __init__(self):
    self.kernel = (3, 3)

  def run(self):
    self.tst_filter()
    self.tst_masked_filter()

  def generate_image(self, xsize, ysize):
    from scitbx.array_family import flex
    image = flex.random_double(xsize * ysize)
    image.reshape(flex.grid(ysize, xsize))
    return image

  def generate_mask(self, xsize, ysize):
    from scitbx.array_family import flex
    mask = flex.random_bool(xsize * ysize, 0.9)
    mask.reshape(flex.grid(ysize, xsize))
    return mask

  def tst_filter(self):
    from numpy import median
    from dials.algorithms.image.filter import median_filter

    xsize = 200
    ysize = 300

    image = self.generate_image(xsize, ysize)

    result = median_filter(image, self.kernel)

    eps = 1e-7

    for j in range(self.kernel[0], ysize-self.kernel[0]):
      for i in range(self.kernel[1], xsize-self.kernel[1]):
        j0 = j - self.kernel[0]
        j1 = j + self.kernel[0] + 1
        i0 = i - self.kernel[1]
        i1 = i + self.kernel[1] + 1
        pixels = image[j0:j1,i0:i1]
        value = median(pixels.as_numpy_array())
        assert(abs(value - result[j,i]) < eps)


  def tst_masked_filter(self):

    from numpy import median
    from dials.algorithms.image.filter import median_filter

    xsize = 200
    ysize = 300

    image = self.generate_image(xsize, ysize)
    mask = self.generate_mask(xsize, ysize)
    result = median_filter(image, mask, self.kernel)
    eps = 1e-7

    for j in range(self.kernel[0], ysize-self.kernel[0]):
      for i in range(self.kernel[1], xsize-self.kernel[1]):
        if mask[j,i]:
          j0 = j - self.kernel[0]
          j1 = j + self.kernel[0] + 1
          i0 = i - self.kernel[1]
          i1 = i + self.kernel[1] + 1
          pixels = image[j0:j1,i0:i1]
          pmask = mask[j0:j1,i0:i1]
          pixels = pixels.as_1d().select(pmask.as_1d())
          if len(pixels) & 1:
            value = median(pixels.as_numpy_array())
          else:
            pixels = sorted(list(pixels))
            value = pixels[len(pixels) // 2]
          assert(abs(value - result[j,i]) < eps)
        else:
          assert(abs(result[j,i] - 0.0) < eps)


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
