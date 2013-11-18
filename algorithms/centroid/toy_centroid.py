# toy centroid implementation: this is deliberately bad to illustrate how the
# interface works and how things could be improved
from __future__ import division
from dials.interfaces.centroid.centroid_interface_prototype import \
    centroid_interface_prototype as centroid_interface
from dials.interfaces.centroid.centroid_interface_prototype import \
    CentroidException

class toy_centroid(centroid_interface):
  def __init__(self, reflections):

    centroid_interface.__init__(self, reflections)

    return

  def compute_shoebox_centroid(self, shoebox):
    from dials.algorithms.image.centroid import centroid_points
    from scitbx.array_family import flex
    import math

    f_size, r_size, c_size = shoebox.all()

    # build the list of pixels - let's be dumb and just have a literal
    # list - and assign density of a pixel to the centre of the
    # volume...

    points = flex.vec3_double()
    pixels = flex.double()

    for i in shoebox:
      if i < 0:
        raise CentroidException, 'negative pixels in cube'

    try:
      for f in range(f_size):
        for r in range(r_size):
          for c in range(c_size):
            points.append((f + 0.5, r + 0.5, c + 0.5))
            pixels.append(shoebox[f, r, c])

    except IndexError:
      raise CentroidException, 'data outside range'

    centroid = centroid_points(pixels, points)
    _f, _r, _c = centroid.mean()
    _sf, _sr, _sc = centroid.unbiased_variance()
    _d = centroid.sum_pixels()

    return _f, _r, _c, _sf, _sr, _sc, _d
