from __future__ import division

class Analyse(object):

  def __init__(self, integrate_file):
    self.integrate_file = integrate_file

  def __call__(self):

    from iotbx.xds import integrate_hkl
    from scipy.interpolate import griddata
    from scitbx import matrix
    from matplotlib import pylab
    import numpy

    # Read the file
    print "Read INTEGRATE.HKL"
    handle = integrate_hkl.reader()
    handle.read_file(self.integrate_file)

    # Get the data
    xyzcal = handle.xyzcal
    xyzobs = handle.xyzobs
    iobs = handle.iobs
    sigma = handle.sigma

    width, height = handle.detector_size

    print "Get Diff arrays"
    diff = []
    x = []
    y = []
    for c, o, i, sig in zip(xyzcal, xyzobs, iobs, sigma):

      o = matrix.col(o)
      c = matrix.col(c)

      if o.length() > 0 and c[2] < 10:

        # Calculate the difference
        diff.append((c - o).length())
        x.append(c[0])
        y.append(c[1])

    print "Create grid array"
    xp = numpy.arange(width * height, dtype=numpy.int32) % width
    yp = numpy.arange(width * height, dtype=numpy.int32) / width
    points = numpy.zeros(shape=(width * height, 2), dtype=numpy.int32)
    points[:,0] = xp
    points[:,1] = yp

    print "Grid data"
    grid = griddata((x, y), diff, points, 'cubic', 0)
    grid.shape = (height, width)

    pylab.imshow(grid)
    pylab.scatter(x, y)
    pylab.show()

if __name__ == '__main__':
  import sys
  analyse = Analyse(sys.argv[1])
  analyse()
