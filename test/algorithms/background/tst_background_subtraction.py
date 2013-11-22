from __future__ import division

class TestForStrong:

  def __init__(self, path):
    from iotbx.xds import integrate_hkl
    from dials.model.serialize import load
    import os

    # Load the sweep
    sweep = load.sweep(os.path.join(path, 'sweep.json'))

    # Read the INTEGRATE.HKL file
    reader = integrate_hkl.reader()
    reader.read_file(os.path.join(path, 'INTEGRATE.HKL'))

    # Sort the reflections by I/SIGMA
    iobs = reader.iobs
    sigma = reader.sigma
    ios = [ i / s for i, s in zip(iobs, sigma) ]
    index = sorted(list(range(len(iobs))), key=lambda i: ios[i], reverse=True)

    # Number to select
    n_select = len(index)

    # Get the profiles of the top n reflections
    self.profiles = []
    data = sweep.to_array()
    for i in index[0:n_select]:
      x, y, z = reader.xyzcal[i]
      x0, x1 = int(x) - 5, int(x) + 5
      y0, y1 = int(y) - 5, int(y) + 5
      z0, z1 = int(z) - 5, int(z) + 5
      if z0 < 0: z0 = 0
      if z1 > 9: z1 = 9
      profile = data[z0:z1, y0:y1, x0:x1]
      self.profiles.append(profile)

  def run(self):

    from dials.algorithms.background import XdsSubtractorAlgorithm
    from scitbx.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    create_background = XdsSubtractorAlgorithm(10, 3)

    for prof1 in self.profiles:

      #prof1 = self.profiles[0]
      mask = flex.int(flex.grid(prof1.all()), MaskCode.Background | MaskCode.Valid)
      background = create_background(prof1.as_double(), mask)
      print background

#            from matplotlib import pylab
#            from scitbx.array_family import flex
#    #        vmin, vmax = flex.min(self.profiles[0]), flex.max(self.profiles[0])
#
#            p = prof1.as_numpy_array()
#            m = mask.as_numpy_array()
#            v = p * (m == 3)
#            print flex.mean(flex.int(v).as_double())
#            print flex.mean((prof1 * (mask == 3).as_int()).as_double())
#            for i in range(prof1.all()[0]):
#                im = flex.int(prof1.as_numpy_array()[i])
#                ma = mask.as_numpy_array()[i]
#                print ma
#                ma = flex.int(ma)
#                for j in range(len(ma)):
#                    if ma[j] & MaskCode.Foreground:
#                        im[j] = int(background)

  #            print im.as_numpy_array()
  #            print ma.as_numpy_array()
  ##            pylab.imshow(im, vmin=vmin, vmax=vmax, interpolation='none')
  ##            pylab.show()

class Test:

  def __init__(self):
    import libtbx.load_env
    import os
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    path = os.path.join(dials_regression, 'background_test_data')

    self.test_for_strong = TestForStrong(path)

  def run(self):

    self.test_for_strong.run()



if __name__ == '__main__':
  test = Test()
  test.run()
