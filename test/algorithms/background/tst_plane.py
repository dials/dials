from __future__ import division

class Test(object):

  def __init__(self):

    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      return

    # The directory path
    path = join(
      dials_regression,
      "integration_test_data",
      "i04-weak-data",
      "jmp_mosflm_test")

    # The input files
    self.reflection_filename = join(path, "mosflm_reflections.pickle")
    self.shoebox_filename = join(path, "shoeboxes.pickle")

  def run(self):
    from dials.array_family import flex
    import cPickle as pickle
    from math import floor, ceil, sqrt
    from dials.algorithms.shoebox import MaskCode

    # Read the data
    rtable = flex.reflection_table.from_pickle(self.reflection_filename)
    shoeboxes, masks = pickle.load(open(self.shoebox_filename))
    assert(len(rtable) == len(shoeboxes))
    assert(len(rtable) == len(masks))

    # Compute the background for each reflection and check against the values
    # read from the mosflm.lp file. Currently this fails for 1 strange
    # reflection whose pixel values in the mosflm file do not match those
    # extracted from the images.
    count = 0
    VAR1 = []
    VAR2 = []
    DIFF = []
    for i in range(len(rtable)):
      from dials.algorithms.background import PlaneModel
      xdet, ydet = rtable[i]["xy"]
      nx = rtable[i]['nx']
      ny = rtable[i]['ny']
      nc = rtable[i]['nc']
      nrx = rtable[i]['nrx']
      nry = rtable[i]['nry']
      bbox = rtable[i]['bbox']
      I = rtable[i]['intensity.sum.value']
      Ivar = rtable[i]['intensity.sum.variance']
      lp = rtable[i]['lp']
      data = shoeboxes[i]
      mask = masks[i]
      fraction = 1.0
      nsigma = 4
      try:
        model = PlaneModel(data, mask, fraction, nsigma)
      except Exception:
        count += 1
        continue
      n = model.noutlier()
      a1 = model.a()
      b1 = model.b()
      c1 = model.c()
      a2 = rtable[i]['background'][0]
      b2 = rtable[i]['background'][1]
      c2 = rtable[i]['background'][2]
      try:
        assert(abs(a1 - b2) < 0.01)
        assert(abs(b1 + a2) < 0.01)
        assert(abs(c1 - c2) < 0.1)
      except Exception:
        count += 1
        continue
        #print "BG %d:(%.2f, %.2f, %.1f), (%.2f, %.2f, %.1f): %d" % \
          #(i, a1, b1, c1, a2, b2, c2, n)
        #print "X, Y: ", xdet, ydet
        #print "NX, NY: ", nx, ny
        #print "NRX, NRY, NC", nrx, nry, nc
        #print int(floor(xdet + 0.5)) - nx // 2, int(floor(ydet + 0.5)) - ny // 2
        #print "BBOX: ", bbox
        #print "N Outliers: ", model.noutlier()
        #print "N Background: ", model.nbackground()
        #print "Max DIff: ", model.maxdiff()
        #print data.as_numpy_array().transpose()[::-1,::-1]
        #print mask.as_numpy_array().transpose()[::-1,::-1]
        #raise

      background = data.as_double()
      hy = background.all()[0] // 2
      hx = background.all()[1] // 2
      for jj in range(background.all()[0]):
        for ii in range(background.all()[1]):
          x = ii - hx
          y = jj - hy
          background[jj,ii] = a1 * x + b1 * y + c1

      # Test the summation results. Edge reflections use profile fitted
      # intensity in MOSFLM. Therefore ignore these. There also appears to be a
      # some discrepancy with very low <= 0 reflections where an extra 0.5 is
      # added. Not sure why this is so ignore these reflections as well.
      from dials.algorithms.integration import integrate_by_summation
      intensity = integrate_by_summation(data.as_double(), background, mask)
      I2 = intensity.intensity()
      Ivar2 = intensity.variance()
      I1 = I
      Ivar1 = Ivar
      if mask.count(0) == 0 and mask.count(2) == 0 and I1 > 0:
        VAR1.append(sqrt(Ivar1))
        VAR2.append(sqrt(Ivar2))
        DIFF.append(sqrt(Ivar1) - sqrt(Ivar2))
        try:
          assert(abs(I1 - I2) < 1.0)
          assert(abs(sqrt(Ivar1) - sqrt(Ivar2)) < 1.0)
        except Exception:
          count += 1
          #import numpy
          #numpy.set_printoptions(precision=4, linewidth=200)
          #print "# %d" % i
          #print "I: %f, %f, %f" % (I2, I1, lp)
          #print "DEBUG: ", c1 * 25
          #print "PF: %f" % rtable[i]['intensity.prf.value']
          #print "BG (%.4f, %.4f, %.4f), (%.2f, %.2f, %.1f): %d" % \
            #(a1, b1, c1, a2, b2, c2, n)
          #print "X, Y: ", xdet, ydet
          #print "NX, NY: ", nx, ny
          #print "NRX, NRY, NC", nrx, nry, nc
          #print int(floor(xdet + 0.5)) - nx // 2, int(floor(ydet + 0.5)) - ny // 2
          #print "BBOX: ", bbox
          #print "N Outliers: ", model.noutlier()
          #print "N Background: ", model.nbackground()
          #print "Max DIff: ", model.maxdiff()
          #temp = (mask == MaskCode.Valid | MaskCode.Foreground).as_1d().as_int()
          #temp.resize(flex.grid(*data.all()))
          #temp = temp.as_double()
          #print data.as_numpy_array().transpose()[::-1,::-1]
          #print (background * temp).as_numpy_array().transpose()[::-1,::-1]
          #print mask.as_numpy_array().transpose()[::-1,::-1]
          #print ((data.as_double() - background) * temp).as_numpy_array().transpose()[::-1,::-1]
          #raise
          continue

    # Only 1 should fail
    assert(count == 1)
    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
