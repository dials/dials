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
    from math import floor

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
    for i in range(len(rtable)):
      from dials.algorithms.background import PlaneModel
      xdet, ydet = rtable[i]["xy"]
      nx = rtable[i]['nx']
      ny = rtable[i]['ny']
      nc = rtable[i]['nc']
      nrx = rtable[i]['nrx']
      nry = rtable[i]['nry']
      bbox = rtable[i]['bbox']
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
    assert(count == 1)
    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
