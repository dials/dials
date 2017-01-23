
from __future__ import division

class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):
    from os.path import join, exists
    from libtbx import easy_run

    assert(exists(join(self.path, "datablock.json")))

    input_filename = join(self.path, "datablock.json")

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.generate_mask',
      input_filename,
    ]).raise_if_errors()

    assert(exists("mask.pickle"))

    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask2.pickle',
      'untrusted.rectangle=100,200,100,200'
    ]).raise_if_errors()
    assert(exists("mask2.pickle"))

    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask3.pickle',
      'untrusted.circle=100,100,10'
    ]).raise_if_errors()
    assert(exists("mask3.pickle"))

    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask4.pickle',
      'resolution_range=2,3',
    ]).raise_if_errors()
    assert(exists("mask4.pickle"))

    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask5.pickle',
      'd_min=3',
      'd_max=2',
    ]).raise_if_errors()
    assert(exists("mask5.pickle"))

    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask6.pickle',
      '\'ice_rings{filter=True;d_min=2}\'',
    ]).raise_if_errors()
    assert(exists("mask6.pickle"))

    print 'OK'

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask3.pickle',
      'untrusted.polygon=100,100,100,200,200,200,200,100'
    ]).raise_if_errors()
    assert(exists("mask3.pickle"))

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
