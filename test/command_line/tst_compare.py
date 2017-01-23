from __future__ import absolute_import, division


class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "integration_test_data/i04-weak-data")

  def run(self):
    from os.path import join
    from libtbx import easy_run

    # Call dials.compare_mosflm_dials
    easy_run.fully_buffered([
      'dials.compare_mosflm_dials',
      join(self.path, 'integrate.mtz'),
      join(self.path, 'integrated.pickle'),
      join(self.path, 'crystal.json'),
      join(self.path, 'sweep.json'),
    ]).raise_if_errors()


    # remember to uncomment the next line
    #assert(len(table) == 361)
    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
