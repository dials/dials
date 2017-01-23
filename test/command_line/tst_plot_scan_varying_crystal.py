from __future__ import absolute_import, division
from dials.array_family import flex # import dependency

class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(
      dials_regression,
      "refinement_test_data",
      "multi_sweep_one_sample",
      "glucose_isomerase",
      "SWEEP1",
      "index",
      "sv_refined_experiments.json")

  def run(self):
    from libtbx import easy_run

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.plot_scan_varying_crystal',
      self.path
    ]).raise_if_errors()

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
