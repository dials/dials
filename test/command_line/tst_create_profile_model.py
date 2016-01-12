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

    self.path = join(dials_regression, "integration_test_data", "i04-weak-data2")

  def run(self):
    from os.path import join
    from libtbx import easy_run
    from dials.algorithms.profile_model.factory import phil_scope
    from libtbx.phil import parse
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory

    # Call dials.create_profile_model
    easy_run.fully_buffered([
      'dials.create_profile_model',
      join(self.path, 'experiments.json'),
      join(self.path, 'indexed.pickle'),
    ]).raise_if_errors()


    experiments =  ExperimentListFactory.from_json_file(
      "experiments_with_profile_model.json",
      check_format=False)
    sigma_b = experiments[0].profile.sigma_b(deg=True)
    sigma_m = experiments[0].profile.sigma_m(deg=True)
    eps = 1e-3
    try:
      assert(abs(sigma_b - 0.02195) < eps)
      assert(abs(sigma_m - 0.06833) < eps)
    except Exception:
      print sigma_b
      print sigma_m
      raise
    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
