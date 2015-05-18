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

    self.path = join(dials_regression, "refinement_test_data", "i04_weak_data")

  def run(self):
    from os.path import join
    from libtbx import easy_run
    from dials.algorithms.profile_model.factory import phil_scope
    from libtbx.phil import parse

    # Call dials.create_profile_model
    easy_run.fully_buffered([
      'dials.create_profile_model',
      join(self.path, 'experiments.json'),
      join(self.path, 'indexed_strong3.pickle'),
    ]).raise_if_errors()

    with open('profile.phil', 'r') as infile:
      text = '\n'.join(infile.readlines())
      params = phil_scope.fetch(parse(text)).extract()
      assert(params.profile.gaussian_rs.scan_varying == False)
      assert(len(params.profile.gaussian_rs.model) == 1)
      sigma_b = params.profile.gaussian_rs.model[0].sigma_b
      sigma_m = params.profile.gaussian_rs.model[0].sigma_m
      eps = 1e-6
      try:
        assert(abs(sigma_b - 0.02262206634) < eps)
        assert(abs(sigma_m - 0.0774177) < eps)
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
