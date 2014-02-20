from __future__ import division


class Test(object):

  def __init__(self):
    from os.path import join
    from dials.array_family import flex
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "refinement_test_data", "i04_weak_data")

  def run(self):
    from os.path import abspath, join
    from libtbx import easy_run
    import os
    from uuid import uuid4

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.create_profile_model',
      join(self.path, 'experiments.json'),
      join(self.path, 'indexed_strong3.pickle'),
    ]).raise_if_errors()

    with open('profile.phil', 'r') as infile:
      for line in infile.readlines():
        tokens = line.split("=")
        if "sigma_b" in tokens[0]:
          sigma_b = float(tokens[1])
        if "sigma_m" in tokens[0]:
          sigma_m = float(tokens[1])
      eps = 1e-7
      assert(abs(sigma_b - 0.02229742649) < eps)
      assert(abs(sigma_m - 0.07971599709) < eps)
    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
