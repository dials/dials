from __future__ import absolute_import, division, print_function

def test(dials_regression, tmpdir):
  from os.path import join
  path = join(dials_regression, "integration_test_data/i04-weak-data")
  tmpdir.chdir()

  from libtbx import easy_run
  # Call dev.dials.compare_mosflm_dials
  easy_run.fully_buffered([
      'dev.dials.compare_mosflm_dials',
      join(path, 'integrate.mtz'),
      join(path, 'integrated.pickle'),
      join(path, 'crystal.json'),
      join(path, 'sweep.json'),
  ]).raise_if_errors()

  # remember to uncomment the next line
  #assert(len(table) == 361)

if __name__ == '__main__':
  import libtbx.load_env
  dials_regression = libtbx.env.dist_path('dials_regression')
  import mock
  tmpdir = mock.Mock()
  from dials.test import cd_auto
  with cd_auto(__file__):
    test(dials_regression, tmpdir)
