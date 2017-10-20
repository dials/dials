from __future__ import absolute_import, division
from dials.array_family import flex # import dependency


class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'SKIP: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")
    self.integration_test_data = join(dials_regression, "integration_test_data")

  def run(self):

    from os.path import join, exists
    from libtbx import easy_run

    # Run a few commands from stdin
    stdin_lines = [
      "import template=%s" % join(self.path, "centroid_####.cbf"),
      "find_spots",
      "discover_better_experimental_model",
      "index",
      "refine_bravais_settings",
      "reindex solution=22",
      "refine",
      "goto 6",
    ]

    easy_run.fully_buffered(
      'idials',
      stdin_lines=stdin_lines).raise_if_errors()
    print 'OK'

    # Check that state works
    stdin_lines = [
      "refine",
      "integrate profile.fitting=False",
      "export ignore_profile_fitting=True keep_partials=True include_partials=True",
      "goto 7",
      "integrate profile.fitting=False",
      "export ignore_profile_fitting=True keep_partials=True include_partials=True",
    ]

    easy_run.fully_buffered('idials',
                            stdin_lines=stdin_lines).raise_if_errors()

    print 'OK'

    # Check all the stuff we expect, exists
    assert exists("dials.state")
    assert exists("dials-1")
    assert exists("10_integrated.mtz")
    assert exists("12_integrated.mtz")
    assert exists("dials-1/1_import")
    assert exists("dials-1/2_find_spots")
    assert exists("dials-1/3_discover_better_experimental_model")
    assert exists("dials-1/4_index")
    assert exists("dials-1/5_refine_bravais_settings")
    assert exists("dials-1/6_reindex")
    assert exists("dials-1/7_refine")
    assert exists("dials-1/8_refine")
    assert exists("dials-1/9_integrate")
    assert exists("dials-1/10_export")
    assert exists("dials-1/11_integrate")
    assert exists("dials-1/12_export")

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
