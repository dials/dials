from __future__ import absolute_import, division, print_function

import os

def test(dials_regression, tmpdir):
  tmpdir.chdir()

  from libtbx import easy_run

  # Run a few commands from stdin
  stdin_lines = [
    "import template=" + os.path.join(dials_regression, "centroid_test_data", "centroid_####.cbf"),
    "find_spots",
    "search_beam_position",
    "index",
    "refine_bravais_settings",
    "reindex solution=22",
    "refine",
    "goto 6",
  ]

  easy_run.fully_buffered(
    'idials',
    stdin_lines=stdin_lines).raise_if_errors()

  # Check that state works
  stdin_lines = [
    "refine",
    "integrate profile.fitting=False",
    "export intensity='sum'",
    "goto 7",
    "integrate profile.fitting=False",
    "export intensity='sum'",
  ]

  easy_run.fully_buffered('idials',
                          stdin_lines=stdin_lines).raise_if_errors()

  # Check all the stuff we expect, exists
  assert os.path.exists("dials.state")
  assert os.path.exists("dials-1")
  assert os.path.exists("10_integrated.mtz")
  assert os.path.exists("12_integrated.mtz")
  assert os.path.exists("dials-1/1_import")
  assert os.path.exists("dials-1/2_find_spots")
  assert os.path.exists("dials-1/3_search_beam_position")
  assert os.path.exists("dials-1/4_index")
  assert os.path.exists("dials-1/5_refine_bravais_settings")
  assert os.path.exists("dials-1/6_reindex")
  assert os.path.exists("dials-1/7_refine")
  assert os.path.exists("dials-1/8_refine")
  assert os.path.exists("dials-1/9_integrate")
  assert os.path.exists("dials-1/10_export")
  assert os.path.exists("dials-1/11_integrate")
  assert os.path.exists("dials-1/12_export")
