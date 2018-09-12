from __future__ import absolute_import, division, print_function

import json
import os
from libtbx import easy_run

def test_stereographic_projectsion(dials_regression, run_in_tmpdir):
  path = os.path.join(dials_regression, "experiment_test_data")
  cmd = " ".join((
    "dials.stereographic_projection",
    "%s/experiment_1.json" %path,
    "hkl_limit=4",
    "plot.show=False",
    "plot.filename=proj.png",
    "json.filename=proj.json"))
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("projections.txt")
  assert os.path.exists("proj.png")
  assert os.path.exists("proj.json")

  with open("proj.json", "rb") as f:
    d = json.load(f)
    assert d.keys() == ['data', 'layout']
    assert d['data'][0]['name'] == 'stereographic_projections'
    assert len(d['data'][0]['x']) == 578
