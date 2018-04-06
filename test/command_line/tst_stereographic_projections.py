from __future__ import absolute_import, division

def run():
  import os
  import libtbx.load_env
  from libtbx import easy_run
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(dials_regression, "experiment_test_data")
  cmd = " ".join((
    "dials.stereographic_projection",
    "%s/experiment_1.json" %path,
    "hkl_limit=4",
    "plot.show=False",
    "plot.filename=proj.png",
    "json.filename=proj.json"))
  print cmd
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert os.path.exists("projections.txt")
  assert os.path.exists("proj.png")
  assert os.path.exists("proj.json")

  import json
  with open("proj.json", "rb") as f:
    d = json.load(f)
    assert d.keys() == ['data', 'layout']
    assert d['data'][0]['name'] == 'stereographic_projections'
    assert len(d['data'][0]['x']) == 578

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
