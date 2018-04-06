from __future__ import absolute_import, division, print_function

import os
import procrunner

def test_spots_xds(tmpdir):
  tmpdir.chdir()

  xds_input = 'SPOT.XDS'
  output_pickle = "spot.pickle"

  with open(xds_input, "wb") as f:
    f.write("""\
 2411.40 1000.70 25.00 16384. 0 0 0
 1328.60 2170.40 20.57 7326. 0 0 0
 177.56 2191.30 24.94 6779. 0 0 0
 1231.34 1220.04 24.99 1952. 0 0 0
 1227.07 1230.56 24.81 349. 0 0 0
 1341.63 1243.25 5.64 321. 2 -2 11
 1125.23 1197.72 12.14 231. -1 2 -10
 1317.52 1171.59 19.28 120. 6 -4 6
 1260.25 1300.55 13.67 116. -4 2 6
 1090.27 1199.47 41.49 114. -2 3 -13
""")

  result = procrunner.run_process([
      "dials.import_xds",
      xds_input, #xparm_file,
      "input.method=reflections",
      "output.filename="+output_pickle,
      "remove_invalid=True",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists(output_pickle)

  import cPickle as pickle
  with open(output_pickle, "rb") as f:
    reflections = pickle.load(f)
  assert len(reflections) == 5

  os.remove(xds_input)
  assert not os.path.exists(xds_input)

  # now test we can export it again
  result = procrunner.run_process([
      "dials.export",
      "format=xds",
      output_pickle,
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists(os.path.join("xds", "SPOT.XDS"))

  with open(os.path.join("xds", "SPOT.XDS"), "rb") as f:
    txt = f.read()
  assert [line.strip() for line in txt.split('\n')] == \
         [line.strip() for line in """\
 1341.63 1243.25 5.64 321.00  2 -2 11
 1125.23 1197.72 12.14 231.00  -1 2 -10
 1317.52 1171.59 19.28 120.00  6 -4 6
 1260.25 1300.55 13.67 116.00  -4 2 6
 1090.27 1199.47 41.49 114.00  -2 3 -13
""".split('\n')]


def test_export_xds(dials_regression, tmpdir):
  tmpdir.chdir()

  result = procrunner.run_process([
      "dials.find_spots",
      os.path.join(dials_regression, "centroid_test_data", "datablock.json"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('strong.pickle')

  result = procrunner.run_process([
      "dials.export",
      "format=xds",
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
      "strong.pickle",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("xds/XDS.INP")
  assert os.path.exists("xds/XPARM.XDS")
  assert os.path.exists("xds/SPOT.XDS")

  os.remove("xds/XDS.INP")
  os.remove("xds/XPARM.XDS")
  assert not os.path.exists("xds/XDS.INP")
  assert not os.path.exists("xds/XPARM.XDS")
  result = procrunner.run_process([
      "dials.export",
      "format=xds",
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("xds/XDS.INP")
  assert os.path.exists("xds/XPARM.XDS")
