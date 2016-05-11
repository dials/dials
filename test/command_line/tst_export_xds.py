from __future__ import division
from dials.array_family import flex # import dependency

import os
import cPickle as pickle
from libtbx.test_utils import open_tmp_directory, show_diff
from libtbx import easy_run

def exercise_spots_xds():
  txt = """\
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
"""

  tmp_dir = os.path.abspath(open_tmp_directory())
  f = open(os.path.join(tmp_dir, "SPOT.XDS"), mode="wb")
  f.write(txt)
  f.close()

  output_pickle = "%s.pickle" %f.name[:-4]
  args = ["dials.import_xds", f.name, #xparm_file,
          "input.method=reflections",
          "output.filename='%s'" %output_pickle,
          "remove_invalid=True"]
  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists(output_pickle)
  with open(output_pickle, "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 5
  os.remove(os.path.join(tmp_dir, "SPOT.XDS"))
  assert not os.path.exists(os.path.join(tmp_dir, "SPOT.XDS"))

  # now test we can export it
  args = ["dials.export format=xds", output_pickle]
  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists(os.path.join("xds", "SPOT.XDS"))
  with open(os.path.join("xds", "SPOT.XDS"), "rb") as f:
    txt = f.read()
    assert not show_diff(
      "\n".join([l.rstrip() for l in txt.split("\n")]), """\
 1341.63 1243.25 5.64 321.00  2 -2 11
 1125.23 1197.72 12.14 231.00  -1 2 -10
 1317.52 1171.59 19.28 120.00  6 -4 6
 1260.25 1300.55 13.67 116.00  -4 2 6
 1090.27 1199.47 41.49 114.00  -2 3 -13
""")


def export_xds():
  from libtbx import easy_run
  import libtbx.load_env

  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print "skipping exercise_export_xds: dials_regression not available"
    return

  cwd = os.path.abspath(os.curdir)
  tmp_dir = os.path.abspath(open_tmp_directory(suffix="export_xds"))
  os.chdir(tmp_dir)

  args = ["dials.find_spots",
          os.path.join(dials_regression, "centroid_test_data",
                       "datablock.json")]

  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  pickle_path = os.path.join(tmp_dir, "strong.pickle")
  assert os.path.exists(pickle_path)

  args = ["dials.export", "format=xds",
          os.path.join(dials_regression, "centroid_test_data",
                       "experiments.json"),
          pickle_path]
  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists("xds/XDS.INP")
  assert os.path.exists("xds/XPARM.XDS")
  assert os.path.exists("xds/SPOT.XDS")

  os.remove("xds/XDS.INP")
  os.remove("xds/XPARM.XDS")
  assert not os.path.exists("xds/XDS.INP")
  assert not os.path.exists("xds/XPARM.XDS")
  args = ["dials.export", "format=xds",
          os.path.join(dials_regression, "centroid_test_data",
                       "experiments.json")]
  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()
  assert os.path.exists("xds/XDS.INP")
  assert os.path.exists("xds/XPARM.XDS")

  os.chdir(cwd)


def run():
  export_xds()
  exercise_spots_xds()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
