from __future__ import division

def run():

  import os
  import shutil
  import fileinput
  import libtbx.load_env
  from libtbx import easy_run
  from libtbx.test_utils import show_diff
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(
    dials_regression, "centroid_test_data/centroid_####.cbf")

  print os.getcwd()

  cmd = "dials.import template=%s" %path
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.find_spots datablock.json"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.index datablock.json strong.pickle space_group=P422"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.integrate experiments.json indexed.pickle"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.export integrated_experiments.json integrated.pickle format=best"
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  assert os.path.exists("best.dat")
  with open("best.dat", "rb") as f:
    lines = ''.join(f.readlines()[:10])
    assert not show_diff(lines, """\
  191.5439       0.00       0.06
   63.8485       1.97       1.39
   38.3098       1.96       1.35
   27.3648       1.59       1.41
   21.2845       1.52       1.41
   17.4154       1.84       2.86
   14.7368       1.86       1.50
   12.7727       1.88       1.47
   11.2708       1.91       2.04
   10.0852       1.85       1.39
""")

  assert os.path.exists("best.hkl")
  with open("best.hkl", "rb") as f:
    lines = ''.join(f.readlines()[:10])
    assert not show_diff(lines, """\
 -36   -4   17      -2.32      12.41
 -36   -4   18      -6.05      12.34
 -36   -3   16      11.93      14.98
 -36   -2   15      -2.02      15.18
 -36   -1   13      -7.12      15.63
 -36    0   11     -22.86      15.35
 -36    0   12      -3.67      13.66
 -36    1    9      23.31      15.71
 -36    1   10      -2.44      14.55
 -35   -4   18      30.60      16.10
""")

  assert os.path.exists("best.par")
  with open("best.par", "rb") as f:
    lines = f.read()
    assert not show_diff(lines, """\
# parameter file for BEST
TITLE          From DIALS
DETECTOR       PILA
SITE           Not set
DIAMETER       434.64
PIXEL          0.172
ROTAXIS        -0.00 0.00 1.00 FAST
POLAXIS        0.00 1.00 0.00
GAIN               1.00
CMOSAIC            0.23
PHISTART           0.00
PHIWIDTH           0.20
DISTANCE         191.01
WAVELENGTH      0.97950
POLARISATION    0.99900
SYMMETRY       P422
UB                  0.01      0.02      0.00
                    0.01      0.00     -0.02
                   -0.02      0.01     -0.00
CELL              42.19    42.19    39.68  90.00  90.00  90.00
RASTER           13  13   7   3   4
SEPARATION      2.960  2.960
BEAM            219.864  212.610
# end of parameter file for BEST
""")
  return


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"
