"""Test combination of multiple datablocks and strong spots files."""

# python imports
from __future__ import absolute_import, division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from dxtbx.datablock import DataBlockFactory
from dials.array_family import flex

def test1():

  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(
    dials_regression, "centroid_test_data/centroid_####.cbf")

  # example combined two different spot-finding settings for the same dataset
  # e.d. for comparison with the reciprocal lattice viewer.
  cmd = "dials.import template={0}".format(path)
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = "dials.find_spots datablock.json output.reflections=strong1.pickle"
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = ("dials.find_spots datablock.json sigma_strong=5 "
         "output.reflections=strong2.pickle")
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  cmd = ("dev.dials.combine_datablocks datablock.json datablock.json "
         "strong1.pickle strong2.pickle")
  result = easy_run.fully_buffered(cmd).raise_if_errors()

  # load results
  comb_db = DataBlockFactory.from_json_file('combined_datablocks.json')[0]
  comb_strong = flex.reflection_table.from_pickle("combined_strong.pickle")

  # load reference models and reflections
  db = DataBlockFactory.from_json_file('datablock.json')[0]
  ref_detector = db.unique_detectors()[0]
  ref_beam = db.unique_beams()[0]
  ref_scan = db.unique_scans()[0]
  ref_goniometer = db.unique_goniometers()[0]
  strong1 = flex.reflection_table.from_pickle("strong1.pickle")
  strong2 = flex.reflection_table.from_pickle("strong2.pickle")

  # check the models have not been modified
  for imset in comb_db.extract_imagesets():
    assert imset.get_detector() == ref_detector
    assert imset.get_beam() == ref_beam
    assert imset.get_scan() == ref_scan
    assert imset.get_goniometer() == ref_goniometer
  print "OK"

  # check the reflections are unaffected, except for the change in id
  s1 = comb_strong.select(comb_strong['id'] == 0)
  s2 = comb_strong.select(comb_strong['id'] == 1)
  s2['id'] = flex.size_t(len(s2), 0)
  for r1, r2 in zip(s1, strong1):
    assert r1 == r2
  print "OK"
  for r1, r2 in zip(s2, strong2):
    assert r1 == r2
  print "OK"

  return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  test1()

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    from libtbx.utils import show_times_at_exit
    show_times_at_exit()
    run()
