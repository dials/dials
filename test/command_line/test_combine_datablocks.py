from __future__ import absolute_import, division, print_function

import os

from dials.array_family import flex
from dxtbx.datablock import DataBlockFactory
from libtbx import easy_run

def test_combination_of_multiple_datablocks_and_strong_spots_files(dials_regression, tmpdir):
  tmpdir.chdir()

  path = os.path.join(dials_regression, "centroid_test_data/centroid_####.cbf")

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

  # check the reflections are unaffected, except for the change in id
  s1 = comb_strong.select(comb_strong['id'] == 0)
  s2 = comb_strong.select(comb_strong['id'] == 1)
  s2['id'] = flex.size_t(len(s2), 0)
  for r1, r2 in zip(s1, strong1):
    assert r1 == r2
  for r1, r2 in zip(s2, strong2):
    assert r1 == r2
