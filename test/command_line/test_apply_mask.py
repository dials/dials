from __future__ import absolute_import, division, print_function

import os
import procrunner

def test_datablocks(dials_regression, run_in_tmpdir):
  input_filename = os.path.join(dials_regression, "centroid_test_data", "datablock.json")
  mask_filename = os.path.join(dials_regression, "centroid_test_data", "lookup_mask.pickle")
  output_filename = "output_datablock.json"

  result = procrunner.run_process([
      'dials.apply_mask',
      'input.datablock=%s' % input_filename,
      'input.mask=%s' % mask_filename,
      'output.datablock=%s' % output_filename,
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  from dials.array_family import flex  # import dependency
  from dxtbx.datablock import DataBlockFactory
  datablocks = DataBlockFactory.from_json_file(output_filename)

  assert len(datablocks) == 1
  imagesets = datablocks[0].extract_imagesets()
  assert len(imagesets) == 1
  imageset = imagesets[0]
  assert imageset.external_lookup.mask.filename == mask_filename

def test_experimentss(dials_regression, run_in_tmpdir):
  input_filename = os.path.join(dials_regression, "centroid_test_data", "experiments.json")
  mask_filename = os.path.join(dials_regression, "centroid_test_data", "lookup_mask.pickle")
  output_filename = "output_experiments.json"

  result = procrunner.run_process([
      'dials.apply_mask',
      'input.experiments=%s' % input_filename,
      'input.mask=%s' % mask_filename,
      'output.experiments=%s' % output_filename,
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  from dials.array_family import flex  # import dependency
  from dxtbx.model.experiment_list import ExperimentListFactory
  experiments = ExperimentListFactory.from_json_file(output_filename)

  assert len(experiments) == 1
  imageset = experiments[0].imageset
  assert imageset.external_lookup.mask.filename == mask_filename
