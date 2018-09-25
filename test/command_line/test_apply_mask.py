from __future__ import absolute_import, division, print_function

import os
import procrunner

def test(dials_regression, run_in_tmpdir):
  input_filename = os.path.join(dials_regression, "centroid_test_data",
                                "datablock.json")
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
  imagesets = experiments.imagesets()
  assert len(imagesets) == 1
  imageset = imagesets[0]
  assert imageset.external_lookup.mask.filename == mask_filename
