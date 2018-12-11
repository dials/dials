from __future__ import absolute_import, division, print_function

import os

import dials.util.phil
import mock

# Modules the phil parser uses, that we want to mock
from dials.array_family import flex

@mock.patch("os.path.exists", mock.Mock(return_value=True))
@mock.patch("dxtbx.model.experiment_list.ExperimentListFactory")
@mock.patch("dxtbx.datablock.DataBlockFactory")
def test(DataBlockFactory, ExperimentListFactory, dials_regression):
  # Only use these filenames for verification
  path = "centroid_test_data"
  datablock_path = os.path.join(dials_regression, path, "datablock.json")
  experiments_path = os.path.join(dials_regression, path, "experiments.json")
  reflections_path1 = os.path.join(dials_regression, path, "integrated.pickle")
  reflections_path2 = os.path.join(dials_regression, path, "integrated.mpack")

  phil_scope = dials.util.phil.parse('''
    input {
      reflections1 = %s
        .type = reflection_table
      reflections2 = %s
        .type = reflection_table
      datablock = %s
        .type = datablock
      experiments = %s
        .type = experiment_list
    }
  ''' % (reflections_path1,
         reflections_path2,
         datablock_path,
         experiments_path))

  params = phil_scope.extract()

  # Check the right filenames were parsed
  assert(params.input.reflections1.filename == reflections_path1)
  assert(params.input.reflections2.filename == reflections_path2)
  assert(params.input.datablock.filename == datablock_path)
  assert(params.input.experiments.filename == experiments_path)
  # Check that we got the expected objects back
  assert isinstance(params.input.datablock.data, mock.Mock)
  assert isinstance(params.input.experiments.data, mock.Mock)
  # Check we had the correct calls made
  assert DataBlockFactory.from_json_file.call_args[0] == (datablock_path,)
  assert ExperimentListFactory.from_json_file.call_args[0] == (experiments_path,)
