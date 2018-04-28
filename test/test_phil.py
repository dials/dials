from __future__ import absolute_import, division, print_function

import os

import dials.util.phil
import mock

# Modules the phil parser uses, that we want to mock
import dials.array_family.flex

@mock.patch("os.path.exists", mock.Mock(return_value=True))
@mock.patch("dials.array_family.flex")
@mock.patch("dxtbx.model.experiment_list.ExperimentListFactory")
@mock.patch("dxtbx.datablock.DataBlockFactory")
def test(DataBlockFactory, ExperimentListFactory, flex):
  # Only use these filenames for verification
  path = "centroid_test_data"
  datablock_path = os.path.join(path, "datablock.json")
  experiments_path = os.path.join(path, "experiments.json")
  reflections_path = os.path.join(path, "integrated.pickle")


  phil_scope = dials.util.phil.parse('''
    input {
      reflections = %s
        .type = reflection_table
      datablock = %s
        .type = datablock
      experiments = %s
        .type = experiment_list
    }
  ''' % (reflections_path,
         datablock_path,
         experiments_path))

  params = phil_scope.extract()
  # Check the right filenames were parsed
  assert(params.input.reflections.filename == reflections_path)
  assert(params.input.datablock.filename == datablock_path)
  assert(params.input.experiments.filename == experiments_path)
  # Check that we got the expected objects back
  assert isinstance(params.input.reflections.data, mock.Mock)
  assert isinstance(params.input.datablock.data, mock.Mock)
  assert isinstance(params.input.experiments.data, mock.Mock)
  # Check we had the correct calls made
  flex.reflection_table.from_pickle.assert_called_once_with(reflections_path)
  assert DataBlockFactory.from_json_file.call_args[0] == (datablock_path,)
  assert ExperimentListFactory.from_json_file.call_args[0] == (experiments_path,)
