from __future__ import absolute_import, division, print_function

import os

import dials.util.phil
import mock

@mock.patch("os.path.exists", mock.Mock(return_value=True))
@mock.patch("dxtbx.model.experiment_list.ExperimentListFactory")
def test(ExperimentListFactory, dials_regression):
  # Only use these filenames for verification
  path = "centroid_test_data"
  experiments_path = os.path.join(dials_regression, path, "experiments.json")
  reflections_path1 = os.path.join(dials_regression, path, "integrated.pickle")
  reflections_path2 = os.path.join(dials_regression, path, "integrated.mpack")

  phil_scope = dials.util.phil.parse('''
    input {
      reflections1 = %s
        .type = reflection_table
      reflections2 = %s
        .type = reflection_table
      experiments = %s
        .type = experiment_list
    }
  ''' % (reflections_path1,
         reflections_path2,
         experiments_path))

  params = phil_scope.extract()

  # Check the right filenames were parsed
  assert(params.input.reflections1.filename == reflections_path1)
  assert(params.input.reflections2.filename == reflections_path2)
  assert(params.input.experiments.filename == experiments_path)
  # Check that we got the expected objects back
  assert isinstance(params.input.experiments.data, mock.Mock)
  # Check we had the correct calls made
  assert ExperimentListFactory.from_json_file.call_args[0] == (experiments_path,)
