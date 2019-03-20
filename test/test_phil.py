from __future__ import absolute_import, division, print_function

from dials.array_family import flex
import dials.util.phil
import mock


@mock.patch("dxtbx.model.experiment_list.ExperimentListFactory")
def test(ExperimentListFactory, dials_data):
    # Only use these filenames for verification
    path = dials_data("centroid_test_data")
    experiments_path = path.join("experiments.json").strpath
    reflections_path1 = path.join("integrated.pickle").strpath
    reflections_path2 = path.join("integrated.mpack").strpath

    phil_scope = dials.util.phil.parse(
        """
    input {
      reflections1 = %s
        .type = reflection_table
      reflections2 = %s
        .type = reflection_table
      experiments = %s
        .type = experiment_list
    }
  """
        % (reflections_path1, reflections_path2, experiments_path)
    )

    params = phil_scope.extract()

    # Check the right filenames were parsed
    assert params.input.reflections1.filename == reflections_path1
    assert params.input.reflections2.filename == reflections_path2
    assert params.input.experiments.filename == experiments_path
    # Check that we got the expected objects back
    assert isinstance(params.input.experiments.data, mock.Mock)
    assert isinstance(params.input.reflections1.data, flex.reflection_table)
    assert isinstance(params.input.reflections2.data, flex.reflection_table)
    # Check we had the correct calls made
    assert ExperimentListFactory.from_json_file.call_args[0] == (experiments_path,)
