from __future__ import annotations

from unittest import mock

import dials.util.phil
from dials.array_family import flex


@mock.patch("dials.util.phil.ExperimentListFactory")
def test(ExperimentListFactory, dials_data):
    # Only use these filenames for verification
    path = dials_data("centroid_test_data", pathlib=True)
    experiments_path = path / "experiments.json"
    reflections_path1 = path / "integrated.pickle"
    reflections_path2 = path / "integrated.refl"

    phil_scope = dials.util.phil.parse(
        f"""
    input {{
      reflections1 = {reflections_path1}
        .type = reflection_table
      reflections2 = {reflections_path2}
        .type = reflection_table
      experiments = {experiments_path}
        .type = experiment_list
    }}
  """
    )

    params = phil_scope.extract()

    # Check the right filenames were parsed
    assert params.input.reflections1.filename == str(reflections_path1)
    assert params.input.reflections2.filename == str(reflections_path2)
    assert params.input.experiments.filename == str(experiments_path)
    # Check that we got the expected objects back
    assert isinstance(params.input.experiments.data, mock.Mock)
    assert isinstance(params.input.reflections1.data, flex.reflection_table)
    assert isinstance(params.input.reflections2.data, flex.reflection_table)
    # Check we had the correct calls made
    assert ExperimentListFactory.from_json_file.call_args[0] == (str(experiments_path),)
