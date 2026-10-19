from __future__ import annotations

from unittest import mock

import dials.util.phil
from dials.array_family import flex


@mock.patch("dials.util.phil.ExperimentListFactory")
def test(ExperimentListFactory, dials_data):
    # Only use these filenames for verification
    path = dials_data("centroid_test_data")
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


def test_files_are_read_lazily(dials_data):
    # libtbx.phil converts words to objects on every extract and fetch_diff, so
    # the converters must not read the file until the data are used
    reflections_path = dials_data("centroid_test_data") / "integrated.refl"
    master = dials.util.phil.parse(
        """
    input {
      reflections = None
        .type = reflection_table
    }
  """
    )
    user = dials.util.phil.parse(f"input.reflections = {reflections_path}")

    with mock.patch.object(
        flex.reflection_table, "from_file", wraps=flex.reflection_table.from_file
    ) as from_file:
        working = master.fetch(source=user)
        assert master.fetch_diff(source=working).as_str() != ""
        params = working.extract()
        assert params.input.reflections.filename == str(reflections_path)
        from_file.assert_not_called()

        data = params.input.reflections.data
        assert isinstance(data, flex.reflection_table)
        assert params.input.reflections.data is data
        from_file.assert_called_once_with(str(reflections_path))

    # It can still be used like the FilenameDataWrapper namedtuple
    filename, unpacked = params.input.reflections
    assert filename == str(reflections_path)
    assert unpacked is data
    assert params.input.reflections[1] is data
