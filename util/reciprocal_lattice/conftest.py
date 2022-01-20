import os

import pytest

from dxtbx.serialize import load

from dials.array_family import flex


@pytest.fixture
def multi_sequence_data(dials_regression):
    experiments = load.experiment_list(
        os.path.join(
            dials_regression, "indexing_test_data", "multi_sweep", "experiments.json"
        ),
        check_format=False,
    )
    reflections = flex.reflection_table.from_file(
        os.path.join(
            dials_regression, "indexing_test_data", "multi_sweep", "indexed.pickle"
        )
    )
    return {"reflections": reflections, "experiments": experiments}


@pytest.fixture
def centroid_test_data(dials_data):
    experiments = load.experiment_list(
        dials_data("centroid_test_data").join("experiments.json").strpath,
        check_format=False,
    )
    reflections = flex.reflection_table.from_file(
        dials_data("centroid_test_data").join("integrated.pickle").strpath
    )
    return {"reflections": reflections, "experiments": experiments}
