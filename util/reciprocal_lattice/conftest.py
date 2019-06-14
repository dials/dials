from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.serialize import load

from dials.array_family import flex


@pytest.fixture
def multi_sweep_data(dials_regression):
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
def centroid_test_data(dials_regression):
    experiments = load.experiment_list(
        os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
        check_format=False,
    )
    reflections = flex.reflection_table.from_file(
        os.path.join(dials_regression, "centroid_test_data", "integrated.pickle")
    )
    return {"reflections": reflections, "experiments": experiments}
