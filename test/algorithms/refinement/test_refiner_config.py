"""Test Refiners can be constructed with various configurations"""

from __future__ import absolute_import, division, print_function

import os
import pytest
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.refinement.refiner import phil_scope
from dials.algorithms.refinement import RefinerFactory
from dials.array_family import flex


@pytest.mark.parametrize(
    "detector_parameterisation_choice",
    [
        "automatic",
        pytest.param("single", marks=pytest.mark.xfail()),
        "multiple",
        "hierarchical",
    ],
)
def test_multi_panel_parameterisations(
    dials_regression, detector_parameterisation_choice
):

    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "cspad_refinement"
    )
    exp_file = os.path.join(data_dir, "cspad_refined_experiments_step6_level2_300.json")
    ref_file = os.path.join(data_dir, "cspad_reflections_step7_300.pickle")

    reflections = flex.reflection_table.from_pickle(ref_file)
    experiments = ExperimentListFactory.from_json_file(exp_file, check_format=False)

    # Set refinement parameters
    params = phil_scope.fetch(source=phil.parse("")).extract()
    params.refinement.parameterisation.detector.panels = (
        detector_parameterisation_choice
    )

    # Construct refiner
    refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments
    )
    assert refiner.experiment_type == "stills"
