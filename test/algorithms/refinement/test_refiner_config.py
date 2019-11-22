"""Test Refiners can be constructed with various configurations"""

from __future__ import absolute_import, division, print_function

import os

import pytest
from dials.algorithms.refinement.refiner import phil_scope
from dials.algorithms.refinement import DialsRefineConfigError, RefinerFactory
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx import phil
from dials.util.slice import slice_reflections
from dials.algorithms.refinement.refiner import _trim_scans_to_observations
from copy import deepcopy


@pytest.mark.parametrize(
    "detector_parameterisation_choice",
    ["automatic", "single", "multiple", "hierarchical"],
)
def test_multi_panel_parameterisations(
    dials_regression, detector_parameterisation_choice
):

    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "cspad_refinement"
    )
    exp_file = os.path.join(data_dir, "cspad_refined_experiments_step6_level2_300.json")
    ref_file = os.path.join(data_dir, "cspad_reflections_step7_300.pickle")

    reflections = flex.reflection_table.from_file(ref_file)
    experiments = ExperimentListFactory.from_json_file(exp_file, check_format=False)

    # Set refinement parameters
    params = phil_scope.fetch(source=phil.parse("")).extract()
    params.refinement.parameterisation.detector.panels = (
        detector_parameterisation_choice
    )

    # Construct refiner

    if detector_parameterisation_choice == "single":
        with pytest.raises(DialsRefineConfigError):
            # Cannot create a single panel parameterisation for a multi-panel detector
            RefinerFactory.from_parameters_data_experiments(
                params, reflections, experiments
            )
    else:
        refiner = RefinerFactory.from_parameters_data_experiments(
            params, reflections, experiments
        )
        assert refiner.experiment_type == "stills"


def test_trim_scans_to_observations(dials_data):

    # Use 4 scan data for this test
    data_dir = dials_data("l_cysteine_dials_output")
    experiments = ExperimentListFactory.from_json_file(
        (data_dir / "indexed.expt").strpath, check_format=False
    )
    reflections = flex.reflection_table.from_file((data_dir / "indexed.refl").strpath)

    # Check the image range is what we expect
    image_ranges = [e.scan.get_image_range() for e in experiments]
    for a, b in zip(image_ranges, [(1, 1700), (1, 1700), (1, 1700), (1, 1800)]):
        assert a == b

    # If image range unchanged, nothing should happen
    trim_expt = _trim_scans_to_observations(deepcopy(experiments), reflections)
    new_ranges = [e.scan.get_image_range() for e in trim_expt]
    for a, b in zip(image_ranges, new_ranges):
        assert a == b

    # Slice 20 images off head and tail
    sliced_ranges = [(r[0] + 20, r[1] - 20) for r in image_ranges]
    sliced = slice_reflections(reflections, sliced_ranges)

    # Now trimmed scans should have ranges array ranges equal to their min, max
    # shoebox z coords
    trim_expt = _trim_scans_to_observations(deepcopy(experiments), sliced)
    new_ranges = [e.scan.get_image_range() for e in trim_expt]

    for i, e in enumerate(trim_expt):
        refs = sliced.select(sliced["id"] == i)
        bb = refs["shoebox"].bounding_boxes()
        z_min, z_max = bb.parts()[4:]
        assert new_ranges[i] == (min(z_min), max(z_max))

    # Now delete shoebox data. Trimmed scans will be wider than the observed
    # range by >0.5 deg at each end
    del sliced["shoebox"]
    trim_expt = _trim_scans_to_observations(deepcopy(experiments), sliced)
    new_ranges = [e.scan.get_image_range() for e in trim_expt]
    for i, e in enumerate(trim_expt):
        refs = sliced.select(sliced["id"] == i)
        z = refs["xyzobs.px.value"].parts()[2]
        im_width = e.scan.get_oscillation()[1]
        assert ((min(z) - new_ranges[i][0]) / im_width) > 0.5
        assert ((new_ranges[i][1] - max(z)) / im_width) > 0.5
