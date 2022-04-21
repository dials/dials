"""Tests for ΔCC½ algorithms."""

from __future__ import annotations

from unittest import mock

from dxtbx.model import Crystal, Experiment, ExperimentList, Scan

from dials.algorithms.statistics.cc_half_algorithm import CCHalfFromDials
from dials.array_family import flex
from dials.command_line.compute_delta_cchalf import phil_scope


def generated_exp(n=1):
    """Generate an experiment list with two experiments."""
    experiments = ExperimentList()
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [1.0, 0.0, 0.0],
        "real_space_b": [0.0, 1.0, 0.0],
        "real_space_c": [0.0, 0.0, 2.0],
        "space_group_hall_symbol": " C 2y",
    }
    for i in range(n):
        experiments.append(
            Experiment(
                scan=Scan(image_range=[1, 25], oscillation=[0.0, 1.0]),
                crystal=Crystal.from_dict(exp_dict),
                identifier=str(i),
            )
        )
    return experiments


def generated_refl():
    """Generate test data."""
    refls = flex.reflection_table()
    refls["intensity.scale.value"] = flex.double(range(50))
    refls["intensity.scale.variance"] = flex.double(range(50))
    refls["inverse_scale_factor"] = flex.double(50, 1.0)
    refls["id"] = flex.int([0] * 25 + [1] * 25)
    refls["xyzobs.px.value"] = flex.vec3_double(
        [(0, 0, i + 0.5) for i in range(25)] * 2
    )
    vals = [0, 1, 2, 3, 48, 49]
    # set first two of first sweep, and last two of last sweep as outliers
    outliers = flex.bool(50, False)
    for v in vals:
        outliers[v] = True
    refls.set_flags(outliers, refls.flags.outlier_in_scaling)
    refls.experiment_identifiers()[0] = "0"
    refls.experiment_identifiers()[1] = "1"
    return refls


def test_setup_of_CCHalfFromDials():
    """Test the correct setup in image group mode.

    Test for the case of outliers at the end of images, and image ranges not
    equaling a multiple of the grouping."""

    params = phil_scope.extract()
    params.mode = "image_group"
    expts = generated_exp(n=2)
    refls = generated_refl()

    # Expected behaviour is that the outliers will not be included in the
    # image range, and that all groups will have at least 10 images in.
    script = CCHalfFromDials(params, expts, refls)
    assert script.group_to_datasetid_and_range == {
        0: ("0", (5, 14)),
        1: ("0", (15, 25)),
        2: ("1", (1, 10)),
        3: ("1", (11, 23)),
    }
    assert script.datasetid_to_groups == {"0": [0, 1], "1": [2, 3]}


def test_exclusion_in_CCHalfFromDials():
    """Test the exclusion of image groups."""
    # Same input as above, but mock DeltaCCHalf algorithm to just test
    # interpretation of results and setting of excluded regions. With the
    # input, test that outlier edges are correctly removed.

    params = phil_scope.extract()
    params.mode = "image_group"
    expts = generated_exp(n=2)
    refls = generated_refl()

    def mock_algorithm(*_):
        """Mock a result from DeltaCCHalf"""
        algo = mock.Mock()
        algo.run.return_value = None
        algo.results_summary = {
            "per_dataset_delta_cc_half_values": {
                "delta_cc_half_values": [-5.0, -2.0, 4.0, -5.0],
                "datasets": [0, 1, 2, 3],
            },
            "dataset_removal": {"cutoff_value": -1.0},
        }
        return algo

    with mock.patch(
        "dials.algorithms.statistics.cc_half_algorithm.DeltaCCHalf",
        side_effect=mock_algorithm,
    ):
        script = CCHalfFromDials(params, expts, refls)
        script.run()

        assert script.datasetid_to_groups == {"0": [], "1": [2]}  # all but 3 removed
        expts = script.experiments

        assert list(expts.identifiers()) == ["1"]
        assert expts[0].scan.get_valid_image_ranges(expts.identifiers()[0]) == [(1, 10)]
        assert script.results_summary["dataset_removal"][
            "experiment_ids_fully_removed"
        ] == [0]
        assert script.results_summary["dataset_removal"][
            "experiments_fully_removed"
        ] == ["0"]
