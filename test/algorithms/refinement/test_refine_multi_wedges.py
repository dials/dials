"""
Test refinement of multiple narrow sweeps.
"""

from __future__ import absolute_import, division, print_function

import os

import procrunner
from scitbx import matrix
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.refinement.engine import Journal
from dials.array_family import flex  # noqa: F401
import pytest


def test(dials_regression, run_in_tmpdir):
    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "multi_narrow_wedges"
    )

    selection = (2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 17, 18, 19, 20)

    # Combine all the separate sweeps

    result = procrunner.run(
        [
            "dials.combine_experiments",
            "reference_from_experiment.beam=0",
            "reference_from_experiment.goniometer=0",
            "reference_from_experiment.detector=0",
        ]
        + [
            "experiments={}/data/sweep_%03d/experiments.json".format(data_dir) % n
            for n in selection
        ]
        + [
            "reflections={}/data/sweep_%03d/reflections.pickle".format(data_dir) % n
            for n in selection
        ]
    )
    assert not result.returncode and not result.stderr

    # Do refinement and load the results

    # turn off outlier rejection so that test takes about 4s rather than 10s
    # set close_to_spindle_cutoff to old default
    result = procrunner.run(
        [
            "dials.refine",
            "combined.expt",
            "combined.refl",
            "scan_varying=false",
            "outlier.algorithm=null",
            "close_to_spindle_cutoff=0.05",
        ]
    )
    assert not result.returncode and not result.stderr

    refined_experiments = ExperimentListFactory.from_json_file(
        "refined.expt", check_format=False
    )

    # Check results are as expected

    regression_experiments = ExperimentListFactory.from_json_file(
        os.path.join(data_dir, "regression_experiments.json"), check_format=False
    )

    for e1, e2 in zip(refined_experiments, regression_experiments):
        assert e1.crystal.is_similar_to(e2.crystal)
        # FIXME need is_similar_to for detector that checks geometry
        # assert e1.detector == e2.detector
        s0_1 = matrix.col(e1.beam.get_unit_s0())
        s0_2 = matrix.col(e1.beam.get_unit_s0())
        assert s0_1.accute_angle(s0_2, deg=True) < 0.0057  # ~0.1 mrad


def test_order_invariance(dials_regression, run_in_tmpdir):
    """Check that the order that datasets are included in refinement does not
    matter"""

    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "multi_narrow_wedges"
    )
    selection1 = (2, 3, 4, 5, 6)
    selection2 = (2, 3, 4, 6, 5)

    # First run
    result = procrunner.run(
        [
            "dials.combine_experiments",
            "reference_from_experiment.beam=0",
            "reference_from_experiment.goniometer=0",
            "reference_from_experiment.detector=0",
        ]
        + [
            "experiments={0}/data/sweep_%03d/experiments.json".format(data_dir) % n
            for n in selection1
        ]
        + [
            "reflections={0}/data/sweep_%03d/reflections.pickle".format(data_dir) % n
            for n in selection1
        ]
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        [
            "dials.refine",
            "combined.expt",
            "combined.refl",
            "scan_varying=false",
            "outlier.algorithm=tukey",
            "history=history1.json",
            "output.experiments=refined1.expt",
            "output.reflections=refined1.refl",
        ]
    )
    assert not result.returncode and not result.stderr

    # Second run
    result = procrunner.run(
        [
            "dials.combine_experiments",
            "reference_from_experiment.beam=0",
            "reference_from_experiment.goniometer=0",
            "reference_from_experiment.detector=0",
        ]
        + [
            "experiments={0}/data/sweep_%03d/experiments.json".format(data_dir) % n
            for n in selection2
        ]
        + [
            "reflections={0}/data/sweep_%03d/reflections.pickle".format(data_dir) % n
            for n in selection2
        ]
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        [
            "dials.refine",
            "combined.expt",
            "combined.refl",
            "scan_varying=false",
            "outlier.algorithm=tukey",
            "history=history2.json",
            "output.experiments=refined2.expt",
            "output.reflections=refined2.refl",
        ]
    )
    assert not result.returncode and not result.stderr

    # Load results
    refined_experiments1 = ExperimentListFactory.from_json_file(
        "refined1.expt", check_format=False
    )
    refined_experiments2 = ExperimentListFactory.from_json_file(
        "refined2.expt", check_format=False
    )

    history1 = Journal.from_json_file("history1.json")
    history2 = Journal.from_json_file("history1.json")

    # Compare RMSDs
    rmsd1 = history1["rmsd"]
    rmsd2 = history2["rmsd"]
    for a, b in zip(rmsd1, rmsd2):
        assert a == pytest.approx(b)

    # Compare crystals
    crystals1 = [exp.crystal for exp in refined_experiments1]
    crystals2 = [exp.crystal for exp in refined_experiments2[0:8]]
    crystals2.extend([exp.crystal for exp in refined_experiments2[13:16]])
    crystals2.extend([exp.crystal for exp in refined_experiments2[8:13]])
    for a, b in zip(crystals1, crystals2):
        assert a.is_similar_to(b)
