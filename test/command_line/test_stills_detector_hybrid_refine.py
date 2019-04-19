"""
Test dials.stills_detector_hybrid_refine by running a short job
"""

from __future__ import absolute_import, division, print_function

import os
import procrunner
import pytest
from dials.command_line.stills_detector_hybrid_refine import run as run_hybrid_refine


@pytest.mark.parametrize("averaged_reference_detector", [True, False])
def test(dials_regression, run_in_tmpdir, averaged_reference_detector):
    # use 20 indexed pickles from CXI for this test
    data_dir = os.path.join(
        dials_regression, "stills_test_data", "cspad_indexing_results", "cxid9114_r0097"
    )

    experiments_paths = [
        "idx-20140615231407812_refined_experiments.json",
        "idx-20140615231408153_refined_experiments.json",
        "idx-20140615231408912_refined_experiments.json",
        "idx-20140615231409020_refined_experiments.json",
        "idx-20140615231409470_refined_experiments.json",
        "idx-20140615231410120_refined_experiments.json",
        "idx-20140615231411153_refined_experiments.json",
        "idx-20140615231411170_refined_experiments.json",
        "idx-20140615231411970_refined_experiments.json",
        "idx-20140615231412103_refined_experiments.json",
        "idx-20140615231412495_refined_experiments.json",
        "idx-20140615231413370_refined_experiments.json",
        "idx-20140615231413878_refined_experiments.json",
        "idx-20140615231414128_refined_experiments.json",
        "idx-20140615231414461_refined_experiments.json",
        "idx-20140615231415353_refined_experiments.json",
        "idx-20140615231415761_refined_experiments.json",
        "idx-20140615231415819_refined_experiments.json",
        "idx-20140615231416328_refined_experiments.json",
        "idx-20140615231416694_refined_experiments.json",
    ]

    reflections_paths = [
        "idx-20140615231407812_indexed.pickle",
        "idx-20140615231408153_indexed.pickle",
        "idx-20140615231408912_indexed.pickle",
        "idx-20140615231409020_indexed.pickle",
        "idx-20140615231409470_indexed.pickle",
        "idx-20140615231410120_indexed.pickle",
        "idx-20140615231411153_indexed.pickle",
        "idx-20140615231411170_indexed.pickle",
        "idx-20140615231411970_indexed.pickle",
        "idx-20140615231412103_indexed.pickle",
        "idx-20140615231412495_indexed.pickle",
        "idx-20140615231413370_indexed.pickle",
        "idx-20140615231413878_indexed.pickle",
        "idx-20140615231414128_indexed.pickle",
        "idx-20140615231414461_indexed.pickle",
        "idx-20140615231415353_indexed.pickle",
        "idx-20140615231415761_indexed.pickle",
        "idx-20140615231415819_indexed.pickle",
        "idx-20140615231416328_indexed.pickle",
        "idx-20140615231416694_indexed.pickle",
    ]

    cmd = []
    for exp in experiments_paths:
        exp = os.path.join(data_dir, exp)
        cmd.append("experiments={0}".format(exp))
    for ref in reflections_paths:
        ref = os.path.join(data_dir, ref)
        cmd.append("reflections={0}".format(ref))

    if averaged_reference_detector:
        # specify hierarchy_level=0
        cmd.append(
            "detector_phase.refinement.parameterisation.detector.hierarchy_level=0"
        )
        cmd.append("reference_detector=average")
    else:
        # specify hierarchy_level=1
        cmd.append(
            "detector_phase.refinement.parameterisation.detector.hierarchy_level=1"
        )
    result = run_hybrid_refine(args=cmd)
    # TODO: Test used to assert empty stderr. Should do this again once
    #       we have access to contextlib.redirect_stderr
    assert not result
