# Test multiple stills refinement.

from __future__ import absolute_import, division, print_function

import os

from dxtbx.model.experiment_list import ExperimentListFactory
import procrunner


def test1(dials_regression, run_in_tmpdir):
    """
    Refinement test of 300 CSPAD images, testing auto_reduction, parameter
    fixing, constraints, SparseLevMar, and sauter_poon outlier rejection. See
    README in the regression folder for more details.
    """
    from scitbx import matrix

    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "cspad_refinement"
    )

    result = procrunner.run(
        [
            "dials.refine",
            os.path.join(data_dir, "cspad_refined_experiments_step6_level2_300.json"),
            os.path.join(data_dir, "cspad_reflections_step7_300.pickle"),
            os.path.join(data_dir, "refine.phil"),
        ]
    )
    assert not result.returncode and not result.stderr

    # load results
    reg_exp = ExperimentListFactory.from_json_file(
        os.path.join(data_dir, "regression_experiments.json"), check_format=False
    )
    ref_exp = ExperimentListFactory.from_json_file("refined.expt", check_format=False)

    # compare results
    tol = 1e-5
    for b1, b2 in zip(reg_exp.beams(), ref_exp.beams()):
        assert b1.is_similar_to(
            b2,
            wavelength_tolerance=tol,
            direction_tolerance=tol,
            polarization_normal_tolerance=tol,
            polarization_fraction_tolerance=tol,
        )
        s0_1 = matrix.col(b1.get_unit_s0())
        s0_2 = matrix.col(b2.get_unit_s0())
        assert s0_1.accute_angle(s0_2, deg=True) < 0.0057  # ~0.1 mrad
    for c1, c2 in zip(reg_exp.crystals(), ref_exp.crystals()):
        assert c1.is_similar_to(c2)

    for d1, d2 in zip(reg_exp.detectors(), ref_exp.detectors()):
        assert d1.is_similar_to(
            d2,
            fast_axis_tolerance=1e-4,
            slow_axis_tolerance=1e-4,
            origin_tolerance=1e-2,
        )
