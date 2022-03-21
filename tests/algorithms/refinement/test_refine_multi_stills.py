from __future__ import annotations

import os
from pathlib import Path

import procrunner
import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx import phil
from scitbx import matrix

from dials.algorithms.refinement import RefinerFactory
from dials.array_family import flex


def test1(dials_regression, tmp_path):
    data_dir = Path(dials_regression) / "refinement_test_data" / "multi_stills"

    result = procrunner.run(
        [
            "dials.refine",
            data_dir / "combined_experiments.json",
            data_dir / "combined_reflections.pickle",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    # load results
    reg_exp = ExperimentListFactory.from_json_file(
        data_dir / "regression_experiments.json", check_format=False
    )
    ref_exp = ExperimentListFactory.from_json_file(
        tmp_path / "refined.expt", check_format=False
    )

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


@pytest.mark.skipif(
    os.name == "nt",
    reason="Multiprocessing error on Windows: 'This class cannot be instantiated from Python'",
)
def test_multi_process_refinement_gives_same_results_as_single_process_refinement(
    dials_regression, tmp_path
):
    data_dir = Path(dials_regression) / "refinement_test_data" / "multi_stills"
    cmd = [
        "dials.refine",
        data_dir / "combined_experiments.json",
        data_dir / "combined_reflections.pickle",
        "outlier.algorithm=null",
        "engine=LBFGScurvs",
        "output.reflections=None",
    ]
    result = procrunner.run(
        cmd + ["output.experiments=refined_nproc4.expt", "nproc=4"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    result = procrunner.run(
        cmd + ["output.experiments=refined_nproc1.expt", "nproc=1"],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    # load results
    nproc1 = ExperimentListFactory.from_json_file(
        tmp_path / "refined_nproc1.expt", check_format=False
    )
    nproc4 = ExperimentListFactory.from_json_file(
        tmp_path / "refined_nproc4.expt", check_format=False
    )

    # compare results
    for b1, b2 in zip(nproc1.beams(), nproc4.beams()):
        assert b1.is_similar_to(b2)
    for c1, c2 in zip(nproc1.crystals(), nproc4.crystals()):
        assert c1.is_similar_to(c2)
    for d1, d2 in zip(nproc1.detectors(), nproc4.detectors()):
        assert d1.is_similar_to(
            d2,
            fast_axis_tolerance=5e-5,
            slow_axis_tolerance=5e-5,
            origin_tolerance=5e-5,
        )


def test_restrained_refinement_with_fixed_parameterisations(dials_regression, tmp_path):
    # Avoid a regression to https://github.com/dials/dials/issues/1142 by
    # testing that refinement succeeds when some parameterisations are fixed
    # by parameter auto reduction code, but restraints are requested for
    # those parameterisations.

    # The phil scope
    from dials.algorithms.refinement.refiner import phil_scope

    user_phil = phil.parse(
        """
refinement {
  parameterisation {
    auto_reduction {
      min_nref_per_parameter = 90
      action = fail *fix remove
    }
    crystal {
      unit_cell {
        restraints {
          tie_to_target {
            values = 95 95 132 90 90 120
            sigmas = 1 1 1 0 0 0
            id = 0 1 2 3 4 5 6 7 8 9
          }
        }
      }
    }
  }
}
"""
    )

    working_phil = phil_scope.fetch(source=user_phil)
    working_params = working_phil.extract()

    # use the multi stills test data
    data_dir = Path(dials_regression) / "refinement_test_data" / "multi_stills"
    experiments_path = data_dir / "combined_experiments.json"
    pickle_path = data_dir / "combined_reflections.pickle"

    experiments = ExperimentListFactory.from_json_file(
        experiments_path, check_format=False
    )
    reflections = flex.reflection_table.from_file(pickle_path)

    refiner = RefinerFactory.from_parameters_data_experiments(
        working_params, reflections, experiments
    )

    history = refiner.run()
    rmsd_limits = (0.2044, 0.2220, 0.0063)
    for a, b in zip(history["rmsd"][-1], rmsd_limits):
        assert a < b
