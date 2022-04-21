from __future__ import annotations

from pathlib import Path

import procrunner
import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx import phil

from dials.algorithms.refinement import RefinerFactory
from dials.algorithms.refinement.engine import Journal
from dials.algorithms.refinement.refiner import phil_scope as refiner_phil_scope
from dials.array_family import flex


def test_joint_refinement(dials_regression, tmp_path):
    """A basic test of joint refinement of the CS-PAD detector at hierarchy level 2
    with 300 crystals."""

    bevington = pytest.importorskip("scitbx.examples.bevington")
    if not hasattr(bevington, "non_linear_ls_eigen_wrapper"):
        pytest.skip("Skipping test as SparseLevMar engine not available")

    data_dir = Path(dials_regression) / "refinement_test_data" / "xfel_metrology"

    # Do refinement and load the history
    result = procrunner.run(
        [
            "dials.refine",
            data_dir / "benchmark_level2d.json",
            data_dir / "benchmark_level2d.pickle",
            data_dir / "refine.phil",
            "history=history.json",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    # there are plenty of things we could do with the refinement history, but
    # here just check that final RMSDs are low enough
    history = Journal.from_json_file(tmp_path / "history.json")
    final_rmsd = history["rmsd"][-1]
    assert final_rmsd[0] < 0.0354
    assert final_rmsd[1] < 0.0406
    assert final_rmsd[2] < 0.0018

    # also check that the used_in_refinement flag got set correctly
    rt = flex.reflection_table.from_file(tmp_path / "refined.refl")
    uir = rt.get_flags(rt.flags.used_in_refinement)
    assert uir.count(True) == history["num_reflections"][-1]


def test_constrained_refinement(dials_regression, tmp_path):
    """Do constrained refinement, checking that a panel group with no data
    on it still moves with its partners in the constraint.
    See https://github.com/dials/dials/issues/990"""

    bevington = pytest.importorskip("scitbx.examples.bevington")
    if not hasattr(bevington, "non_linear_ls_eigen_wrapper"):
        pytest.skip("Skipping test as SparseLevMar engine not available")

    data_dir = Path(dials_regression) / "refinement_test_data" / "xfel_metrology"

    # Load experiments and reflections
    refl = flex.reflection_table.from_file(data_dir / "benchmark_level2d.pickle")
    expt = ExperimentListFactory.from_json_file(data_dir / "benchmark_level2d.json")

    # There are zero reflections on some panels, so these will only move via constraints
    for i in [8, 10, 11, 26, 27, 40, 42, 43, 56, 58, 59]:
        assert (refl["panel"] == i).count(True) == 0

    # Get parameters, combining refine.phil with constraints that enforce distances to move in lockstep
    refine_phil = phil.parse(data_dir.joinpath("refine.phil").read_text())
    constraint_phil = phil.parse(
        """
refinement {
  parameterisation {
    detector {
      fix_list=Tau2,Tau3
      constraints {
        parameter=Dist
      }
    }
  }
}
"""
    )
    params = refiner_phil_scope.fetch(sources=[refine_phil, constraint_phil]).extract()

    detector = expt.detectors()[0]
    initial_distances = [p.get_distance() for p in detector]

    # Set up a refiner
    refiner = RefinerFactory.from_parameters_data_experiments(params, refl, expt)

    refiner.run()
    detector = refiner.get_experiments().detectors()[0]

    final_distances = [p.get_distance() for p in detector]

    # The shifts between initial and final distances should all be equal
    dist_diff = [a - b for a, b in zip(final_distances, initial_distances)]
    for d in dist_diff[1:]:
        assert d == pytest.approx(dist_diff[0])
