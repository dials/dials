from __future__ import annotations

""" Test the situation that led to https://github.com/dials/dials/issues/423.
In that case instantiating a Refiner for an experiment list with an I23
detector model caused the panel origins to move before any refinement took
place. This occurred because for the input models.expt the root frame for
the hierarchical detector is on source side of the laboratory frame origin, not
on the detector side. Prior to the fix this resulted in incorrect calculation
of the offsets of all panels from the root frame.
"""


def test_run(dials_data):
    from dxtbx.model.experiment_list import ExperimentListFactory
    from libtbx import phil

    from dials.algorithms.refinement import RefinerFactory
    from dials.algorithms.refinement.refiner import phil_scope
    from dials.array_family import flex

    data_dir = dials_data("refinement_test_data", pathlib=True)
    exp_file = data_dir / "dials-423.json"
    ref_file = data_dir / "dials-423.pickle"

    reflections = flex.reflection_table.from_file(ref_file)
    experiments = ExperimentListFactory.from_json_file(exp_file, check_format=False)
    """Test that the detector remains similar after refiner construction"""

    params = phil_scope.fetch(source=phil.parse("")).extract()

    # disable outlier rejection for speed of refiner construction
    params.refinement.reflections.outlier.algorithm = "null"

    refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments
    )

    d1 = experiments[0].detector
    d2 = refiner.get_experiments()[0].detector

    assert d1.is_similar_to(d2)
