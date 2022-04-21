from __future__ import annotations

import os

import procrunner

from dxtbx.model.experiment_list import ExperimentListFactory


def test_sequence_to_stills(dials_regression, tmpdir):
    path = os.path.join(
        dials_regression, "refinement_test_data", "radiation_damaged_thaumatin"
    )
    input_experiments = os.path.join(path, "refined_experiments_P42.json")
    input_reflections = os.path.join(path, "indexed.pickle")
    result = procrunner.run(
        [
            "dials.sequence_to_stills",
            input_experiments,
            input_reflections,
            "domain_size_ang=500",
            "half_mosaicity_deg=0.1",
            "max_scan_points=10",
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr

    assert tmpdir.join("stills.expt").check(file=1)
    assert tmpdir.join("stills.refl").check(file=1)

    experiments = ExperimentListFactory.from_json_file(
        tmpdir.join("stills.expt").strpath, check_format=False
    )
    assert len(experiments) == 10
