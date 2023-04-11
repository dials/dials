from __future__ import annotations

import procrunner

from dxtbx.model.experiment_list import ExperimentListFactory


def test_sequence_to_stills(dials_data, tmpdir):
    data_dir = dials_data("insulin_processed", pathlib=True)
    input_experiments = data_dir / "integrated.expt"
    input_reflections = data_dir / "refined.refl"
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
