from __future__ import annotations

import pathlib

import procrunner

from dxtbx.model.experiment_list import ExperimentListFactory


def test_experiments(dials_data, tmp_path):
    input_filename = (
        dials_data("centroid_test_data", pathlib=True) / "imported_experiments.json"
    )
    mask_filename = (
        dials_data("centroid_test_data", pathlib=True) / "lookup_mask.pickle"
    )
    output_filename = tmp_path / "output.expt"

    result = procrunner.run(
        [
            "dials.apply_mask",
            f"input.experiments={input_filename}",
            f"input.mask={mask_filename}",
            f"output.experiments={output_filename}",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    experiments = ExperimentListFactory.from_json_file(output_filename)

    assert len(experiments) == 1
    imageset = experiments[0].imageset
    assert pathlib.Path(imageset.external_lookup.mask.filename) == mask_filename
