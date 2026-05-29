from __future__ import annotations

import pathlib
import shutil
import subprocess

from dxtbx.model.experiment_list import ExperimentListFactory


def test_experiments(dials_data, tmp_path):
    input_filename = dials_data("centroid_test_data") / "imported_experiments.json"
    mask_filename = dials_data("centroid_test_data") / "lookup_mask.pickle"
    output_filename = tmp_path / "output.expt"

    result = subprocess.run(
        [
            shutil.which("dials.apply_mask"),
            f"input.experiments={input_filename}",
            f"input.mask={mask_filename}",
            f"output.experiments={output_filename}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    experiments = ExperimentListFactory.from_json_file(output_filename)

    assert len(experiments) == 1
    imageset = experiments[0].imageset
    assert pathlib.Path(imageset.external_lookup.mask.filename) == mask_filename
