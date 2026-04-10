from __future__ import annotations

import json
import shutil
import subprocess
from os.path import join

from dials.array_family import flex


def test_tof_integration(dials_data, tmp_path):
    image_file = join(dials_data("isis_sxd_example_data"), "sxd_nacl_run.nxs")

    # Update .expt file with correct image path
    init_expt_path = join(dials_data("isis_sxd_nacl_processed"), "refined.expt")
    expt_path = join(tmp_path, "refined.expt")
    with open(init_expt_path, "rb") as g:
        expt_json = json.load(g)

    expt_json["imageset"][0]["template"] = image_file
    with open(expt_path, "w") as g:
        json.dump(expt_json, g)

    reflections = flex.reflection_table.from_msgpack_file(
        join(dials_data("isis_sxd_nacl_processed"), "refined.refl")
    )

    # Reduce number of reflections
    reflections = reflections.select(reflections["intensity.sum.value"] > 30000)

    refl_path = join(tmp_path, "refined.refl")
    reflections.as_msgpack_file(refl_path)

    # Shoebox summation
    result = subprocess.run(
        [shutil.which("dials.tof_integrate"), expt_path, refl_path],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # Shoebox summation with seed_skewness mask
    result = subprocess.run(
        [
            shutil.which("dials.tof_integrate"),
            expt_path,
            refl_path,
            "mask=seed_skewness",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # 1D profile fitting
    result = subprocess.run(
        [shutil.which("dials.tof_integrate"), expt_path, refl_path, "method=profile1d"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # 3D profile fitting
    result = subprocess.run(
        [shutil.which("dials.tof_integrate"), expt_path, refl_path, "method=profile3d"],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
