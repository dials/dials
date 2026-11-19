from __future__ import annotations

import shutil
import subprocess

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load

from dials.array_family.flex import reflection_table


@pytest.fixture(scope="module")
def tof_indexing_input(dials_data, tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("tof_indexing_input")
    image_file = dials_data("isis_sxd_example_data") / "sxd_nacl_run.nxs"
    experiments = ExperimentListFactory.from_filenames([image_file])
    experiments_file = tmp_path / "imported.expt"
    experiments.as_file(experiments_file)

    # Filter reflections
    reflections = reflection_table.from_msgpack_file(
        dials_data("isis_sxd_nacl_processed") / "strong.refl"
    )
    _, _, _, _, z0, z1 = reflections["bbox"].parts()
    reflections = reflections.select((z1 - z0) > 20)
    reflections_file = tmp_path / "strong.refl"
    reflections.as_msgpack_file(reflections_file)

    return experiments_file, reflections_file


# Each method is a separate test so that they can run in parallel under xdist
@pytest.mark.parametrize(
    "method", ["fft1d", "fft3d", "real_space_grid_search", "low_res_spot_match"]
)
def test_tof_indexing(tof_indexing_input, tmp_path, method):
    """
    Run indexing methods expected to work with ToF data
    """
    experiments_file, reflections_file = tof_indexing_input
    subprocess.run(
        (
            shutil.which("dials.index"),
            experiments_file,
            reflections_file,
            "unit_cell=5.64,5.64,5.64,90,90,90",
            "space_group=225",
            f"indexing.method={method}",
        ),
        cwd=tmp_path,
        capture_output=True,
    ).check_returncode()

    experiments = load.experiment_list(tmp_path / "indexed.expt")
    assert len(experiments) == 1
