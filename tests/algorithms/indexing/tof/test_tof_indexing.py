from __future__ import annotations

import shutil
import subprocess
from os.path import join

from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load

from dials.array_family.flex import reflection_table


def test_tof_indexing(dials_data, tmp_path):
    """
    Run indexing methods expected to work with ToF data
    """

    image_file = join(
        dials_data("isis_sxd_example_data", pathlib=True), "sxd_nacl_run.nxs"
    )
    experiments = ExperimentListFactory.from_filenames([image_file])
    experiments_file = join(tmp_path, "imported.expt")
    experiments.as_file(experiments_file)

    reflections_file = join(
        dials_data("isis_sxd_nacl_processed", pathlib=True), "strong.refl"
    )

    # Filter reflections
    reflections = reflection_table.from_msgpack_file(reflections_file)
    _, _, _, _, z0, z1 = reflections["bbox"].parts()
    reflections = reflections.select((z1 - z0) > 20)
    reflections_file = join(tmp_path, "strong.refl")
    reflections.as_msgpack_file(reflections_file)

    methods = ["fft1d", "fft3d", "real_space_grid_search", "low_res_spot_match"]

    for method in methods:
        expt_path = join(tmp_path, f"{method}.expt")
        refl_path = join(tmp_path, f"{method}.refl")
        subprocess.run(
            (
                shutil.which("dials.index"),
                experiments_file,
                reflections_file,
                "unit_cell=5.64,5.64,5.64,90,90,90",
                "space_group=225",
                f"indexing.method={method}",
                f"output.experiments={expt_path}",
                f"output.reflections={refl_path}",
            ),
            cwd=tmp_path,
            capture_output=True,
        ).check_returncode()

        experiments = load.experiment_list(expt_path)
        assert len(experiments) == 1
