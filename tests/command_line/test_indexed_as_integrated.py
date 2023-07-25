from __future__ import annotations

import shutil
import subprocess

import iotbx.merging_statistics


def test_indexed_as_integrated(dials_data, tmp_path):
    data_dir = dials_data("insulin_processed", pathlib=True)
    refl = data_dir / "indexed.refl"
    expt = data_dir / "indexed.expt"

    command = [
        shutil.which("dials.indexed_as_integrated"),
        refl,
        expt,
        "output.reflections=output.refl",
    ]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "output.refl").is_file()

    # now make sure we can run dials.symmetry and dials.scale successfully
    sym_command = [shutil.which("dials.symmetry"), tmp_path / "output.refl", expt]
    result = subprocess.run(sym_command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "symmetrized.refl").is_file()
    assert (tmp_path / "symmetrized.expt").is_file()

    scale_command = [
        shutil.which("dials.scale"),
        tmp_path / "symmetrized.refl",
        tmp_path / "symmetrized.expt",
        "unmerged_mtz=scaled.mtz",
        "d_min=2.0",
    ]
    result = subprocess.run(scale_command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "scaled.refl").is_file()
    assert (tmp_path / "scaled.expt").is_file()
    assert (tmp_path / "scaled.mtz").is_file()

    i_obs = iotbx.merging_statistics.select_data(str(tmp_path / "scaled.mtz"), None)
    result = iotbx.merging_statistics.dataset_statistics(
        i_obs=i_obs,
        use_internal_variance=False,
        eliminate_sys_absent=False,
    )

    assert result.overall.cc_one_half > 0.8
