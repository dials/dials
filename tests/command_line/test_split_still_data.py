from __future__ import annotations

import os
import pathlib

import pytest

from dxtbx.serialize import load

import dials.command_line.split_still_data as split


@pytest.mark.parametrize("use_yaml", [True, False])
def test_split_still_data(dials_data, run_in_tmp_path, use_yaml):
    data = dials_data("cunir_serial_processed", pathlib=True)
    args = [
        os.fspath(data / "integrated.expt"),
        os.fspath(data / "integrated.refl"),
        "nproc=1",
    ]
    if use_yaml:
        images = os.fspath(dials_data("cunir_serial", pathlib=True))
        yml = f"""
---
metadata:
  timepoint:
    {images}/merlin0047_#####.cbf : 'repeat=2'
grouping:
  group_by:
    values:
      - timepoint
"""
        test_yaml = run_in_tmp_path / "tmp.yaml"
        with open(test_yaml, "w") as f:
            f.write(yml)
        args.append(f"grouping={test_yaml}")
    else:
        args.append("series_repeat=2")
    split.run(args=args)
    assert pathlib.Path("group_0_0.expt").is_file()
    assert pathlib.Path("group_0_0.refl").is_file()
    assert pathlib.Path("group_1_0.expt").is_file()
    assert pathlib.Path("group_1_0.refl").is_file()
    expts1 = load.experiment_list("group_0_0.expt", check_format=False)
    assert len(expts1) == 3
    # not old style elist datastructures (no scan, single imageset)
    assert expts1[0].imageset.get_path(0).split("_")[-1].rstrip(".cbf") == "17000"
    assert expts1[1].imageset.get_path(0).split("_")[-1].rstrip(".cbf") == "17002"
    assert expts1[2].imageset.get_path(0).split("_")[-1].rstrip(".cbf") == "17004"
    expts2 = load.experiment_list("group_1_0.expt", check_format=False)
    assert len(expts2) == 2
    assert expts2[0].imageset.get_path(0).split("_")[-1].rstrip(".cbf") == "17001"
    assert expts2[1].imageset.get_path(0).split("_")[-1].rstrip(".cbf") == "17003"
