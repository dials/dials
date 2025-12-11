from __future__ import annotations

import pickle
import shutil
import subprocess

import pytest

from dials.algorithms.background.gmodel import StaticBackgroundModel
from dials.array_family import flex


@pytest.fixture
def model(tmp_path):
    ysize = 2527
    xsize = 2463
    data = flex.double(flex.grid(ysize, xsize), 1)
    model = StaticBackgroundModel()
    model.add(data)

    model_file = tmp_path / "model.pickle"
    with model_file.open("wb") as fh:
        pickle.dump(model, fh, pickle.HIGHEST_PROTOCOL)
    return model_file


def test_simple(dials_data, model, tmp_path):
    experiments = dials_data("centroid_test_data") / "experiments.json"

    reflns_simple = tmp_path / "simple" / "observations.refl"
    reflns_g_simple = tmp_path / "gmodel_simple" / "observations.refl"
    reflns_simple.parent.mkdir()
    reflns_g_simple.parent.mkdir()

    result = subprocess.run(
        [
            shutil.which("dials.integrate"),
            "nproc=1",
            experiments,
            "profile.fitting=False",
            "background.algorithm=simple",
            "background.simple.outlier.algorithm=null",
            f"output.reflections={reflns_simple}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert reflns_simple.is_file()

    result = subprocess.run(
        [
            shutil.which("dials.integrate"),
            "nproc=1",
            experiments,
            "profile.fitting=False",
            "background.algorithm=gmodel",
            "background.gmodel.robust.algorithm=False",
            "background.gmodel.model=model.pickle",
            f"output.reflections={reflns_g_simple}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert reflns_g_simple.is_file()

    reflections1 = flex.reflection_table.from_file(reflns_simple)
    reflections3 = flex.reflection_table.from_file(reflns_g_simple)
    assert len(reflections1) == len(reflections3)

    flag = flex.reflection_table.flags.integrated_sum
    integrated1 = reflections1.select(reflections1.get_flags(flag, all=True))
    integrated3 = reflections3.select(reflections3.get_flags(flag, all=True))

    assert len(integrated1) > 0
    assert len(integrated1) == len(integrated3)

    mean_bg1 = integrated1["background.mean"]
    mean_bg3 = integrated3["background.mean"]
    scale3 = integrated3["background.scale"]

    diff1 = flex.abs(mean_bg1 - mean_bg3)
    assert (scale3 > 0).count(False) == 0
    assert (diff1 < 1e-5).count(False) == 0


def test_robust(dials_data, model, tmp_path):
    experiments = dials_data("centroid_test_data") / "experiments.json"

    reflns_robust = tmp_path / "robust" / "observations.refl"
    reflns_g_robust = tmp_path / "gmodel_robust" / "observations.refl"
    reflns_robust.parent.mkdir()
    reflns_g_robust.parent.mkdir()

    result = subprocess.run(
        [
            shutil.which("dials.integrate"),
            "nproc=1",
            experiments,
            "profile.fitting=False",
            "background.algorithm=glm",
            f"output.reflections={reflns_robust}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert reflns_robust.is_file()

    result = subprocess.run(
        [
            shutil.which("dials.integrate"),
            "nproc=1",
            experiments,
            "profile.fitting=False",
            "background.algorithm=gmodel",
            "background.gmodel.robust.algorithm=True",
            "background.gmodel.model=model.pickle",
            f"output.reflections={reflns_g_robust}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert reflns_g_robust.is_file()

    reflections2 = flex.reflection_table.from_file(reflns_robust)
    reflections4 = flex.reflection_table.from_file(reflns_g_robust)
    assert len(reflections2) == len(reflections4)

    flag = flex.reflection_table.flags.integrated_sum
    integrated2 = reflections2.select(reflections2.get_flags(flag, all=True))
    integrated4 = reflections4.select(reflections4.get_flags(flag, all=True))

    assert len(integrated2) > 0
    assert len(integrated2) == len(integrated4)

    scale4 = integrated4["background.scale"]
    assert (scale4 > 0).count(False) == 0
