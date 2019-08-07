from __future__ import absolute_import, division, print_function

import procrunner
import pytest


@pytest.fixture
def model(tmpdir):
    from dials.array_family import flex
    from dials.algorithms.background.gmodel import StaticBackgroundModel
    import six.moves.cPickle as pickle

    ysize = 2527
    xsize = 2463
    data = flex.double(flex.grid(ysize, xsize), 1)
    model = StaticBackgroundModel()
    model.add(data)

    model_file = tmpdir.join("model.pickle")
    with model_file.open("wb") as fh:
        pickle.dump(model, fh, pickle.HIGHEST_PROTOCOL)
    return model_file


def test_simple(dials_data, model, tmpdir):
    path = dials_data("centroid_test_data")
    experiments = path.join("experiments.json")

    reflns_simple = tmpdir.join("simple").join("observations.refl")
    reflns_g_simple = tmpdir.join("gmodel_simple").join("observations.refl")
    reflns_simple.dirpath().ensure(dir=1)
    reflns_g_simple.dirpath().ensure(dir=1)

    result = procrunner.run(
        [
            "dials.integrate",
            experiments.strpath,
            "profile.fitting=False",
            "background.algorithm=simple",
            "background.simple.outlier.algorithm=null",
            "output.reflections=" + reflns_simple.strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert reflns_simple.check()

    result = procrunner.run(
        [
            "dials.integrate",
            experiments.strpath,
            "profile.fitting=False",
            "background.algorithm=gmodel",
            "background.gmodel.robust.algorithm=False",
            "background.gmodel.model=model.pickle",
            "output.reflections=" + reflns_g_simple.strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert reflns_g_simple.check()

    from dials.array_family import flex

    reflections1 = flex.reflection_table.from_file(reflns_simple.strpath)
    reflections3 = flex.reflection_table.from_file(reflns_g_simple.strpath)
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


def test_robust(dials_data, model, tmpdir):
    path = dials_data("centroid_test_data")
    experiments = path.join("experiments.json")

    reflns_robust = tmpdir.join("robust").join("observations.refl")
    reflns_g_robust = tmpdir.join("gmodel_robust").join("observations.refl")
    reflns_robust.dirpath().ensure(dir=1)
    reflns_g_robust.dirpath().ensure(dir=1)

    result = procrunner.run(
        [
            "dials.integrate",
            experiments.strpath,
            "profile.fitting=False",
            "background.algorithm=glm",
            "output.reflections=" + reflns_robust.strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert reflns_robust.check()

    result = procrunner.run(
        [
            "dials.integrate",
            experiments.strpath,
            "profile.fitting=False",
            "background.algorithm=gmodel",
            "background.gmodel.robust.algorithm=True",
            "background.gmodel.model=model.pickle",
            "output.reflections=" + reflns_g_robust.strpath,
        ],
        working_directory=tmpdir.strpath,
    )
    assert not result.returncode and not result.stderr
    assert reflns_g_robust.check()

    from dials.array_family import flex

    reflections2 = flex.reflection_table.from_file(reflns_robust.strpath)
    reflections4 = flex.reflection_table.from_file(reflns_g_robust.strpath)
    assert len(reflections2) == len(reflections4)

    flag = flex.reflection_table.flags.integrated_sum
    integrated2 = reflections2.select(reflections2.get_flags(flag, all=True))
    integrated4 = reflections4.select(reflections4.get_flags(flag, all=True))

    assert len(integrated2) > 0
    assert len(integrated2) == len(integrated4)

    mean_bg2 = integrated2["background.mean"]
    mean_bg4 = integrated4["background.mean"]
    scale4 = integrated4["background.scale"]
    assert (scale4 > 0).count(False) == 0
