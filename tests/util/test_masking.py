from __future__ import annotations

import json
import os

import procrunner
import pytest

import libtbx
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load

import dials.util.masking
from dials.algorithms.shadowing.filter import filter_shadowed_reflections
from dials.array_family import flex


@pytest.mark.parametrize(
    "path,count_only_shadow,count_mask_shadow,count_mask_no_shadow",
    [
        (
            "shadow_test_data/DLS_I04_SmarGon/Th_3_O45_C45_P48_1_0500.cbf",
            pytest.approx(426758, abs=1e2),
            pytest.approx(917940, abs=1e2),
            pytest.approx(528032, abs=1e2),
        ),
        (
            "shadow_test_data/DLS_I03_SmarGon/protk_2_0600.cbf",
            pytest.approx(519100, abs=1e2),
            pytest.approx(1002068, abs=1e2),
            pytest.approx(527314, abs=1e2),
        ),
    ],
)
def test_dynamic_shadowing(
    path, count_only_shadow, count_mask_shadow, count_mask_no_shadow, dials_regression
):
    path = os.path.join(dials_regression, path)
    assert os.path.exists(path), path
    for shadowing in (libtbx.Auto, True, False):
        format_kwargs = {"dynamic_shadowing": shadowing}
        experiments = ExperimentListFactory.from_filenames(
            [path], format_kwargs=format_kwargs
        )
        imageset = experiments.imagesets()[0]
        detector = imageset.get_detector()
        scan = imageset.get_scan()
        masker = imageset.masker()
        if shadowing is not False:
            assert masker is not None
            mask = masker.get_mask(detector, scan.get_oscillation()[0])
            assert len(mask) == len(detector)
            # only shadowed pixels masked
            assert mask[0].count(False) == count_only_shadow, (
                mask[0].count(False),
                count_only_shadow,
            )
        mask = imageset.get_mask(0)

        # dead pixels, pixels in gaps, etc also masked
        if shadowing is libtbx.Auto or shadowing is True:
            assert mask[0].count(False) == count_mask_shadow, (
                mask[0].count(False),
                count_mask_shadow,
            )
        else:
            assert mask[0].count(False) == count_mask_no_shadow, (
                mask[0].count(False),
                count_mask_no_shadow,
            )


@pytest.mark.xfail(reason="Failing due to deprecation warning in output")
def test_shadow_plot(dials_data, tmp_path):
    result = procrunner.run(
        (
            "dials.import",
            dials_data("image_examples", pathlib=True) / "DLS_I03_smargon_0001.cbf.gz",
        ),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    result = procrunner.run(
        ("dials.shadow_plot", "imported.expt", "json=shadow.json"),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("scan_shadow_plot.png").is_file()
    assert tmp_path.joinpath("shadow.json").is_file()
    d = json.loads(tmp_path.joinpath("shadow.json").read_text())
    assert set(d) == {"fraction_shadowed", "scan_points"}
    assert d["fraction_shadowed"] == pytest.approx([0.003016, 0.003141], 2e-4)
    result = procrunner.run(
        ("dials.shadow_plot", "imported.expt", "mode=2d", "plot=shadow_2d.png"),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("shadow_2d.png").is_file()


def test_filter_shadowed_reflections(dials_regression):
    experiments_json = os.path.join(
        dials_regression, "shadow_test_data", "DLS_I04_SmarGon", "experiments.json"
    )
    predicted_pickle = os.path.join(
        dials_regression, "shadow_test_data", "DLS_I04_SmarGon", "predicted.pickle"
    )

    experiments = load.experiment_list(experiments_json, check_format=True)
    predicted = flex.reflection_table.from_file(predicted_pickle)

    for experiment_goniometer in (True, False):
        shadowed = filter_shadowed_reflections(
            experiments, predicted, experiment_goniometer=experiment_goniometer
        )
        assert shadowed.count(True) == 17
        assert shadowed.count(False) == 674


def test_lru_equality_cache_basic():
    callargs = []
    result = [None]

    def _callappend(*arg):
        callargs.append(arg)
        return result[-1]

    fun = dials.util.masking.lru_equality_cache(maxsize=20)(_callappend)
    hits, misses, maxsize, currsize = fun.cache_info()
    assert hits == 0
    assert misses == 0
    assert maxsize == 20
    assert currsize == 0

    for i in range(100):
        fun(i)
        fun(i)
    hits, misses, maxsize, currsize = fun.cache_info()
    assert hits == 100
    assert misses == 100
    assert maxsize == 20
    assert currsize == 20
    assert len(callargs) == 100
    assert len(set(callargs)) == 100
    for i in range(80, 100):
        fun(i)
    assert fun.cache_info() == (120, 100, 20, 20)
    assert len(callargs) == 100

    # Change the return value
    assert fun(1) is None
    result.append("test_result")
    assert fun(1) is None
    assert fun.__wrapped__(1) == "test_result"
    assert fun(2) == "test_result"


def test_lru_equality_cache_id():
    callargs = []

    def _callappend(*arg):
        callargs.append(arg)

    fun = dials.util.masking.lru_equality_cache(maxsize=1)(_callappend)

    class EqTester:
        def __init__(self, a):
            self.a = a

        def __eq__(self, other):
            return other.a == self.a

        def __hash__(self):
            return hash(id(self))

    a, b = EqTester(-424242), EqTester(-424242)
    assert a == b
    assert hash(a) != hash(b)

    fun(a)
    fun(b)
    assert fun.cache_info() == (1, 1, 1, 1)


def test_generate_mask(dials_data):
    imageset = load.imageset(
        dials_data("centroid_test_data", pathlib=True) / "sweep.json"
    )
    params = dials.util.masking.phil_scope.extract()
    params.border = 2
    params.d_min = 1.5
    params.d_max = 40
    params.untrusted[0].rectangle = (500, 600, 700, 800)
    params.untrusted[0].circle = (1000, 1200, 100)
    params.untrusted[0].polygon = (1500, 1500, 1600, 1600, 1700, 1500, 1500, 1500)
    params.untrusted[0].pixel = (1183, 1383)
    params.resolution_range = [(3.1, 3.2), (5.1, 5.2)]
    params.ice_rings.filter = True
    mask = dials.util.masking.generate_mask(imageset, params)
    assert len(mask) == len(imageset.get_detector())
    for m, im in zip(mask, imageset.get_raw_data(0)):
        assert m.all() == im.all()
    assert mask[0].count(False) == 4060817


@pytest.mark.parametrize(
    "disable_parallax_correction,expected", [(False, 1427394), (True, 1432002)]
)
def test_generate_mask_parallax_correction(
    disable_parallax_correction, expected, dials_data
):
    expts = ExperimentListFactory.from_filenames(
        sorted(dials_data("centroid_test_data", pathlib=True).glob("*.cbf"))
    )
    imageset = expts[0].imageset
    params = dials.util.masking.phil_scope.extract()
    params.disable_parallax_correction = disable_parallax_correction
    params.d_min = 1.2
    params.d_max = 40
    params.untrusted = []
    mask = dials.util.masking.generate_mask(imageset, params)
    assert len(mask) == len(imageset.get_detector())
    for m, im in zip(mask, imageset.get_raw_data(0)):
        assert m.all() == im.all()
    assert mask[0].count(False) == expected
