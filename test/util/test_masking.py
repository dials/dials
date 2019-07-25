from __future__ import absolute_import, division, print_function

import json
import os
import pytest

import libtbx
from libtbx.easy_run import fully_buffered
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.serialize import load

from dials.array_family import flex
from dials.algorithms.shadowing.filter import filter_shadowed_reflections


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


def test_shadow_plot(dials_regression, run_in_tmpdir):
    path = os.path.join(
        dials_regression, "shadow_test_data/DLS_I04_SmarGon/Th_3_O45_C45_P48_1_0500.cbf"
    )

    fully_buffered("dials.import %s" % path).raise_if_errors()
    fully_buffered("dials.shadow_plot imported.expt json=shadow.json").raise_if_errors()
    assert os.path.exists("scan_shadow_plot.png")
    assert os.path.exists("shadow.json")
    with open("shadow.json", "rb") as f:
        d = json.load(f)
        assert set(d) == {"fraction_shadowed", "scan_points"}
        assert d["fraction_shadowed"] == pytest.approx([0.06856597327776767], 2e-4)

    fully_buffered(
        "dials.shadow_plot imported.expt mode=2d plot=shadow_2d.png"
    ).raise_if_errors()
    assert os.path.exists("shadow_2d.png")


def test_filter_shadowed_reflections(dials_regression):
    experiments_json = os.path.join(
        dials_regression, "shadow_test_data/DLS_I04_SmarGon/experiments.json"
    )
    predicted_pickle = os.path.join(
        dials_regression, "shadow_test_data/DLS_I04_SmarGon/predicted.pickle"
    )

    experiments = load.experiment_list(experiments_json, check_format=True)
    predicted = flex.reflection_table.from_file(predicted_pickle)

    for experiment_goniometer in (True, False):
        shadowed = filter_shadowed_reflections(
            experiments, predicted, experiment_goniometer=experiment_goniometer
        )
        assert shadowed.count(True) == 17
        assert shadowed.count(False) == 674
