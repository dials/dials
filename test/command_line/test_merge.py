"""Tests for dials.merge command line program."""

from __future__ import absolute_import, division, print_function

import procrunner
import pytest
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from iotbx import mtz


def validate_mtz(mtz_file, expected_labels, unexpected_labels=None):
    assert mtz_file.check()
    m = mtz.object(mtz_file.strpath)

    assert m.as_miller_arrays()[0].info().wavelength == pytest.approx(0.6889)
    labels = set()
    for ma in m.as_miller_arrays(merge_equivalents=False):
        labels.update(ma.info().labels)
    for l in expected_labels:
        assert l in labels
    if unexpected_labels:
        for l in unexpected_labels:
            assert l not in labels


@pytest.mark.parametrize("anomalous", [True, False])
@pytest.mark.parametrize("truncate", [True, False])
def test_merge(dials_data, tmpdir, anomalous, truncate):
    """Test the command line script with LCY data"""
    # Main options: truncate on/off, anomalous on/off

    mean_labels = ["IMEAN", "SIGIMEAN"]
    anom_labels = ["I(+)", "I(-)", "SIGI(+)", "SIGI(-)"]
    amp_labels = ["F", "SIGF"]
    anom_amp_labels = ["F(+)", "SIGF(+)", "F(-)", "SIGF(-)"]

    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl")
    expts = location.join("scaled_20_25.expt")

    mtz_file = tmpdir.join("merge-%s-%s.mtz" % (anomalous, truncate))

    command = [
        "dials.merge",
        refls,
        expts,
        "truncate=%s" % truncate,
        "anomalous=%s" % anomalous,
        "output.mtz=%s" % mtz_file.strpath,
        "project_name=ham",
        "crystal_name=jam",
        "dataset_name=spam",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    expected_labels = mean_labels
    unexpected_labels = []

    if truncate:
        expected_labels += amp_labels
    else:
        unexpected_labels += amp_labels

    if anomalous:
        expected_labels += anom_labels
    else:
        unexpected_labels += anom_labels

    if anomalous and truncate:
        expected_labels += anom_amp_labels
    else:
        unexpected_labels += anom_amp_labels

    validate_mtz(mtz_file, expected_labels, unexpected_labels)


def test_merge_multi_wavelength(dials_data, tmpdir):
    """Test that merge handles multi-wavelength data suitably - should be
    exported into an mtz with seprate columns for each wavelength."""

    mean_labels = ["%sIMEAN_WAVE%s" % (pre, i) for i in [1, 2] for pre in ["", "SIG"]]
    anom_labels = [
        "%sI_WAVE%s(%s)" % (pre, i, sgn)
        for i in [1, 2]
        for pre in ["", "SIG"]
        for sgn in ["+", "-"]
    ]
    amp_labels = ["%sF_WAVE%s" % (pre, i) for i in [1, 2] for pre in ["", "SIG"]]
    anom_amp_labels = [
        "%sF_WAVE%s(%s)" % (pre, i, sgn)
        for i in [1, 2]
        for pre in ["", "SIG"]
        for sgn in ["+", "-"]
    ]

    location = dials_data("l_cysteine_4_sweeps_scaled")
    refl1 = location.join("scaled_30.refl").strpath
    expt1 = location.join("scaled_30.expt").strpath
    refl2 = location.join("scaled_35.refl").strpath
    expt2 = location.join("scaled_35.expt").strpath
    expts1 = ExperimentListFactory.from_json_file(expt1, check_format=False)
    expts1[0].beam.set_wavelength(0.5)
    expts2 = ExperimentListFactory.from_json_file(expt2, check_format=False)
    expts1.extend(expts2)

    tmp_expt = tmpdir.join("tmp.expt").strpath
    expts1.as_json(tmp_expt)

    reflections1 = flex.reflection_table.from_pickle(refl1)
    reflections2 = flex.reflection_table.from_pickle(refl2)
    # first need to resolve identifiers - usually done on loading
    reflections2["id"] = flex.int(reflections2.size(), 1)
    del reflections2.experiment_identifiers()[0]
    reflections2.experiment_identifiers()[1] = "3"
    reflections1.extend(reflections2)

    tmp_refl = tmpdir.join("tmp.refl").strpath
    reflections1.as_pickle(tmp_refl)

    # Can now run after creating our 'fake' multiwavelength dataset
    command = ["dials.merge", tmp_refl, tmp_expt, "truncate=True", "anomalous=True"]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("merged.mtz").check()
    m = mtz.object(tmpdir.join("merged.mtz").strpath)
    labels = []
    for ma in m.as_miller_arrays(merge_equivalents=False):
        labels.extend(ma.info().labels)
    assert all(x in labels for x in mean_labels)
    assert all(x in labels for x in anom_labels)
    assert all(x in labels for x in amp_labels)
    assert all(x in labels for x in anom_amp_labels)

    # 5 miller arrays for each dataset
    assert m.as_miller_arrays()[0].info().wavelength == pytest.approx(0.5)
    assert m.as_miller_arrays()[5].info().wavelength == pytest.approx(0.6889)


def test_suitable_exit_for_bad_input_from_single_dataset(dials_data, tmpdir):
    location = dials_data("vmxi_proteinase_k_sweeps")

    command = [
        "dials.merge",
        location.join("experiments_0.json"),
        location.join("reflections_0.pickle"),
    ]

    # unscaled data
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode
    assert result.stderr.startswith(b"Sorry")


def test_suitable_exit_for_bad_input_with_more_than_one_reflection_table(
    dials_data, tmpdir
):
    location = dials_data("vmxi_proteinase_k_sweeps")

    command = [
        "dials.merge",
        location.join("experiments_0.json"),
        location.join("reflections_0.pickle"),
        location.join("experiments_1.json"),
        location.join("reflections_1.pickle"),
    ]

    # more than one reflection table.
    result = procrunner.run(command, working_directory=tmpdir)
    assert result.returncode
    assert result.stderr.startswith(b"Sorry")
