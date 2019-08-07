from __future__ import absolute_import, division, print_function

"""Tests for dials.split_experiments when experiment ids are set"""

import procrunner
from dials.array_family import flex
from dxtbx.model import Beam, Experiment, ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory


def generate_exp(wavelength=1):
    """Generate an experiment containing a beam with a given wavelength."""
    beam = Beam(direction=(0.0, 0.0, 1.0), wavelength=wavelength)
    exp = Experiment(beam=beam)
    return exp


def test_split_by_wavelength(tmpdir):
    """Test the split_by_wavelength option of dials.split_experiments"""
    experiments = ExperimentList()
    exp = generate_exp(wavelength=1.0)
    exp.identifier = "0"
    experiments.append(exp)
    exp = generate_exp(wavelength=0.5)
    exp.identifier = "1"
    experiments.append(exp)

    reflections = flex.reflection_table()
    reflections["id"] = flex.int([0, 1])
    reflections["intensity"] = flex.double([100.0, 200.0])
    reflections.experiment_identifiers()[0] = "0"
    reflections.experiment_identifiers()[1] = "1"

    experiments.as_json(tmpdir.join("tmp.expt").strpath)
    reflections.as_pickle(tmpdir.join("tmp.refl").strpath)

    result = procrunner.run(
        ["dials.split_experiments", "tmp.expt", "tmp.refl", "by_wavelength=True"],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    for i, (wl, ids, intensity) in enumerate(
        zip([0.5, 1.0], ["1", "0"], [200.0, 100.0])
    ):
        assert tmpdir.join("split_%d.expt" % i).check()
        assert tmpdir.join("split_%d.refl" % i).check()
        exp_single = ExperimentListFactory.from_json_file(
            tmpdir.join("split_%d.expt" % i).strpath, check_format=False
        )
        ref_single = flex.reflection_table.from_pickle(
            tmpdir.join("split_%d.refl" % i).strpath
        )
        assert exp_single[0].beam.get_wavelength() == wl
        assert exp_single[0].identifier == ids
        id_ = ref_single["id"][0]
        assert ref_single.experiment_identifiers()[id_] == ids
        assert list(ref_single["intensity"]) == [intensity]

    # Now test for successful error handling if no identifiers set.
    experiments[0].identifier = ""
    experiments[1].identifier = ""
    experiments.as_json(tmpdir.join("tmp.expt").strpath)
    result = procrunner.run(
        ["dials.split_experiments", "tmp.expt", "tmp.refl", "by_wavelength=True"],
        working_directory=tmpdir,
    )
    assert result.returncode == 1
    assert result.stderr.startswith("Sorry")

    experiments[0].identifier = "0"
    experiments[1].identifier = "1"
    del reflections.experiment_identifiers()[0]
    del reflections.experiment_identifiers()[1]
    experiments.as_json(tmpdir.join("tmp.expt").strpath)
    reflections.as_pickle(tmpdir.join("tmp.refl").strpath)
    result = procrunner.run(
        ["dials.split_experiments", "tmp.expt", "tmp.refl", "by_wavelength=True"],
        working_directory=tmpdir,
    )
    assert result.returncode == 1
    assert result.stderr.startswith("Sorry")
