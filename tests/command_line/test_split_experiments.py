"""Tests for dials.split_experiments when experiment ids are set"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from dxtbx.model import Beam, Experiment, ExperimentList
from dxtbx.serialize import load

from dials.array_family import flex


def generate_exp(wavelength=1):
    """Generate an experiment containing a beam with a given wavelength."""
    beam = Beam(direction=(0.0, 0.0, 1.0), wavelength=wavelength)
    exp = Experiment(beam=beam)
    return exp


@pytest.mark.parametrize("with_identifiers", ["True", "False"])
@pytest.mark.parametrize("option", ["chunk_size=2", 'chunk_sizes="2 2 1"'])
def test_split_chunk_sizes(tmp_path, option, with_identifiers):
    ids = list(range(0, 5))
    experiments = ExperimentList()
    reflections = flex.reflection_table()
    reflections["id"] = flex.int(ids)
    reflections["intensity"] = flex.double([(i + 1) * 100.0 for i in ids])

    for i in ids:
        exp = generate_exp()
        if with_identifiers:
            exp.identifier = str(i)
            reflections.experiment_identifiers()[i] = str(i)
        experiments.append(exp)

    experiments.as_json(tmp_path / "tmp.expt")
    reflections.as_file(tmp_path / "tmp.refl")

    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            tmp_path / "tmp.expt",
            tmp_path / "tmp.refl",
            option,
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for j, n, intensities in zip(
        [0, 1, 2], [2, 2, 1], [[100.0, 200.0], [300.0, 400.0], [500.0]]
    ):
        assert (tmp_path / f"split_{j}.refl").is_file()
        assert (tmp_path / f"split_{j}.expt").is_file()
        expts = load.experiment_list(tmp_path / f"split_{j}.expt", check_format=False)
        assert len(expts) == n
        refls = flex.reflection_table.from_file(tmp_path / f"split_{j}.refl")
        assert list(set(refls["id"])) == list(range(0, n))
        assert list(refls["intensity"]) == intensities
        refls.assert_experiment_identifiers_are_consistent(expts)


def test_split_by_wavelength(tmp_path):
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

    experiments.as_json(tmp_path / "tmp.expt")
    reflections.as_file(tmp_path / "tmp.refl")

    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            "tmp.expt",
            "tmp.refl",
            "by_wavelength=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    for i, (wl, ids, intensity) in enumerate(
        zip([0.5, 1.0], ["1", "0"], [200.0, 100.0])
    ):
        assert (tmp_path / f"split_{i}.expt").is_file()
        assert (tmp_path / f"split_{i}.refl").is_file()
        exp_single = load.experiment_list(
            tmp_path / f"split_{i}.expt", check_format=False
        )
        ref_single = flex.reflection_table.from_file(tmp_path / f"split_{i}.refl")
        assert exp_single[0].beam.get_wavelength() == wl
        assert exp_single[0].identifier == ids
        id_ = ref_single["id"][0]
        assert ref_single.experiment_identifiers()[id_] == ids
        assert list(ref_single["intensity"]) == [intensity]

    # Now test for successful error handling if no identifiers set.
    experiments[0].identifier = ""
    experiments[1].identifier = ""
    experiments.as_json(tmp_path / "tmp.expt")
    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            "tmp.expt",
            "tmp.refl",
            "by_wavelength=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode == 1
    assert result.stderr.startswith(b"Sorry")

    experiments[0].identifier = "0"
    experiments[1].identifier = "1"
    del reflections.experiment_identifiers()[0]
    del reflections.experiment_identifiers()[1]
    experiments.as_json(tmp_path / "tmp.expt")
    reflections.as_file(tmp_path / "tmp.refl")
    result = subprocess.run(
        [
            shutil.which("dials.split_experiments"),
            "tmp.expt",
            "tmp.refl",
            "by_wavelength=True",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode == 1
    assert result.stderr.startswith(b"Sorry")
