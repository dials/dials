from __future__ import annotations

import json
import pathlib
import shutil
import subprocess

import pytest

from dxtbx.serialize import load

from dials.array_family import flex
from dials.command_line.ssx_integrate import run_integration, working_phil
from dials.util.options import ArgumentParser

# Note that tests are grouped and run serially, to stop many processes trying to
# extract data from images at same time, which appears to lead to race
# conditions in some CI jobs.


@pytest.mark.xdist_group(name="group1")
def test_ssx_integrate_fullprocess(dials_data, tmp_path):
    # Download data set and the internally referenced images
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    dials_data("cunir_serial", pathlib=True)
    pathlib.Path.mkdir(tmp_path / "nuggets")
    result = subprocess.run(
        [
            shutil.which("dials.ssx_integrate"),
            ssx / "indexed.refl",
            ssx / "indexed.expt",
            "nproc=1",
            "batch_size=3",
            "output.json=data.json",
            "output.nuggets=nuggets",
            "algorithm=stills",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("integrated_1.refl").is_file()
    assert tmp_path.joinpath("integrated_1.expt").is_file()
    assert tmp_path.joinpath("integrated_2.refl").is_file()
    assert tmp_path.joinpath("integrated_2.expt").is_file()
    assert tmp_path.joinpath("dials.ssx_integrate.html").is_file()
    assert tmp_path.joinpath("data.json").is_file()
    for i in range(1, 6):
        assert tmp_path.joinpath(f"nuggets/nugget_integrated_{i}.json").is_file()


expected_simple1 = {
    "likelihood": 171374.174649,
    "mosaicity": [0.0003630],
}
expected_s1a1 = {
    "likelihood": 171217.954551,
    "mosaicity": [0.0003951, 6.6198e-06],
}
expected_s1a3 = {
    "likelihood": 171209.435802,
    "mosaicity": [0.0004138, 0.0450596, 1.39279e-05],
}
expected_simple6 = {
    "likelihood": 176234.854949,
    "mosaicity": [0.0001368, 0.0003956, 0.0004691],
}
expected_s6a1 = {
    "likelihood": 177094.368288,
    "mosaicity": [1.9031e-05, 0.0003987, 0.0004636, 0.0324567],
}


@pytest.mark.parametrize(
    "model,expected",
    [
        ("simple1", expected_simple1),
        ("simple6", expected_simple6),
        ("simple1angular1", expected_s1a1),
        ("simple1angular3", expected_s1a3),
        ("simple6angular1", expected_s6a1),
    ],
)
@pytest.mark.xdist_group(name="group1")
def test_ssx_integrate_fullprocess_ellipsoid(dials_data, tmp_path, model, expected):
    # Download data set and the internally referenced images
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    dials_data("cunir_serial", pathlib=True)

    expts = load.experiment_list(ssx / "indexed.expt", check_format=False)
    expts[3:4].as_file(tmp_path / "single.expt")
    refls = flex.reflection_table.from_file(ssx / "indexed.refl")
    refls = refls.select_on_experiment_identifiers([expts[3].identifier])
    refls.as_file(tmp_path / "single.refl")

    result = subprocess.run(
        [
            shutil.which("dials.ssx_integrate"),
            tmp_path / "single.refl",
            tmp_path / "single.expt",
            "nproc=1",
            "algorithm=ellipsoid",
            f"ellipsoid.rlp_mosaicity={model}",
            "n_macro_cycles=2",
            f"output.history={tmp_path / 'history.json'}",
            "max_iter=100",
            "LL_tolerance=1e-3",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert tmp_path.joinpath("integrated_1.refl").is_file()
    assert tmp_path.joinpath("integrated_1.expt").is_file()
    assert tmp_path.joinpath("dials.ssx_integrate.html").is_file()
    expts = load.experiment_list(tmp_path / "integrated_1.expt", check_format=False)
    mosaicity = expts[0].profile.mosaicity()
    assert list(mosaicity.values()) == pytest.approx(expected["mosaicity"], abs=1e-5)
    with (tmp_path / "history.json").open("r") as fh:
        data = json.load(fh)
        assert data["0"]["likelihood_per_iteration"][-1][-1] == pytest.approx(
            expected["likelihood"], abs=1e-5
        )


@pytest.mark.parametrize("algorithm,expected_n_refls", [("stills", 614)])
@pytest.mark.xdist_group(name="group1")
def test_ssx_integrate_algorithms(dials_data, algorithm, expected_n_refls):
    # Download data set and the internally referenced images
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    dials_data("cunir_serial", pathlib=True)

    indexed_refl = flex.reflection_table.from_file(
        ssx / "indexed.refl"
    ).split_by_experiment_id()
    indexed_expts = load.experiment_list(ssx / "indexed.expt", check_format=True)

    parser = ArgumentParser(phil=working_phil, check_format=False)
    params, _ = parser.parse_args(args=[], quick_parse=True)

    params.algorithm = algorithm
    params.nproc = 1
    params.image_range = "1:2"

    results = list(run_integration(indexed_refl, indexed_expts, params))
    assert len(results) == 1  # only one batch expected
    experiments = results[0][0]
    reflections = results[0][1]

    assert len(experiments) == 2
    assert len(reflections) == pytest.approx(expected_n_refls, abs=9)
