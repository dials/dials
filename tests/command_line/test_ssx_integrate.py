import procrunner
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
    result = procrunner.run(
        [
            "dev.dials.ssx_integrate",
            str(ssx / "indexed.refl"),
            str(ssx / "indexed.expt"),
            "nproc=1",
            "batch_size=3",
            "output.json=data.json",
            "algorithm=stills",
        ],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "integrated_1.refl").is_file()
    assert (tmp_path / "integrated_1.expt").is_file()
    assert (tmp_path / "integrated_2.refl").is_file()
    assert (tmp_path / "integrated_2.expt").is_file()
    assert (tmp_path / "dials.ssx_integrate.html").is_file()
    assert (tmp_path / "data.json").is_file()


@pytest.mark.xdist_group(name="group1")
def test_ssx_integrate_algorithms(dials_data):
    # Download data set and the internally referenced images
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    dials_data("cunir_serial", pathlib=True)

    indexed_refl = flex.reflection_table.from_file(
        ssx / "indexed.refl"
    ).split_by_experiment_id()
    indexed_expts = load.experiment_list(ssx / "indexed.expt", check_format=True)

    parser = ArgumentParser(phil=working_phil, check_format=False)
    params, _ = parser.parse_args(args=[], quick_parse=True)

    params.algorithm = "stills"
    params.nproc = 1
    params.image_range = "1:2"

    results = list(run_integration(indexed_refl, indexed_expts, params))
    assert len(results) == 1
    experiments = results[0][0]
    reflections = results[0][1]

    assert len(experiments) == 2
    expected_n_refls = 614
    assert len(reflections) == pytest.approx(expected_n_refls, abs=9)

    params.algorithm = "ellipsoid"
    results = list(run_integration(indexed_refl, indexed_expts, params))
    assert len(results) == 1
    experiments = results[0][0]
    reflections = results[0][1]
    assert len(experiments) == 2
    assert experiments.profiles()[0].name == "ellipsoid"
    expected_n_refls = 1258
    assert len(reflections) == pytest.approx(expected_n_refls, abs=9)
    assert indexed_expts[0].crystal != experiments[0].crystal

    params.profile.ellipsoid.unit_cell.fixed = True
    params.profile.ellipsoid.orientation.fixed = True

    results = list(run_integration(indexed_refl, indexed_expts, params))
    assert len(results) == 1
    experiments = results[0][0]
    reflections = results[0][1]
    assert len(experiments) == 2
    assert experiments.profiles()[0].name == "ellipsoid"
    expected_n_refls = 1266
    assert len(reflections) == pytest.approx(expected_n_refls, abs=9)
    assert indexed_expts[0].crystal == experiments[0].crystal
