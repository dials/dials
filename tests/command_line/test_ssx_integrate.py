from pathlib import Path

import pytest

from dxtbx.serialize import load

from dials.array_family import flex
from dials.command_line.ssx_integrate import run as run_integrate

# Note that tests are grouped and run serially, to stop many processes trying to
# extract data from images at same time, which appears to lead to race
# conditions in some CI jobs.


@pytest.mark.xdist_group(name="group1")
def test_ssx_integrate_stills(dials_data, run_in_tmpdir):
    # Download data set and the internally referenced images
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    dials_data("cunir_serial", pathlib=True)

    indexed_refl = str(ssx / "indexed.refl")
    indexed_expts = str(ssx / "indexed.expt")

    run_integrate(
        [
            indexed_refl,
            indexed_expts,
            "algorithm=stills",
            "nproc=1",
            "json=data.json",
            "image_range=1:2",
        ]
    )

    assert Path("integrated_0.refl").is_file()
    assert Path("integrated_0.expt").is_file()
    assert Path("dials.ssx_integrate.html").is_file()
    assert Path("data.json").is_file()
    experiments = load.experiment_list("integrated_0.expt")
    assert len(experiments) == 2
    refls = flex.reflection_table.from_file("integrated_0.refl")
    expected_n_refls = 614
    assert len(refls) == pytest.approx(expected_n_refls, abs=9)

    run_integrate(
        [
            indexed_refl,
            indexed_expts,
            "algorithm=ellipsoid",
            "nproc=1",
            "image_range=1:2",
        ]
    )

    experiments = load.experiment_list("integrated_0.expt")
    assert len(experiments) == 2
    assert experiments.profiles()[0].name == "ellipsoid"
    refls = flex.reflection_table.from_file("integrated_0.refl")
    expected_n_refls = 1258
    assert len(refls) == pytest.approx(expected_n_refls, abs=9)
    indexed = load.experiment_list(indexed_expts)
    assert indexed[0].crystal != experiments[0].crystal

    # now run with fixing uc and orientation
    run_integrate(
        [
            indexed_refl,
            indexed_expts,
            "algorithm=ellipsoid",
            "nproc=1",
            "image_range=1:2",
            "unit_cell.fixed=True",
            "orientation.fixed=True",
        ]
    )

    experiments = load.experiment_list("integrated_0.expt")
    assert len(experiments) == 2
    assert experiments.profiles()[0].name == "ellipsoid"
    refls = flex.reflection_table.from_file("integrated_0.refl")
    expected_n_refls = 1266
    assert len(refls) == pytest.approx(expected_n_refls, abs=9)
    assert indexed[0].crystal == experiments[0].crystal
