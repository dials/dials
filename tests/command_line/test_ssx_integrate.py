import os

import pytest

from dxtbx.serialize import load

from dials.array_family import flex
from dials.command_line.ssx_integrate import run as run_integrate


@pytest.fixture
def indexed_data(dials_data):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    _ = dials_data("cunir_serial", pathlib=True)  # Make sure the images are downloaded
    return (str(ssx / "indexed.refl"), str(ssx / "indexed.expt"))


def test_ssx_integrate_stills(indexed_data, run_in_tmpdir):
    refls, expts = indexed_data

    run_integrate([refls, expts, "algorithm=stills", "nproc=1", "json=data.json"])

    assert os.path.exists("integrated_0.refl")
    assert os.path.exists("integrated_0.expt")
    assert os.path.exists("dials.ssx_integrate.html")
    assert os.path.exists("data.json")
    experiments = load.experiment_list("integrated_0.expt")
    assert len(experiments) == 5
    refls = flex.reflection_table.from_file("integrated_0.refl")
    assert len(refls) > 2150 and len(refls) < 2200  # 2179 at 18/01/22


@pytest.mark.parametrize("fix_uc_and_orientation", [False, True])
def test_ssx_integrate_potato(indexed_data, run_in_tmpdir, fix_uc_and_orientation):
    refls, expts = indexed_data

    # set batch size to test generation of multiple output files
    args = [refls, expts, "algorithm=potato", "batch_size=3", "nproc=1"]
    if fix_uc_and_orientation:
        args.extend(["unit_cell.fixed=True", "orientation.fixed=True"])
    run_integrate(args)

    assert os.path.exists("integrated_0.refl")
    assert os.path.exists("integrated_0.expt")
    assert os.path.exists("integrated_1.refl")
    assert os.path.exists("integrated_1.expt")
    assert os.path.exists("dials.ssx_integrate.html")
    experiments_0 = load.experiment_list("integrated_0.expt")
    assert len(experiments_0) == 3
    refls_0 = flex.reflection_table.from_file("integrated_0.refl")
    assert len(refls_0) > 3710 and len(refls_0) < 3760  # 3730/3740 at 18/01/22
    experiments_1 = load.experiment_list("integrated_1.expt")
    assert len(experiments_1) == 2
    refls_1 = flex.reflection_table.from_file("integrated_1.refl")
    assert len(refls_1) > 1450 and len(refls_1) < 1500  # 1467/1473 at 18/01/22
    indexed = load.experiment_list(expts)
    if fix_uc_and_orientation:
        assert indexed[0].crystal == experiments_0[0].crystal
    else:
        assert indexed[0].crystal != experiments_0[0].crystal


@pytest.mark.parametrize(
    "rlp_mosaicity", ["simple1", "simple6", "angular2", "angular4"]
)
def test_ssx_integrate_potato_profile_models(
    indexed_data, run_in_tmpdir, rlp_mosaicity
):
    refls, expts = indexed_data

    # set batch size to test generation of multiple output files
    run_integrate(
        [refls, expts, "algorithm=potato", "nproc=1", f"rlp_mosaicity={rlp_mosaicity}"]
    )

    assert os.path.exists("integrated_0.refl")
    assert os.path.exists("integrated_0.expt")
    assert os.path.exists("dials.ssx_integrate.html")
