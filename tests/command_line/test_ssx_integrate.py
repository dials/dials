import os

import pytest

from dxtbx.serialize import load

from dials.array_family import flex
from dials.command_line.ssx_integrate import run as run_integrate

# Note that tests are grouped and run serially, to stop many processes trying to
# extract data from images at same time, which appears to lead to race
# conditions in CI jobs.


@pytest.fixture
def indexed_data(dials_data):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    _ = dials_data("cunir_serial", pathlib=True)  # Make sure the images are downloaded
    return (str(ssx / "indexed.refl"), str(ssx / "indexed.expt"))


@pytest.mark.xdist_group(name="group1")
@pytest.mark.serial
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
@pytest.mark.xdist_group(name="group1")
@pytest.mark.serial
def test_ssx_integrate_ellipsoid(indexed_data, run_in_tmpdir, fix_uc_and_orientation):
    refls, expts = indexed_data

    # set batch size to test generation of multiple output files
    args = [refls, expts, "algorithm=ellipsoid", "batch_size=3", "nproc=1", "-vv"]
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
    refls_1 = flex.reflection_table.from_file("integrated_1.refl")
    print(len(refls_0))
    print(len(refls_1))
    assert (
        len(refls_0) > 3740 - 20 and len(refls_0) < 3740 + 20
    )  # 3730/3740 at 18/01/22
    experiments_1 = load.experiment_list("integrated_1.expt")
    assert len(experiments_1) == 2

    assert (
        len(refls_1) > 1470 - 20 and len(refls_1) < 1470 + 20
    )  # 1467/1473 at 18/01/22
    indexed = load.experiment_list(expts)
    if fix_uc_and_orientation:
        assert indexed[0].crystal == experiments_0[0].crystal
    else:
        assert indexed[0].crystal != experiments_0[0].crystal


@pytest.mark.parametrize(
    "rlp_mosaicity,expected_results",
    [
        (
            "simple1",
            {
                "n_refl": 13376,
                "mosaicity": {"spherical": [0.0239, 0.0214, 0.0341, 0.0209, 0.0214]},
            },
        ),
        (
            "simple6",
            {
                "n_refl": 7245,
                "mosaicity": {
                    "min": [0.0096, 0.0075, 0.0222, 0.0079, 0.0107],
                    "mid": [0.0240, 0.0227, 0.0357, 0.0227, 0.0228],
                    "max": [0.0323, 0.0283, 0.0414, 0.0269, 0.0271],
                },
            },
        ),
        (
            "angular2",
            {
                "n_refl": 12548,
                "mosaicity": {
                    "radial": [0.0285, 0.0241, 0.0347, 0.0246, 0.0247],
                    "angular": [0.0211, 0.0199, 0.0337, 0.0187, 0.0195],
                },
            },
        ),
        (
            "angular4",
            {
                "n_refl": 5203,
                "mosaicity": {
                    "radial": [0.0285, 0.0241, 0.0348, 0.0246, 0.0247],
                    "angular_0": [0.0292, 0.0278, 0.0426, 0.0259, 0.0263],
                    "angular_1": [0.0066, 0.0047, 0.0215, 0.0053, 0.0083],
                },
            },
        ),
    ],
)
@pytest.mark.xdist_group(name="group1")
@pytest.mark.serial
def test_ssx_integrate_ellipsoid_profile_models(
    indexed_data, run_in_tmpdir, rlp_mosaicity, expected_results
):
    refls, expts = indexed_data

    # set batch size to test generation of multiple output files
    run_integrate(
        [
            refls,
            expts,
            "algorithm=ellipsoid",
            "nproc=1",
            f"rlp_mosaicity={rlp_mosaicity}",
            "-vv",
        ]
    )

    assert os.path.exists("integrated_0.refl")
    assert os.path.exists("integrated_0.expt")
    assert os.path.exists("dials.ssx_integrate.html")

    refls_0 = flex.reflection_table.from_file("integrated_0.refl")

    assert len(refls_0) > expected_results["n_refl"] - 20
    assert len(refls_0) < expected_results["n_refl"] + 20

    expts = load.experiment_list("integrated_0.expt")
    for i, p in enumerate(expts.profiles()):
        m = p.mosaicity()
        expected_mosaicity = expected_results["mosaicity"]
        expected = {k: expected_mosaicity[k][i] for k in expected_mosaicity.keys()}
        for k, v in expected.items():
            print(v, m[k], k)
            assert v == pytest.approx(m[k], abs=1e-4)
