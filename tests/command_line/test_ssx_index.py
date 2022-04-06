from __future__ import annotations

import os.path

import procrunner

from dxtbx.serialize import load

from dials.command_line.ssx_index import run


def test_ssx_index_reference_geometry(dials_data, tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = ssx / "imported_with_ref_5.expt"
    refls = ssx / "strong_5.refl"

    result = procrunner.run(
        ["dev.dials.ssx_index", expts, refls],
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    assert (tmp_path / "dials.ssx_index.html").is_file()
    experiments = load.experiment_list(tmp_path / "indexed.expt", check_format=False)
    assert len(experiments) == 4  # only 4 out of the 5 get indexed


def test_ssx_index_no_reference_geometry(dials_data, tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = ssx / "imported_no_ref_5.expt"
    refls = ssx / "strong_5.refl"

    result = procrunner.run(
        ["dev.dials.ssx_index", expts, refls, "-vv"],
        working_directory=tmp_path,
    )

    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    assert (tmp_path / "dials.ssx_index.html").is_file()
    experiments = load.experiment_list(tmp_path / "indexed.expt", check_format=False)
    assert (
        len(experiments) == 3
    )  # only 3 out of the 5 get indexed if no reference geometry


def test_ssx_index_bad_input(dials_data, run_in_tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = str(ssx / "imported_no_ref_5.expt")
    refls = str(ssx / "strong_1.refl")

    run([expts, refls])
    assert os.path.exists("indexed.refl")
    assert os.path.exists("indexed.expt")
    experiments = load.experiment_list("indexed.expt", check_format=False)
    assert len(experiments) == 0


def test_ssx_index_input_unit_cell(dials_data, run_in_tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = str(ssx / "imported_with_ref_5.expt")
    refls = str(ssx / "strong_5.refl")

    # invoke the run function
    run([expts, refls, "max_lattices=2", "unit_cell=96.4,96.4,96.4,90,90,90"])

    assert os.path.exists("indexed.refl")
    assert os.path.exists("indexed.expt")
    assert os.path.exists("dials.ssx_index.html")
    experiments = load.experiment_list("indexed.expt", check_format=False)
    assert (
        len(experiments) == 5
    )  # only 4 out of the 5 images gets indexed, one has two lattices

    # now if just fft1d method:
    # invoke the run function
    run([expts, refls, "method=fft1d", "unit_cell=96.4,96.4,96.4,90,90,90"])
    experiments = load.experiment_list("indexed.expt", check_format=False)
    assert (
        len(experiments) == 2
    )  # only 2 out of the 5 images get indexed without real space grid search
