from __future__ import annotations

import os
from pathlib import Path

import procrunner
import pytest

import dxtbx
from dxtbx.format.FormatCBFCspad import FormatCBFCspadInMemory
from dxtbx.imageset import ImageSet, ImageSetData, MemReader
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.phil import parse

from dials.array_family import flex
from dials.command_line.stills_process import Processor, phil_scope

cspad_cbf_in_memory_phil = """
dispatch.squash_errors = False
spotfinder {
  filter.min_spot_size=2
  threshold.dispersion.gain=25
  threshold.dispersion.global_threshold=100
}
indexing {
  known_symmetry {
    space_group = P6122
    unit_cell = 92.9 92.9 130.4 90 90 120
  }
  refinement_protocol.d_min_start=1.7
  stills.refine_candidates_with_known_symmetry=True
}
"""

sacla_phil = """
dispatch.squash_errors = True
dispatch.coset = True
input.reference_geometry=%s
indexing {
  known_symmetry {
    space_group = P43212
    unit_cell = 78.9 78.9 38.1 90 90 90
  }
  refinement_protocol.d_min_start = 2.2
  stills.refine_candidates_with_known_symmetry=True
}
spotfinder {
  filter.min_spot_size = 2
}
refinement {
  parameterisation {
    detector.fix_list = Dist,Tau1
  }
}
profile {
  gaussian_rs {
    centroid_definition = com
  }
}
output.composite_output = True
"""


@pytest.mark.parametrize("composite_output", [True, False])
def test_cspad_cbf_in_memory(dials_regression, tmp_path, composite_output):
    # Check the data files for this test exist
    image_path = Path(
        dials_regression,
        "image_examples",
        "LCLS_cspad_nexus",
        "idx-20130301060858801.cbf",
    )
    assert image_path.is_file()

    tmp_path.joinpath("process_lcls.phil").write_text(cspad_cbf_in_memory_phil)

    params = phil_scope.fetch(parse(file_name=tmp_path / "process_lcls.phil")).extract()
    params.output.experiments_filename = None
    params.output.composite_output = composite_output
    cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        if composite_output:
            processor = Processor(params, composite_tag="memtest")
        else:
            processor = Processor(params)
        mem_img = dxtbx.load(image_path)
        raw_data = mem_img.get_raw_data()  # cache the raw data to prevent swig errors
        mem_img = FormatCBFCspadInMemory(mem_img._cbf_handle)
        mem_img._raw_data = raw_data
        mem_img._cbf_handle = None  # drop the file handle to prevent swig errors
        imgset = ImageSet(ImageSetData(MemReader([mem_img]), None))
        imgset.set_beam(mem_img.get_beam())
        imgset.set_detector(mem_img.get_detector())
        experiments = ExperimentListFactory.from_imageset_and_crystal(imgset, None)
        processor.process_experiments(
            "20130301060858801", experiments
        )  # index/integrate the image
        if composite_output:
            processor.finalize()
            result = "idx-memtest_integrated.refl"
        else:
            result = "idx-20130301060858801_integrated.refl"
    finally:
        os.chdir(cwd)

    n_refls = list(
        range(140, 152)
    )  # large ranges to handle platform-specific differences
    table = flex.reflection_table.from_file(tmp_path / result)
    assert len(table) in n_refls, len(table)
    assert "id" in table
    assert (table["id"] == 0).count(False) == 0


@pytest.mark.parametrize("use_mpi", [True, False])
def test_sacla_h5(dials_data, tmp_path, use_mpi, in_memory=False):
    # Only allow MPI tests if we've got MPI capabilities
    if use_mpi:
        pytest.importorskip("mpi4py")

    # Check the data files for this test exist
    sacla_path = dials_data("image_examples", pathlib=True)
    image_path = sacla_path / "SACLA-MPCCD-run266702-0-subset.h5"
    assert image_path.is_file()

    geometry_path = (
        sacla_path / "SACLA-MPCCD-run266702-0-subset-refined_experiments_level1.json"
    )
    assert geometry_path.is_file()

    # Write the .phil configuration to a file
    tmp_path.joinpath("process_sacla.phil").write_text(sacla_phil % geometry_path)

    # Call dials.stills_process
    if use_mpi:
        command = [
            "mpirun",
            "-n",
            "4",
            "dials.stills_process",
            "mp.method=mpi mp.composite_stride=4 output.logging_dir=.",
        ]
    else:
        command = ["dials.stills_process"]
    command += [image_path, "process_sacla.phil"]
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr

    def test_refl_table(result_filename, ranges):
        table = flex.reflection_table.from_file(result_filename)
        for expt_id, n_refls in enumerate(ranges):
            subset = table.select(table["id"] == expt_id)
            assert len(subset) in n_refls, (result_filename, expt_id, len(table))
        assert "id" in table
        assert set(table["id"]) == {0, 1, 2, 3}

    # large ranges to handle platform-specific differences
    test_refl_table(
        tmp_path / "idx-0000_integrated.refl",
        [
            list(range(140, 160)),
            list(range(575, 600)),
            list(range(420, 445)),
            list(range(485, 510)),
        ],
    )

    test_refl_table(
        tmp_path / "idx-0000_coset6.refl",
        [
            list(range(145, 160)),
            list(range(545, 570)),
            list(range(430, 455)),
            list(range(490, 515)),
        ],
    )


def test_pseudo_scan(dials_data, tmp_path):
    result = procrunner.run(
        (
            "dials.stills_process",
            dials_data("centroid_test_data", pathlib=True) / "centroid_000[1-2].cbf",
            "convert_sequences_to_stills=True",
            "squash_errors=False",
            "composite_output=True",
        ),
        working_directory=tmp_path,
    )
    assert not result.returncode and not result.stderr

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "idx-0000_refined.expt", check_format=False
    )
    assert len(experiments) == 2
