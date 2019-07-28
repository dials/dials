from __future__ import absolute_import, division, print_function

import os

import pytest

import dxtbx
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.format.FormatCBFCspad import FormatCBFCspadInMemory
from dxtbx.imageset import ImageSet, ImageSetData, MemReader
from libtbx import easy_run
from libtbx.phil import parse

from dials.command_line.stills_process import phil_scope, Processor
from dials.array_family import flex


def test_cspad_cbf_in_memory(dials_regression, run_in_tmpdir):
    # Check the data files for this test exist
    image_path = os.path.join(
        dials_regression,
        "image_examples",
        "LCLS_cspad_nexus",
        "idx-20130301060858801.cbf",
    )
    assert os.path.isfile(image_path)

    with open("process_lcls.phil", "w") as f:
        f.write(
            """
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
        )

    params = phil_scope.fetch(parse(file_name="process_lcls.phil")).extract()
    params.output.experiments_filename = None
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

    result = "idx-20130301060858801_integrated.refl"
    n_refls = list(
        range(140, 152)
    )  # large ranges to handle platform-specific differences
    table = flex.reflection_table.from_file(result)
    assert len(table) in n_refls, len(table)
    assert "id" in table
    assert (table["id"] == 0).count(False) == 0


@pytest.mark.parametrize("use_mpi", [True, False])
def test_sacla_h5(dials_regression, run_in_tmpdir, use_mpi, in_memory=False):
    # Only allow MPI tests if we've got MPI capabilities
    if use_mpi:
        pytest.importorskip("mpi4py")

    # Check the data files for this test exist
    sacla_path = os.path.join(dials_regression, "image_examples", "SACLA_MPCCD_Cheetah")
    image_path = os.path.join(sacla_path, "run266702-0-subset.h5")
    assert os.path.isfile(image_path)

    geometry_path = os.path.join(sacla_path, "refined_experiments_level1.json")
    assert os.path.isfile(geometry_path)

    # Write the .phil configuration to a file
    with open("process_sacla.phil", "w") as f:
        f.write(
            """
      dispatch.squash_errors = True
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
        threshold {
          dispersion {
            gain = 5.46 # from dials.estimate_gain run266702-0-subset.h5 max_images=4
            global_threshold = 50
          }
        }
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
      """
            % geometry_path
        )

    # Call dials.stills_process
    if use_mpi:
        command = ["mpirun", "-n", "4", "dials.stills_process", "mp.method=mpi"]
    else:
        command = ["dials.stills_process"]
    command += [image_path, "process_sacla.phil"]
    result = easy_run.fully_buffered(command).raise_if_errors()
    result.show_stdout()

    for result_filename, n_refls in zip(
        [
            "idx-run266702-0-subset_00000_integrated.refl",
            "idx-run266702-0-subset_00001_integrated.refl",
            "idx-run266702-0-subset_00003_integrated.refl",
        ],
        [list(range(205, 225)), list(range(565, 580)), list(range(475, 500))],
    ):  # large ranges to handle platform-specific differences
        table = flex.reflection_table.from_file(result_filename)
        assert len(table) in n_refls, (result_filename, len(table))
        assert "id" in table
        assert (table["id"] == 0).count(False) == 0
