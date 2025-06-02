from __future__ import annotations

import shutil
import subprocess

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.phil import parse

from dials.array_family import flex


def test_kapton(tmp_path, dials_data):
    """Test script for kapton correction applied to integrated data.
    Currently only testing kapton 2019 correction on rayonix-340 at LCLS
    """
    dd = dials_data("lcls_rayonix_kapton", pathlib=True)
    image_file = dd / "hit-20181213155134902.cbf"
    mask_file = dd / "mask_rayonix340mx_4x4.pickle"
    geom_file = dd / "experiments_000.json"

    # Create phil files for the two situations being tests
    #  a. without kapton
    #  b. with kapton

    stills_process_input = parse(
        f"""spotfinder.lookup.mask={mask_file}\n
                          integration.lookup.mask={mask_file}\n
                          input.reference_geometry={geom_file}\n
                          spotfinder.filter.min_spot_size=2\n
                          spotfinder.filter.d_min=2\n
                          spotfinder.filter.d_max=18\n
                          spotfinder.threshold.dispersion.gain=0.46\n
                          spotfinder.threshold.dispersion.global_threshold=100\n
                          indexing.known_symmetry.space_group='P 21 21 21'\n
                          indexing.known_symmetry.unit_cell='41.9 75.7 102 90 90 90'\n
                          indexing.refinement_protocol.d_min_start=2\n
                          integration.debug.output=True\n
                          integration.debug.separate_files=False\n
                          integration.debug.delete_shoeboxes=True\n
                          profile.gaussian_rs.centroid_definition=com\n """
    )

    kapton_input = parse(
        """ integration {
                       absorption_correction {
                         apply=True
                         algorithm=kapton_2019
                         fuller_kapton {
                           xtal_height_above_kapton_mm {
                               value=0.04
                             }
                           rotation_angle_deg {
                             value=0.55
                             }
                           kapton_half_width_mm {
                             value=0.665
                             }
                           kapton_thickness_mm {
                             value=0.025
                             }
                           smart_sigmas=True
                           }
                         }
                       }"""
    )

    with_kapton_phil = tmp_path / "params_with_kapton.phil"
    without_kapton_phil = tmp_path / "params_without_kapton.phil"

    without_kapton_phil.write_text(
        stills_process_input.as_str()
        + "output.integrated_filename=without_kapton.mpack\noutput.integrated_experiments_filename=without_kapton.expt"
    )
    with_kapton_phil.write_text(
        stills_process_input.as_str()
        + kapton_input.as_str()
        + "output.integrated_filename=with_kapton.mpack\noutput.integrated_experiments_filename=with_kapton.expt"
    )

    command_without_kapton = (
        shutil.which("dials.stills_process"),
        image_file,
        "params_without_kapton.phil",
    )
    command_with_kapton = (
        shutil.which("dials.stills_process"),
        image_file,
        "params_with_kapton.phil",
    )
    subprocess.run(command_without_kapton, cwd=tmp_path, capture_output=True)
    subprocess.run(command_with_kapton, cwd=tmp_path, capture_output=True)

    # Now compare the 2 experimental results
    # Currently just comparing the median values to get a sense of the effect if the kapton and whether it is being applied correctly
    expt_without_kapton = ExperimentListFactory.from_json_file(
        tmp_path / "without_kapton.expt", check_format=False
    )
    refl_without_kapton = flex.reflection_table.from_file(
        tmp_path / "without_kapton.mpack"
    )
    expt_with_kapton = ExperimentListFactory.from_json_file(
        tmp_path / "with_kapton.expt", check_format=False
    )
    refl_with_kapton = flex.reflection_table.from_file(tmp_path / "with_kapton.mpack")

    without_kapton_medians = []
    with_kapton_medians = []
    count = 0
    for experiments, reflections in zip(
        (expt_without_kapton, expt_with_kapton), (refl_without_kapton, refl_with_kapton)
    ):
        all_x, all_y, all_i = flex.double(), flex.double(), flex.double()
        for expt_id, experiment in enumerate(experiments):
            refls = reflections.select(reflections["id"] == expt_id)
            for panel_id, panel in enumerate(experiment.detector):
                panel_refls = refls.select(refls["panel"] == panel_id)
                x, y, z = panel_refls["xyzobs.px.value"].parts()
                for i in range(len(panel_refls)):
                    lab_x, lab_y, lab_z = panel.get_pixel_lab_coord((x[i], y[i]))
                    all_x.append(lab_x)
                    all_y.append(lab_y)
                    all_i.append(panel_refls["intensity.sum.value"][i])

        for sel in all_x <= 0, all_x > 0, all_y <= 0, all_y > 0:
            if count == 0:
                without_kapton_medians.append(flex.median(all_i.select(sel)))
            if count == 1:
                with_kapton_medians.append(flex.median(all_i.select(sel)))
        count += 1

    # Now compare results between uncorrected and corrected data

    # x < 0 where the kapton shadow is
    assert without_kapton_medians[0] < with_kapton_medians[0]
    # x > 0 where no kapton shadow present
    assert without_kapton_medians[1] == pytest.approx(with_kapton_medians[1], abs=0.1)
    # y < 0; kapton correction should average out but should be slightly higher
    assert without_kapton_medians[2] == pytest.approx(with_kapton_medians[2], abs=5.0)
    assert without_kapton_medians[2] < with_kapton_medians[2]
    # y < 0; kapton correction should average out but should be slightly higher
    assert without_kapton_medians[3] == pytest.approx(with_kapton_medians[3], abs=5.0)
    assert without_kapton_medians[3] < with_kapton_medians[3]
