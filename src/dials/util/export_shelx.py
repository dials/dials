from __future__ import annotations

import logging
from math import isclose

from cctbx import crystal, miller
from iotbx.shelx.write_ins import LATT_SYMM

from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.array_family import flex
from dials.util import Sorry
from dials.util.filter_reflections import filter_reflection_table

logger = logging.getLogger(__name__)


def export_shelx(scaled_data, experiment_list, params):
    """Export scaled data corresponding to experiment_list to
    a SHELX HKL formatted text file."""

    # Handle requesting profile intensities (default via auto) but no column
    if "profile" in params.intensity and "intensity.prf.value" not in scaled_data:
        raise Sorry(
            "Requested profile intensity data but only summed present. Use intensity=sum."
        )

    # use supplied best unit cell or that determined from experiment list to define d in reflection table.
    best_unit_cell = params.mtz.best_unit_cell
    if best_unit_cell is None:
        best_unit_cell = determine_best_unit_cell(experiment_list)
    else:
        logger.info("Using supplied unit cell across experiments : %s", best_unit_cell)
    scaled_data["d"] = best_unit_cell.d(scaled_data["miller_index"])

    # Clean up reflection table with mtz defaults (as in export_xds_ascii)
    scaled_data = filter_reflection_table(
        scaled_data,
        intensity_choice=params.intensity,
        partiality_threshold=params.mtz.partiality_threshold,
        combine_partials=params.mtz.combine_partials,
        min_isigi=params.mtz.min_isigi,
        filter_ice_rings=params.mtz.filter_ice_rings,
        d_min=params.mtz.d_min,
    )

    # Check that all experiments have the same space group
    if len({x.crystal.get_space_group().make_tidy() for x in experiment_list}) != 1:
        raise ValueError("Experiments do not have a unique space group")

    # Create miller set with space group from 1st crystal in experiment list and best unit cell
    miller_set = miller.set(
        crystal_symmetry=crystal.symmetry(
            unit_cell=best_unit_cell,
            space_group=experiment_list[0].crystal.get_space_group(),
        ),
        indices=scaled_data["miller_index"],
        anomalous_flag=False,
    )

    intensity_choice = (
        params.intensity[0] if params.intensity[0] != "profile" else "prf"
    )
    intensities, variances = (
        scaled_data["intensity." + intensity_choice + ".value"],
        scaled_data["intensity." + intensity_choice + ".variance"],
    )
    assert variances.all_gt(0)
    i_obs = miller.array(miller_set, data=intensities)
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(flex.sqrt(variances))

    # If requested scale up data to maximise use of scale_range
    if params.shelx.scale:
        max_val = max(
            i_obs.data().min_max_mean().max, i_obs.sigmas().min_max_mean().max
        )
        min_val = min(
            i_obs.data().min_max_mean().min, i_obs.sigmas().min_max_mean().min
        )
        min_scale = params.shelx.scale_range[0] / min_val
        max_scale = params.shelx.scale_range[1] / max_val
        scale = min(min_scale, max_scale)
        i_obs = i_obs.apply_scaling(factor=scale)

    # write the SHELX HKL file
    hkl_file = params.shelx.hklout
    with open(hkl_file, "w") as hkl_file_object:
        i_obs.export_as_shelx_hklf(
            file_object=hkl_file_object,
            scale_range=params.shelx.scale_range,
            normalise_if_format_overflow=True,
        )
    logger.info(f"Written {i_obs.size()} relflections to {hkl_file}")

    # and a stub of an .ins file with information from the .expt file
    _write_ins(
        experiment_list,
        best_unit_cell=params.mtz.best_unit_cell,
        ins_file=params.shelx.ins,
    )
    logger.info(f"Written {params.shelx.ins}")


def _write_ins(experiment_list, best_unit_cell, ins_file):
    sg = experiment_list[0].crystal.get_space_group()
    unit_cells = []
    wavelengths = []

    # Check for single wavelength
    for exp in experiment_list:
        wl = exp.beam.get_wavelength()
        if not any([isclose(wl, w, abs_tol=1e-4) for w in wavelengths]):
            wavelengths.append(wl)
    if len(wavelengths) > 1:
        raise ValueError("Experiments have more than one wavelength")
    else:
        wl = wavelengths[0]

    # if user has supplied best_unit_cell use it
    if best_unit_cell is not None:
        uc = best_unit_cell
        uc_sd = None
    else:
        for exp in experiment_list:
            unit_cells.append(
                exp.crystal.get_recalculated_unit_cell() or exp.crystal.get_unit_cell()
            )

        if len(unit_cells) > 1:
            if (
                len({uc.parameters() for uc in unit_cells}) > 1
            ):  # have different cells so no esds
                uc = determine_best_unit_cell(experiment_list)
                uc_sd = None
            else:  # identical (recalculated?) unit cell with esds
                uc = (
                    experiment_list[0].crystal.get_recalculated_unit_cell()
                    or experiment_list[0].crystal.get_unit_cell()
                )
                uc_sd = (
                    experiment_list[0].crystal.get_recalculated_cell_parameter_sd()
                    or experiment_list[0].crystal.get_cell_parameter_sd()
                )
        else:  # single unit cell
            uc = (
                experiment_list[0].crystal.get_recalculated_unit_cell()
                or experiment_list[0].crystal.get_unit_cell()
            )
            uc_sd = (
                experiment_list[0].crystal.get_recalculated_cell_parameter_sd()
                or experiment_list[0].crystal.get_cell_parameter_sd()
            )

    with open(ins_file, "w") as f:
        f.write(
            f"TITL {sg.type().number()} in {sg.type().lookup_symbol().replace(' ','')}\n"
        )
        f.write(
            "CELL {:8.5f} {:8.4f} {:8.4f} {:8.4f} {:8.3f} {:8.3f} {:8.3f}\n".format(
                wl, *uc.parameters()
            )
        )
        if uc_sd:
            f.write(
                "ZERR {:8.3f} {:8.4f} {:8.4f} {:8.4f} {:8.3f} {:8.3f} {:8.3f}\n".format(
                    sg.order_z(), *uc_sd
                )
            )
        LATT_SYMM(f, sg)
