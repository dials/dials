from __future__ import annotations

import logging
import re
from math import isclose
from typing import Tuple

from cctbx import crystal, miller
from cctbx.eltbx import chemical_elements
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

    # and a stub of an .ins file with information from the .expt file
    _write_ins(
        experiment_list,
        best_unit_cell=params.mtz.best_unit_cell,
        composition=params.shelx.composition,
        ins_file=params.shelx.ins,
    )

    logger.info(f"Written {params.shelx.ins}")

    # If requested scale up data to maximise use of scale_range
    scale = None
    sigmas = flex.sqrt(variances)
    if params.shelx.scale:
        max_val = max(max(sigmas), max(intensities))
        min_val = min(min(sigmas), min(intensities))

        min_scale = params.shelx.scale_range[0] / min_val
        max_scale = params.shelx.scale_range[1] / max_val
        scale = min(min_scale, max_scale)

    if experiment_list.all_laue() or experiment_list.all_tof():
        sigmas = flex.sqrt(variances)
        if scale is not None:
            intensities *= scale
            sigmas *= scale
        batch_numbers = _generate_laue_batch_numbers(scaled_data)
        return _export_as_shelx_hklf2(
            filename=params.shelx.hklout,
            miller_indices=scaled_data["miller_index"],
            intensities=intensities,
            sigmas=sigmas,
            wavelengths=scaled_data["wavelength_cal"],
            batch_numbers=batch_numbers,
            normalise_if_format_overflow=True,
            scale_range=params.shelx.scale_range,
            scale_factor=scale,
        )

    i_obs = miller.array(miller_set, data=intensities)
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(flex.sqrt(variances))
    if params.shelx.scale:
        i_obs = i_obs.apply_scaling(factor=scale)

    # write the SHELX HKL file
    hkl_file = params.shelx.hklout
    with open(hkl_file, "w") as hkl_file_object:
        i_obs.export_as_shelx_hklf(
            file_object=hkl_file_object,
            scale_range=params.shelx.scale_range,
            normalise_if_format_overflow=True,
        )
    logger.info(f"Written {i_obs.size()} reflections to {hkl_file}")


def _generate_laue_batch_numbers(reflections: flex.reflection_table) -> flex.int:
    if "imageset_id" in reflections:
        idx_id = reflections["imageset_id"]
    else:
        idx_id = reflections["id"]

    return idx_id


def _export_as_shelx_hklf2(
    filename: str,
    miller_indices: flex.miller_index,
    intensities: flex.double,
    sigmas: flex.double,
    batch_numbers: flex.int,
    wavelengths: flex.double,
    normalise_if_format_overflow: bool = False,
    full_dynamic_range: bool = False,
    scale_range: Tuple[float, float] | None = None,
    scale_factor: float | None = None,
) -> None:
    """
    Writes an HKL file using the HKLF 2 format:
    H | K | L | Fo2 | σ(Fo2) | Batch | λ
    """

    ## Formatters to ensure .hkl file has correct columns

    def fmt_3i4(h: int, k: int, l: int) -> str:
        result = f"{int(h):4d}{int(k):4d}{int(l):4d}"
        if len(result) != 12:
            raise RuntimeError(f"SHELXL hkl file 3I4 format overflow: {result}")
        return result

    def fmt_f8(v: float) -> str:
        result = f"{v:8.2f}"
        if len(result) != 8:
            result = f"{v:8.1f}"
            if len(result) != 8:
                result = f"{round(v):7d}."
        return result

    def fmt_fullrange_data(v: float) -> str:
        if v < 0.0:
            if abs(v) < 1.0:
                result = f"{v:8.5f}"
            else:
                result = f"{v:8.6g}"
                if "." not in result:
                    result = f"{round(v):7d}."
        else:
            if v < 1.0:
                result = f"{v:8.6f}"
            else:
                result = f"{v:8.7g}"
            if "." not in result:
                result = f"{round(v):7d}."
        return result

    def fmt_fullrange_sigma(v: float) -> str:
        if abs(v) >= 1.0:
            result = f"{v:8.7g}"
            if "." not in result:
                result = f"{round(v):7d}."
        elif abs(v) < 1e-5:
            result = f"{0.00001:8.5g}"
        else:
            result = f"{v:8.6g}"
        return result

    ## Sanity check

    assert len(intensities) == len(sigmas)
    assert len(intensities) == len(miller_indices)
    assert len(intensities) == len(batch_numbers)
    assert len(intensities) == len(wavelengths)

    ## Scaling (same approach as cctbx miller_array_export_as_shelx_hklf)
    # https://github.com/cctbx/cctbx_project/blob/09bc17a40132e4204b1142a093baed33929f6212/iotbx/shelx/hklf.py#L7
    if scale_factor is not None:
        intensities *= scale_factor
        sigmas *= scale_factor

    min_val = min(min(intensities), min(sigmas))
    max_val = max(max(intensities), max(sigmas))
    min_sc = 1.0
    max_sc = 1.0
    scale = 1.0

    if full_dynamic_range:
        # Use same values as miller_array_export_as_shelx_hklf
        max_abs = 999999.0
        max_val = max(abs(max_val), abs(min_val))
        if max_val > max_abs:
            scale = max_abs / max_val
    else:
        if scale_range is None:
            # Use same values as miller_array_export_as_shelx_hklf
            scale_range = (-999999.0, 9999999.0)
        if min_val < scale_range[0]:
            if not normalise_if_format_overflow:
                raise RuntimeError(f"SHELX hkl file F8.2 format overflow: {min_val}")
            min_sc = scale_range[0] / min_val
        if max_val > scale_range[1]:
            if not normalise_if_format_overflow:
                raise RuntimeError(f"SHELX hkl file F8.2 format overflow: {max_val}")
            max_sc = scale_range[1] / max_val
        scale = min(min_sc, max_sc)

    ## Write file

    with open(filename, "w") as g:
        for hkl, I, s, b, wl in zip(
            miller_indices, intensities, sigmas, batch_numbers, wavelengths
        ):
            h, k, l = hkl
            I_scaled = I * scale
            s_scaled = s * scale

            if full_dynamic_range:
                line = (
                    fmt_3i4(h, k, l)
                    + fmt_fullrange_data(I_scaled)
                    + fmt_fullrange_sigma(s_scaled)
                )
            else:
                line = fmt_3i4(h, k, l) + fmt_f8(I_scaled) + fmt_f8(s_scaled)

            line += f"{int(b):4d}{float(wl):8.4f}\n"
            g.write(line)

        g.write(f"{0:4d}{0:4d}{0:4d}{0.0:8.2f}{0.0:8.2f}{0:4d}{0.0:8.4f}\n")


def _write_ins(experiment_list, best_unit_cell, composition, ins_file):
    sg = experiment_list[0].crystal.get_space_group()
    unit_cells = []
    wavelengths = []

    if experiment_list.all_laue() or experiment_list.all_tof():
        wl = 1.0
    else:
        # Check for single wavelength
        for exp in experiment_list:
            wl = exp.beam.get_wavelength()
            if not any(isclose(wl, w, abs_tol=1e-4) for w in wavelengths):
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
            f"TITL {sg.type().number()} in {sg.type().lookup_symbol().replace(' ', '')}\n"
        )
        f.write(
            "CELL {:7.5f} {:9.5f} {:9.5f} {:9.5f} {:8.4f} {:8.4f} {:8.4f}\n".format(
                wl, *uc.parameters()
            )
        )
        if uc_sd:
            f.write(
                "ZERR {:7.2f} {:9.5f} {:9.5f} {:9.5f} {:8.4f} {:8.4f} {:8.4f}\n".format(
                    sg.order_z(), *uc_sd
                )
            )
        LATT_SYMM(f, sg)
        logger.info(
            f"Using composition {composition} to write SFAC & UNIT instructions"
        )
        sorted_composition = _parse_compound(composition)
        f.write("SFAC %s\n" % " ".join(sorted_composition))
        f.write(
            "UNIT %s\n"
            % " ".join(
                [str(sorted_composition[e] * sg.order_z()) for e in sorted_composition]
            )
        )
        f.write("HKLF 4\n")  # Intensities
        f.write("END\n")


def _parse_compound(composition):
    elements = chemical_elements.proper_caps_list()[:94]
    m = re.findall(r"([A-Z][a-z]?)(\d*)", composition)
    if all(e[1] == "" for e in m):
        # Assume user has just supplied list of elements so UNIT instruction cannot be calculated
        result = {element: 0 for element, count in m}
    else:
        result = {element: int(count or 1) for element, count in m}
    # Check for elements without scattering factors in SHELXL
    unmatched = list(set(result).difference(elements))
    if unmatched:
        raise Sorry(f"Element(s) {' '.join(unmatched)} in {composition} not recognised")

    if "C" in result:
        # move C to start of list
        elements.insert(0, elements.pop(5))
    sorted_result = {e: result[e] for e in elements if e in result}
    return sorted_result
