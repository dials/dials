import logging

from cctbx import crystal, miller

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
