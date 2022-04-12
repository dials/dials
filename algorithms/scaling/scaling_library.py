"""
Module of library functions, to perform core scaling operations on reflection
tables and experiment lists. Some functions, such as create_scaling_model and
merging statistics calculations are called from the main dials.scale script,
whereas others are provided as library functions for calling from custom
scripts. The functions defined here should ideally only require reflection
tables and ExperimentList objects (and sometimes phil_scope objects if
necessary), and return common dials objects such as reflection tables and
ExperimentLists.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from unittest.mock import Mock

import numpy as np
import pkg_resources

import iotbx.merging_statistics
from cctbx import crystal, miller, uctbx
from dxtbx import flumpy
from dxtbx.model import Experiment
from dxtbx.util import ersatz_uuid4
from iotbx import cif, mtz, pdb
from libtbx import Auto, phil

from dials.algorithms.scaling.Ih_table import IhTable
from dials.algorithms.scaling.model.model import KBScalingModel, PhysicalScalingModel
from dials.algorithms.scaling.scaling_utilities import (
    DialsMergingStatisticsError,
    calculate_prescaling_correction,
)
from dials.array_family import flex
from dials.util import Sorry
from dials.util.options import ArgumentParser

logger = logging.getLogger("dials")


def set_image_ranges_in_scaling_models(experiments):
    """Set the batch range in scaling models if not already set."""
    for exp in experiments:
        if exp.scan:
            valid_image_ranges = exp.scan.get_valid_image_ranges(exp.identifier)
            if "valid_image_range" not in exp.scaling_model.configdict:
                # only set if not currently set i.e. set initial
                exp.scaling_model.set_valid_image_range(exp.scan.get_image_range())
            if exp.scaling_model.configdict["valid_image_range"] != [
                valid_image_ranges[0][0],
                valid_image_ranges[-1][1],
            ]:
                # first and last values in whole list of tuples
                exp.scaling_model.limit_image_range(
                    (valid_image_ranges[0][0], valid_image_ranges[-1][1])
                )
    return experiments


def choose_initial_scaling_intensities(reflection_table, intensity_choice="profile"):
    """Choose which intensities to initially use for scaling. The LP, QE and
    partiality corrections are also applied. Two new columns are
    added to the reflection table 'intensity' and 'variance', which have
    all corrections applied except an inverse scale factor."""
    if intensity_choice == "profile":
        intensity_choice = "prf"  # rename to allow string matching with refl table
    if "intensity.prf.value" not in reflection_table:
        intensity_choice = "sum"
    elif intensity_choice == "prf":
        if (
            reflection_table.get_flags(reflection_table.flags.integrated_prf).count(
                True
            )
            == 0
        ):
            logger.warning(
                "No profile fitted reflections in this dataset, using summation intensities"
            )
            intensity_choice = "sum"
    reflection_table = calculate_prescaling_correction(reflection_table)
    conv = reflection_table["prescaling_correction"]
    # if prf/sum, use those. If combine, use prf else sum for each refl.
    if intensity_choice == "prf":
        reflection_table["intensity"] = reflection_table["intensity.prf.value"] * conv
        reflection_table["variance"] = (
            reflection_table["intensity.prf.variance"] * conv * conv
        )
    else:
        # first fill in summation intensities.
        if "partiality" in reflection_table:
            inverse_partiality = flex.double(reflection_table.size(), 1.0)
            nonzero_partiality_sel = reflection_table["partiality"] > 0.0
            good_refl = reflection_table.select(reflection_table["partiality"] > 0.0)
            inverse_partiality.set_selected(
                nonzero_partiality_sel.iselection(), 1.0 / good_refl["partiality"]
            )
            reflection_table["intensity"] = (
                reflection_table["intensity.sum.value"] * conv * inverse_partiality
            )
            reflection_table["variance"] = reflection_table[
                "intensity.sum.variance"
            ] * flex.pow2(conv * inverse_partiality)
            if "partiality.inv.variance" in reflection_table:
                reflection_table["variance"] += (
                    reflection_table["intensity.sum.value"]
                    * conv
                    * reflection_table["partiality.inv.variance"]
                )
        else:
            reflection_table["intensity"] = (
                reflection_table["intensity.sum.value"] * conv
            )
            reflection_table["variance"] = (
                reflection_table["intensity.sum.variance"] * conv * conv
            )
        if intensity_choice == "combine":
            # now overwrite prf if we have it.
            sel = reflection_table.get_flags(reflection_table.flags.integrated_prf)
            isel = sel.iselection()
            Iprf = (reflection_table["intensity.prf.value"] * conv).select(sel)
            Vprf = (reflection_table["intensity.prf.variance"] * conv * conv).select(
                sel
            )
            reflection_table["intensity"].set_selected(isel, Iprf)
            reflection_table["variance"].set_selected(isel, Vprf)
    variance_mask = reflection_table["variance"] <= 0.0
    reflection_table.set_flags(
        variance_mask, reflection_table.flags.excluded_for_scaling
    )
    return reflection_table


def scale_against_target(
    reflection_table,
    experiment,
    target_reflection_table,
    target_experiment,
    params=None,
    model="KB",
):
    """Determine scale factors for a single dataset, by scaling against a target
    reflection table. Requires a single reflection table for the reflections to
    scale and the target dataset, and an ExperimentList for both datasets. The
    params option can also be specified, if None then the default scaling
    configuration is used. The scaling model can be specified individually.

    Returns the reflection table, with added columns 'inverse_scale_factor' and
    'inverse_scale_factor_variance'."""

    if not params:
        phil_scope = phil.parse(
            """
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.model.model.model_phil_scope
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
    """,
            process_includes=True,
        )
        parser = ArgumentParser(phil=phil_scope, check_format=False)
        params, _ = parser.parse_args(args=[], quick_parse=True)
        params.model = model

    from dials.algorithms.scaling.scaler_factory import TargetScalerFactory

    reflections = [reflection_table, target_reflection_table]
    experiment.append(target_experiment[0])
    experiments = create_scaling_model(params, experiment, reflections)
    experiments[-1].scaling_model.set_scaling_model_as_scaled()
    scaler = TargetScalerFactory.create(params, experiments, reflections)
    scaler.perform_scaling()
    scaler.expand_scales_to_all_reflections(calc_cov=True)
    return scaler.unscaled_scalers[0].reflection_table


def scale_single_dataset(reflection_table, experiment, params=None, model="physical"):
    """Determine scale factors for a single dataset. Requires a reflection table
    and an ExperimentList with a single experiment. A custom params option can be
    specified, if not the default scaling params option will be used, with default
    configuration options. The model can be individually specified.

    Returns the reflection table, with added columns 'inverse_scale_factor' and
    'inverse_scale_factor_variance'."""

    if not params:
        phil_scope = phil.parse(
            """
      include scope dials.algorithms.scaling.model.model.model_phil_scope
      include scope dials.algorithms.scaling.scaling_options.phil_scope
      include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
    """,
            process_includes=True,
        )
        parser = ArgumentParser(phil=phil_scope, check_format=False)
        params, _ = parser.parse_args(args=[], quick_parse=True)

    params.model = model

    from dials.algorithms.scaling.scaler_factory import SingleScalerFactory

    experiments = create_scaling_model(params, experiment, [reflection_table])
    scaler = SingleScalerFactory.create(params, experiments[0], reflection_table)
    from dials.algorithms.scaling.algorithm import scaling_algorithm

    scaler = scaling_algorithm(scaler)
    return scaler.reflection_table


def create_scaling_model(params, experiments, reflections):
    """Loop through the experiments, creating the scaling models."""
    autos = [None, Auto, "auto", "Auto"]
    use_auto_model = params.model in autos
    # Determine non-auto model to use outside the loop over datasets.
    if not use_auto_model:
        model_class = None
        for entry_point in pkg_resources.iter_entry_points("dxtbx.scaling_model_ext"):
            if entry_point.name == params.model:
                model_class = entry_point.load()
                break
        if not model_class:
            raise ValueError(f"Unable to create scaling model of type {params.model}")

    for expt, refl in zip(experiments, reflections):
        if not expt.scaling_model or params.overwrite_existing_models:
            # need to make a new model
            if use_auto_model:
                if not expt.scan:
                    model = KBScalingModel
                else:  # set model as physical unless scan < 1.0 degree
                    osc_range = expt.scan.get_oscillation_range()
                    abs_osc_range = abs(osc_range[1] - osc_range[0])
                    if abs_osc_range < 1.0:
                        model = KBScalingModel
                    else:
                        model = PhysicalScalingModel
            else:
                model = model_class
            expt.scaling_model = model.from_data(params, expt, refl)
        else:
            # allow for updating of an existing model.
            expt.scaling_model.update(params)
    return experiments


def create_Ih_table(
    experiments, reflections, selections=None, n_blocks=1, anomalous=False
):
    """Create an Ih table from a list of experiments and reflections. Optionally,
    a selection list can also be given, to select data from each reflection table.
    Allow an unequal number of experiments and reflections, as only need to
    extract one space group value (can optionally check all same if many)."""
    if selections:
        assert len(selections) == len(
            reflections
        ), """Must have an equal number of
    reflection tables and selections in the input lists."""
    space_group_0 = experiments[0].crystal.get_space_group()
    for experiment in experiments:
        assert (
            experiment.crystal.get_space_group() == space_group_0
        ), """The space
    groups of all experiments must be equal."""
    input_tables = []
    indices_lists = []
    for i, reflection in enumerate(reflections):
        if "inverse_scale_factor" not in reflection:
            reflection["inverse_scale_factor"] = flex.double(reflection.size(), 1.0)
        if selections:
            input_tables.append(reflection.select(selections[i]))
            indices_lists.append(selections[i].iselection())
        else:
            input_tables.append(reflection)
            indices_lists = None
    Ih_table = IhTable(
        input_tables,
        space_group_0,
        indices_lists,
        nblocks=n_blocks,
        anomalous=anomalous,
    )
    return Ih_table


def scaled_data_as_miller_array(
    reflection_table_list,
    experiments,
    best_unit_cell=None,
    anomalous_flag=False,
    wavelength=None,
):
    """Get a scaled miller array from an experiment and reflection table."""
    if len(reflection_table_list) > 1:
        joint_table = flex.reflection_table()
        for reflection_table in reflection_table_list:
            # better to just create many miller arrays and join them?
            refl_for_joint_table = flex.reflection_table()
            for col in [
                "miller_index",
                "intensity.scale.value",
                "inverse_scale_factor",
                "intensity.scale.variance",
            ]:
                refl_for_joint_table[col] = reflection_table[col]
            good_refl_sel = ~reflection_table.get_flags(
                reflection_table.flags.bad_for_scaling, all=False
            )
            refl_for_joint_table = refl_for_joint_table.select(good_refl_sel)
            joint_table.extend(refl_for_joint_table)
    else:
        reflection_table = reflection_table_list[0]
        good_refl_sel = ~reflection_table.get_flags(
            reflection_table.flags.bad_for_scaling, all=False
        )
        joint_table = reflection_table.select(good_refl_sel)
    # Filter out negative scale factors to avoid merging statistics errors.
    # These are not removed from the output data, as it is likely one would
    # want to do further analysis e.g. delta cc1/2 and rescaling, to exclude
    # certain data and get better scale factors for all reflections.
    pos_scales = joint_table["inverse_scale_factor"] > 0
    if pos_scales.count(False) > 0:
        logger.info(
            """There are %s reflections with non-positive scale factors which
will not be used for calculating merging statistics""",
            pos_scales.count(False),
        )
        joint_table = joint_table.select(pos_scales)

    if best_unit_cell is None:
        best_unit_cell = determine_best_unit_cell(experiments)
    miller_set = miller.set(
        crystal_symmetry=crystal.symmetry(
            unit_cell=best_unit_cell,
            space_group=experiments[0].crystal.get_space_group(),
            assert_is_compatible_unit_cell=False,
        ),
        indices=joint_table["miller_index"],
        anomalous_flag=anomalous_flag,
    )
    i_obs = miller.array(
        miller_set,
        data=joint_table["intensity.scale.value"] / joint_table["inverse_scale_factor"],
    )
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(
        flex.sqrt(joint_table["intensity.scale.variance"])
        / joint_table["inverse_scale_factor"]
    )
    if not wavelength:
        wavelength = np.mean([expt.beam.get_wavelength() for expt in experiments])
    i_obs.set_info(
        miller.array_info(
            source="DIALS",
            source_type="reflection_tables",
            wavelength=wavelength,
        )
    )
    return i_obs


def determine_best_unit_cell(experiments):
    """Set the median unit cell as the best cell, for consistent d-values across
    experiments."""
    uc_params = [flex.double() for i in range(6)]
    for exp in experiments:
        unit_cell = (
            exp.crystal.get_recalculated_unit_cell() or exp.crystal.get_unit_cell()
        )
        for i, p in enumerate(unit_cell.parameters()):
            uc_params[i].append(p)
    best_unit_cell = uctbx.unit_cell(parameters=[flex.median(p) for p in uc_params])
    if len(experiments) > 1:
        logger.info("Using median unit cell across experiments : %s", best_unit_cell)
    return best_unit_cell


def merging_stats_from_scaled_array(
    scaled_miller_array, n_bins=20, use_internal_variance=False, anomalous=True
):
    """Calculate the normal and anomalous merging statistics."""

    if scaled_miller_array.is_unique_set_under_symmetry():
        raise DialsMergingStatisticsError(
            "Dataset contains no equivalent reflections, merging statistics "
            "cannot be calculated."
        )
    try:
        result = iotbx.merging_statistics.dataset_statistics(
            i_obs=scaled_miller_array,
            n_bins=n_bins,
            anomalous=False,
            sigma_filtering=None,
            eliminate_sys_absent=False,
            use_internal_variance=use_internal_variance,
            cc_one_half_significance_level=0.01,
        )
    except (RuntimeError, Sorry) as e:
        raise DialsMergingStatisticsError(
            f"Error encountered during merging statistics calculation:\n{e}"
        )

    anom_result = None

    if anomalous:
        intensities_anom = scaled_miller_array.as_anomalous_array()
        intensities_anom = intensities_anom.map_to_asu().customized_copy(
            info=scaled_miller_array.info()
        )
        if intensities_anom.is_unique_set_under_symmetry():
            logger.warning(
                "Anomalous dataset contains no equivalent reflections, anomalous "
                "merging statistics cannot be calculated."
            )
        else:
            try:
                anom_result = iotbx.merging_statistics.dataset_statistics(
                    i_obs=intensities_anom,
                    n_bins=n_bins,
                    anomalous=True,
                    sigma_filtering=None,
                    cc_one_half_significance_level=0.01,
                    eliminate_sys_absent=False,
                    use_internal_variance=use_internal_variance,
                )
            except (RuntimeError, Sorry) as e:
                logger.warning(
                    "Error encountered during anomalous merging statistics "
                    "calculation:\n%s",
                    e,
                    exc_info=True,
                )

    return result, anom_result


def intensity_array_from_cif_file(cif_file, anomalous_flag=True, d_min=0.4):
    """Return an intensity miller array from a cif file."""
    xray_structure = cif.reader(file_path=cif_file).file_object.xray_structure_simple()
    ic = (
        xray_structure.structure_factors(
            anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct"
        )
        .f_calc()
        .as_intensity_array()
    )
    return ic


def intensity_array_from_pdb_file(
    pdb_file, anomalous_flag=True, d_min=2.0, wavelength=None
):
    xray_structure = pdb.hierarchy.input(pdb_file).xray_structure_simple()
    if wavelength:
        xray_structure.set_inelastic_form_factors(photon=wavelength, table="sasaki")
    ic = (
        xray_structure.structure_factors(
            anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct"
        )
        .f_calc()
        .as_intensity_array()
    )
    return ic


def create_datastructures_for_target_mtz(experiments, mtz_file, anomalous=False):
    """Read a merged mtz file and extract miller indices, intensities and
    variances."""
    m = mtz.object(mtz_file)
    ind = m.extract_miller_indices()
    cols = m.columns()
    col_dict = {c.label(): c for c in cols}
    r_t = flex.reflection_table()

    def _extract_anom(ind, col_dict):
        table = flex.reflection_table()
        r_tplus = flex.reflection_table()
        r_tminus = flex.reflection_table()
        r_tplus["miller_index"] = ind
        r_tplus["intensity"] = col_dict["I(+)"].extract_values().as_double()
        r_tplus["variance"] = flex.pow2(
            col_dict["SIGI(+)"].extract_values().as_double()
        )
        r_tminus["miller_index"] = ind
        r_tminus["intensity"] = col_dict["I(-)"].extract_values().as_double()
        r_tminus["variance"] = flex.pow2(
            col_dict["SIGI(-)"].extract_values().as_double()
        )
        r_tplus.extend(r_tminus)
        r_tplus.set_flags(
            flex.bool(r_tplus.size(), False), r_tplus.flags.bad_for_scaling
        )
        r_tplus = r_tplus.select(r_tplus["variance"] != 0.0)
        Ih_table = create_Ih_table(
            [experiments[0]], [r_tplus], anomalous=True
        ).blocked_data_list[0]
        table["intensity"] = flumpy.from_numpy(Ih_table.Ih_values)
        inv_var = Ih_table.sum_in_groups(Ih_table.weights, output="per_refl")
        table["variance"] = flumpy.from_numpy(1.0 / inv_var)
        table["miller_index"] = Ih_table.asu_miller_index
        return table

    if anomalous:
        try:
            r_t = _extract_anom(ind, col_dict)
        except KeyError:
            pass
    if not r_t:
        if "I" in col_dict:  # nice and simple
            r_t["miller_index"] = ind
            r_t["intensity"] = col_dict["I"].extract_values().as_double()
            r_t["variance"] = flex.pow2(col_dict["SIGI"].extract_values().as_double())
        elif "IMEAN" in col_dict:  # nice and simple
            r_t["miller_index"] = ind
            r_t["intensity"] = col_dict["IMEAN"].extract_values().as_double()
            r_t["variance"] = flex.pow2(
                col_dict["SIGIMEAN"].extract_values().as_double()
            )
        elif "I(+)" in col_dict:  # need to combine I+ and I- together into target Ih
            if col_dict["I(+)"].n_valid_values() == 0:  # use I(-)
                r_t["miller_index"] = ind
                r_t["intensity"] = col_dict["I(-)"].extract_values().as_double()
                r_t["variance"] = flex.pow2(
                    col_dict["SIGI(-)"].extract_values().as_double()
                )
            elif col_dict["I(-)"].n_valid_values() == 0:  # use I(+)
                r_t["miller_index"] = ind
                r_t["intensity"] = col_dict["I(+)"].extract_values().as_double()
                r_t["variance"] = flex.pow2(
                    col_dict["SIGI(+)"].extract_values().as_double()
                )
            else:  # Combine both - add together then use Ih table to calculate I and sigma
                r_t = _extract_anom(ind, col_dict)
        else:
            raise KeyError("Unable to find intensities (tried I, IMEAN, I(+)/I(-))")
    logger.info(f"Extracted {r_t.size()} intensities from target mtz")
    r_t = r_t.select(r_t["variance"] > 0.0)
    if r_t.size() == 0:
        raise ValueError("No reflections with positive sigma remain after filtering")
    r_t["d"] = (
        miller.set(
            crystal_symmetry=crystal.symmetry(
                space_group=m.space_group(), unit_cell=m.crystals()[0].unit_cell()
            ),
            indices=r_t["miller_index"],
        )
        .d_spacings()
        .data()
    )
    r_t.set_flags(flex.bool(r_t.size(), True), r_t.flags.integrated)

    exp = Experiment()
    exp.crystal = deepcopy(experiments[0].crystal)
    exp.identifier = ersatz_uuid4()
    r_t.experiment_identifiers()[len(experiments)] = exp.identifier
    r_t["id"] = flex.int(r_t.size(), len(experiments))

    # create a new KB scaling model for the target and set as scaled to fix scale
    # for targeted scaling.
    params = Mock()
    params.KB.decay_correction.return_value = False
    exp.scaling_model = KBScalingModel.from_data(params, [], [])
    exp.scaling_model.set_scaling_model_as_scaled()  # Set as scaled to fix scale.
    return exp, r_t


def create_datastructures_for_structural_model(reflections, experiments, model_file):
    """Read a cif/pdb file, calculate intensities and scale them to the average
    intensity of the reflections. Return an experiment and reflection table to
    be used for the structural model in scaling."""

    # read model, compute Fc, square to F^2
    if model_file.endswith(".cif"):
        ic = intensity_array_from_cif_file(model_file)
    elif model_file.endswith(".pdb"):
        wavelength = np.mean([expt.beam.get_wavelength() for expt in experiments])
        ic = intensity_array_from_pdb_file(model_file, wavelength=wavelength)
    else:
        raise ValueError("target_model not a recognised format (.pdb/.cif)")

    exp = deepcopy(experiments[0])
    params = Mock()
    params.decay_correction.return_value = False
    exp.scaling_model = KBScalingModel.from_data(params, [], [])
    exp.scaling_model.set_scaling_model_as_scaled()  # Set as scaled to fix scale.

    rt = flex.reflection_table()
    rt["intensity"] = ic.data()
    rt["miller_index"] = ic.indices()

    exp.identifier = ersatz_uuid4()
    rt.experiment_identifiers()[len(experiments)] = exp.identifier
    rt["id"] = flex.int(rt.size(), len(experiments))

    return exp, rt
