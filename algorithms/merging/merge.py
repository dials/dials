"""Merging functions for experiment lists and reflection tables."""
from __future__ import absolute_import, division, print_function
import logging

from dials.array_family import flex
from dials.algorithms.scaling.scaling_library import (
    scaled_data_as_miller_array,
    merging_stats_from_scaled_array,
)
from dials.algorithms.scaling.scaling_utilities import DialsMergingStatisticsError
from dials.algorithms.scaling.Ih_table import (
    _reflection_table_to_iobs,
    map_indices_to_asu,
)
from dials.algorithms.symmetry.absences.run_absences_checks import (
    run_systematic_absences_checks,
)
from dials.util.filter_reflections import filter_reflection_table
from dials.util.export_mtz import MADMergedMTZWriter, MergedMTZWriter
from dials.report.analysis import make_merging_statistics_summary
from dxtbx.model import ExperimentList
from mmtbx.scaling import data_statistics
from six.moves import cStringIO as StringIO

logger = logging.getLogger("dials")


def prepare_merged_reflection_table(
    experiments, reflection_table, d_min=None, d_max=None
):
    """Filter the data and prepare a reflection table with merged data."""
    if (
        "inverse_scale_factor" in reflection_table
        and "intensity.scale.value" in reflection_table
    ):
        logger.info("Performing systematic absence checks on scaled data")
        reflections = filter_reflection_table(
            reflection_table, intensity_choice=["scale"], d_min=d_min
        )
        reflections["intensity"] = reflections["intensity.scale.value"]
        reflections["variance"] = reflections["intensity.scale.variance"]
    elif "intensity.prf.value" in reflection_table:
        logger.info(
            "Performing systematic absence checks on unscaled profile-integrated data"
        )
        reflections = filter_reflection_table(
            reflection_table, intensity_choice=["profile"], d_min=d_min, d_max=d_max
        )
        reflections["intensity"] = reflections["intensity.prf.value"]
        reflections["variance"] = reflections["intensity.prf.variance"]
    else:
        logger.info(
            "Performing systematic absence checks on unscaled summation-integrated data"
        )
        reflections = filter_reflection_table(
            reflection_table, intensity_choice=["sum"], d_min=d_min, d_max=d_max
        )
        reflections["intensity"] = reflections["intensity.sum.value"]
        reflections["variance"] = reflections["intensity.sum.variance"]

    # now merge
    space_group = experiments[0].crystal.get_space_group()
    reflections["asu_miller_index"] = map_indices_to_asu(
        reflections["miller_index"], space_group
    )
    reflections["inverse_scale_factor"] = flex.double(reflections.size(), 1.0)
    merged = (
        _reflection_table_to_iobs(
            reflections, experiments[0].crystal.get_unit_cell(), space_group
        )
        .merge_equivalents(use_internal_variance=False)
        .array()
    )
    merged_reflections = flex.reflection_table()
    merged_reflections["intensity"] = merged.data()
    merged_reflections["variance"] = flex.pow2(merged.sigmas())
    merged_reflections["miller_index"] = merged.indices()
    return merged_reflections


def make_MAD_merged_mtz_file(
    experiments,
    reflections,
    wavelengths,
    d_min=None,
    d_max=None,
    combine_partials=True,
    partiality_threshold=0.4,
    anomalous=True,
    use_internal_variance=False,
    assess_space_group=False,
    truncate_intensities=True,
    n_residues=200,
    n_bins=20,
    crystal_names=None,
    dataset_names=None,
    project_name=None,
):
    """Make a multi wavelength merged mtz file from experiments and reflections."""
    # need to add a crystal to the mtz object
    # now go through data selecting on wavelength - loop over each to get mtz_object
    # Create the mtz file

    mtz_writer = MADMergedMTZWriter(
        experiments[0].crystal.get_space_group(), experiments[0].crystal.get_unit_cell()
    )

    # now add each wavelength.
    if not dataset_names or len(dataset_names) != len(wavelengths):
        logger.info(
            "Unequal number of dataset names and wavelengths, using default naming."
        )
        dataset_names = [None] * len(wavelengths)
    if not crystal_names or len(crystal_names) != len(wavelengths):
        logger.info(
            "Unequal number of crystal names and wavelengths, using default naming."
        )
        crystal_names = [None] * len(wavelengths)

    for dname, cname, (wavelength, exp_nos) in zip(
        dataset_names, crystal_names, wavelengths.items()
    ):
        expids = [experiments[i].identifier for i in exp_nos]
        new_exps = ExperimentList([experiments[i] for i in exp_nos])
        sel_reflections = reflections[0].select_on_experiment_identifiers(expids)

        logger.info("Running merge for wavelength: %s", wavelength)

        # merge and truncate the data
        amplitudes, anomalous_amplitudes = None, None
        merged_array, merged_anomalous_array, stats_summary = merge(
            new_exps,
            [sel_reflections],
            d_min=d_min,
            d_max=d_max,
            combine_partials=combine_partials,
            partiality_threshold=partiality_threshold,
            anomalous=anomalous,
            assess_space_group=assess_space_group,
            n_bins=n_bins,
            use_internal_variance=use_internal_variance,
        )
        if anomalous:
            merged_intensities = merged_anomalous_array
        else:
            merged_intensities = merged_array

        if truncate_intensities:
            amplitudes, anomalous_amplitudes = truncate(merged_intensities)
        show_wilson_scaling_analysis(merged_intensities, n_residues=n_residues)
        if stats_summary:
            logger.info(stats_summary)

        #### Add each wavelength as a new crystal.
        mtz_writer.add_crystal(crystal_name=cname, project_name=project_name)
        mtz_writer.add_empty_dataset(wavelength, name=dname)
        mtz_writer.add_dataset(
            merged_array, merged_anomalous_array, amplitudes, anomalous_amplitudes
        )

    return mtz_writer.mtz_file


def make_merged_mtz_file(
    wavelength,
    merged_array,
    merged_anomalous_array=None,
    amplitudes=None,
    anomalous_amplitudes=None,
    crystal_name=None,
    dataset_name=None,
    project_name=None,
):
    """Make an mtz object for the data, adding the date, time and program."""

    assert merged_array.is_xray_intensity_array()

    mtz_writer = MergedMTZWriter(merged_array.space_group(), merged_array.unit_cell())
    mtz_writer.add_crystal(
        crystal_name=crystal_name, project_name=project_name,
    )
    mtz_writer.add_empty_dataset(wavelength, name=dataset_name)
    mtz_writer.add_dataset(
        merged_array, merged_anomalous_array, amplitudes, anomalous_amplitudes
    )

    return mtz_writer.mtz_file


def merge(
    experiments,
    reflections,
    d_min=None,
    d_max=None,
    combine_partials=True,
    partiality_threshold=0.4,
    anomalous=True,
    use_internal_variance=False,
    assess_space_group=False,
    n_bins=20,
):
    """Filter data, assess space group, run french wilson and Wilson stats."""

    logger.info("\nMerging scaled reflection data\n")
    # first filter bad reflections using dials.util.filter methods
    reflections = filter_reflection_table(
        reflections[0],
        intensity_choice=["scale"],
        d_min=d_min,
        d_max=d_max,
        combine_partials=combine_partials,
        partiality_threshold=partiality_threshold,
    )
    # ^ scale factor has been applied, so now set to 1.0 - okay as not
    # going to output scale factor in merged mtz.
    reflections["inverse_scale_factor"] = flex.double(reflections.size(), 1.0)

    scaled_array = scaled_data_as_miller_array([reflections], experiments)
    # Note, merge_equivalents does not raise an error if data is unique.
    merged = scaled_array.merge_equivalents(
        use_internal_variance=use_internal_variance
    ).array()
    merged_anom = None

    if anomalous:
        anomalous_scaled = scaled_array.as_anomalous_array()
        merged_anom = anomalous_scaled.merge_equivalents(
            use_internal_variance=use_internal_variance
        ).array()

    # Before merge, do some assessment of the space_group
    if assess_space_group:
        merged_reflections = flex.reflection_table()
        merged_reflections["intensity"] = merged.data()
        merged_reflections["variance"] = flex.pow2(merged.sigmas())
        merged_reflections["miller_index"] = merged.indices()
        logger.info("Running systematic absences check")
        run_systematic_absences_checks(experiments, merged_reflections)

    # here return merged arrays?
    # Show merging stats again.
    try:
        stats, anom_stats = merging_stats_from_scaled_array(
            scaled_array, n_bins, use_internal_variance,
        )
    except DialsMergingStatisticsError as e:
        logger.error(e, exc_info=True)
        stats_summary = None
    else:
        if anomalous and anom_stats:
            stats_summary = make_merging_statistics_summary(anom_stats)
        else:
            stats_summary = make_merging_statistics_summary(stats)

    return merged, merged_anom, stats_summary


def show_wilson_scaling_analysis(merged_intensities, n_residues=200):
    if not merged_intensities.space_group().is_centric():
        try:
            wilson_scaling = data_statistics.wilson_scaling(
                miller_array=merged_intensities, n_residues=n_residues
            )
        except (IndexError, RuntimeError) as e:
            logger.error(
                "\n"
                "Error encountered during Wilson statistics calculation:\n"
                "Perhaps there are too few unique reflections.\n"
                "%s",
                e,
                exc_info=True,
            )
        else:
            # Divert output through logger - do with StringIO rather than
            # info_handle else get way too much whitespace in output.
            out = StringIO()
            wilson_scaling.show(out=out)
            logger.info(out.getvalue())


def truncate(merged_intensities):
    logger.info("\nScaling input intensities via French-Wilson Method")
    out = StringIO()
    if merged_intensities.anomalous_flag():
        anom_amplitudes = merged_intensities.french_wilson(log=out)
        n_removed = merged_intensities.size() - anom_amplitudes.size()
        assert anom_amplitudes.is_xray_amplitude_array()
        amplitudes = anom_amplitudes.as_non_anomalous_array()
        amplitudes = amplitudes.merge_equivalents().array()
    else:
        anom_amplitudes = None
        amplitudes = merged_intensities.french_wilson(log=out)
        n_removed = merged_intensities.size() - amplitudes.size()
    logger.info("Total number of rejected intensities %s", n_removed)
    logger.debug(out.getvalue())
    return amplitudes, anom_amplitudes
