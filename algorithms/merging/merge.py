"""Merging functions for experiment lists and reflection tables."""
import logging
from cStringIO import StringIO
from dials.array_family import flex
from dials.command_line.space_group import run_sys_abs_checks
from dials.algorithms.scaling.scaling_library import (
    scaled_data_as_miller_array,
    merging_stats_from_scaled_array,
)
from dials.util.filter_reflections import filter_reflection_table
from dials.util.export_mtz import MADMergedMTZWriter, MergedMTZWriter
from dials.report.analysis import make_merging_statistics_summary
from dxtbx.model import ExperimentList
from mmtbx.scaling import data_statistics

logger = logging.getLogger("dials")


def make_MAD_merged_mtz_file(params, experiments, reflections, wavelengths):
    """Make a multi wavelength merged mtz file from experiments and reflections."""
    # need to add a crystal to the mtz object
    # now go through data selecting on wavelength - loop over each to get mtz_object
    # Create the mtz file

    mtz_writer = MADMergedMTZWriter(
        experiments[0].crystal.get_space_group(), experiments[0].crystal.get_unit_cell()
    )

    # now add each wavelength.
    for wavelength, exp_nos in wavelengths.iteritems():
        expids = []

        new_exps = ExperimentList()
        for i in exp_nos:
            expids.append(experiments[i].identifier)  # string
            new_exps.append(experiments[i])
        refls = reflections[0].select_on_experiment_identifiers(expids)

        logger.info("Running merge for wavelength: %s", wavelength)
        merged, anom, amplitudes, anom_amp = merge_and_truncate(
            params, new_exps, [refls]
        )
        #### Add each wavelength as a new crystal.
        mtz_writer.add_crystal()
        mtz_writer.add_dataset(wavelength, merged, anom, amplitudes, anom_amp)

    return mtz_writer.mtz_file


def make_merged_mtz_file(
    wavelength,
    merged_array,
    merged_anomalous_array=None,
    amplitudes=None,
    anomalous_amplitudes=None,
):
    """Make an mtz object for the data, adding the date, time and program."""

    assert merged_array.is_xray_intensity_array()

    mtz_writer = MergedMTZWriter(merged_array.space_group(), merged_array.unit_cell())
    mtz_writer.add_crystal(crystal_name="DIALS")
    mtz_writer.add_dataset(
        merged_array,
        merged_anomalous_array,
        amplitudes,
        anomalous_amplitudes,
        wavelength,
    )

    return mtz_writer.mtz_file


def merge_and_truncate(params, experiments, reflections):
    """Filter data, assess space group, run french wilson and Wilson stats."""

    logger.info("\nMerging scaled reflection data\n")
    # first filter bad reflections using dials.util.filter methods
    reflections = filter_reflection_table(
        reflections[0],
        intensity_choice=["scale"],
        d_min=params.d_min,
        combine_partials=params.combine_partials,
        partiality_threshold=params.partiality_threshold,
    )
    # ^ scale factor has been applied, so now set to 1.0 - okay as not
    # going to output scale factor in merged mtz.
    reflections["inverse_scale_factor"] = flex.double(reflections.size(), 1.0)

    scaled_array = scaled_data_as_miller_array([reflections], experiments)
    if params.anomalous:
        anomalous_scaled = scaled_array.as_anomalous_array()

    merged = scaled_array.merge_equivalents(
        use_internal_variance=params.merging.use_internal_variance
    ).array()
    merged_anom = None
    if params.anomalous:
        merged_anom = anomalous_scaled.merge_equivalents(
            use_internal_variance=params.merging.use_internal_variance
        ).array()

    # Before merge, do some assessment of the space_group
    if params.assess_space_group:
        merged_reflections = flex.reflection_table()
        merged_reflections["intensity"] = merged.data()
        merged_reflections["variance"] = merged.sigmas() ** 2
        merged_reflections["miller_index"] = merged.indices()
        logger.info("Running systematic absences check")
        run_sys_abs_checks(experiments, merged_reflections)

    # Run the stats on truncating on anomalous or non anomalous?
    if params.anomalous:
        intensities = merged_anom
    else:
        intensities = merged

    assert intensities.is_xray_intensity_array()
    amplitudes = None
    anom_amplitudes = None
    if params.truncate:
        logger.info("\nScaling input intensities via French-Wilson Method")
        out = StringIO()
        if params.anomalous:
            anom_amplitudes = intensities.french_wilson(params=params, log=out)
            n_removed = intensities.size() - anom_amplitudes.size()
            assert anom_amplitudes.is_xray_amplitude_array()
            amplitudes = anom_amplitudes.as_non_anomalous_array()
            amplitudes = amplitudes.merge_equivalents().array()
        else:
            amplitudes = intensities.french_wilson(params=params, log=out)
            n_removed = intensities.size() - amplitudes.size()
        logger.info("Total number of rejected intensities %s", n_removed)
        logger.debug(out.getvalue())

    if params.reporting.wilson_stats:
        if not intensities.space_group().is_centric():
            wilson_scaling = data_statistics.wilson_scaling(
                miller_array=intensities, n_residues=params.n_residues
            )  # XXX default n_residues?
            # Divert output through logger - do with StringIO rather than
            # info_handle else get way too much whitespace in output.
            out = StringIO()
            wilson_scaling.show(out=out)
            logger.info(out.getvalue())

    # Apply wilson B to give absolute scale?

    # Show merging stats again.
    if params.reporting.merging_stats:
        stats, anom_stats = merging_stats_from_scaled_array(
            scaled_array, params.merging.n_bins, params.merging.use_internal_variance
        )
        if params.merging.anomalous:
            logger.info(make_merging_statistics_summary(anom_stats))
        else:
            logger.info(make_merging_statistics_summary(stats))

    return merged, merged_anom, amplitudes, anom_amplitudes
