"""Merging functions for experiment lists and reflection tables."""

from __future__ import annotations

import logging
from io import StringIO

from jinja2 import ChoiceLoader, Environment, PackageLoader

from cctbx import uctbx
from mmtbx.scaling import data_statistics

from dials.algorithms.scaling.Ih_table import (
    _reflection_table_to_iobs,
    map_indices_to_asu,
)
from dials.algorithms.scaling.scaling_library import (
    merging_stats_from_scaled_array,
    scaled_data_as_miller_array,
)
from dials.algorithms.scaling.scaling_utilities import DialsMergingStatisticsError
from dials.algorithms.symmetry.absences.run_absences_checks import (
    run_systematic_absences_checks,
)
from dials.array_family import flex
from dials.report.analysis import make_merging_statistics_summary, table_1_summary
from dials.report.plots import d_star_sq_to_d_ticks
from dials.util import tabulate
from dials.util.export_mtz import MADMergedMTZWriter, MergedMTZWriter
from dials.util.filter_reflections import filter_reflection_table

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


class MTZDataClass:

    """Container class (i.e. Python 3.7 dataclass) for per-wavelength mtz dataset."""

    def __init__(
        self,
        wavelength=0.0,
        project_name="AUTOMATIC",
        dataset_name="NATIVE",
        crystal_name="XTAL",
        merged_array=None,
        merged_anomalous_array=None,
        amplitudes=None,
        anomalous_amplitudes=None,
        dano=None,
        multiplicities=None,
    ):
        self.wavelength = wavelength
        self.project_name = project_name
        self.dataset_name = dataset_name
        self.crystal_name = crystal_name
        self.merged_array = merged_array
        self.merged_anomalous_array = merged_anomalous_array
        self.amplitudes = amplitudes
        self.anomalous_amplitudes = anomalous_amplitudes
        self.dano = dano
        self.multiplicities = multiplicities


def make_merged_mtz_file(mtz_datasets):
    """
    Make an mtz file object for the data, adding the date, time and program.

    For multi-wavelength data, each wavelength is added as a new crystal.

    Args:
        mtz_datasets: A list of MTZDataClass objects, one per wavelength of the
            experiment.

    Returns:
        An iotbx mtz file object.
    """

    if len(mtz_datasets) > 1:
        writer = MADMergedMTZWriter
    else:
        writer = MergedMTZWriter

    mtz_writer = writer(
        mtz_datasets[0].merged_array.space_group(),
        mtz_datasets[0].merged_array.unit_cell(),
    )

    #### Add each wavelength as a new crystal.
    for dataset in mtz_datasets:
        mtz_writer.add_crystal(
            crystal_name=dataset.crystal_name, project_name=dataset.project_name
        )
        mtz_writer.add_empty_dataset(dataset.wavelength, name=dataset.dataset_name)
        mtz_writer.add_dataset(
            dataset.merged_array,
            dataset.merged_anomalous_array,
            dataset.amplitudes,
            dataset.anomalous_amplitudes,
            dataset.dano,
            dataset.multiplicities,
        )

    return mtz_writer.mtz_file


def merge(
    experiments,
    reflections,
    d_min=None,
    d_max=None,
    combine_partials=True,
    partiality_threshold=0.4,
    best_unit_cell=None,
    anomalous=True,
    use_internal_variance=False,
    assess_space_group=False,
    n_bins=20,
):
    """
    Merge reflection table data and generate a summary of the merging statistics.

    This procedure filters the input data, merges the data (normal and optionally
    anomalous), assesses the space group symmetry and generates a summary
    of the merging statistics.

    Returns two merge_equivalents objects and a statistics summary.
    """

    logger.info("\nMerging scaled reflection data\n")
    # first filter bad reflections using dials.util.filter methods
    reflections = filter_reflection_table(
        reflections,
        intensity_choice=["scale"],
        d_min=d_min,
        d_max=d_max,
        combine_partials=combine_partials,
        partiality_threshold=partiality_threshold,
    )
    # ^ scale factor has been applied, so now set to 1.0 - okay as not
    # going to output scale factor in merged mtz.
    reflections["inverse_scale_factor"] = flex.double(reflections.size(), 1.0)
    scaled_array = scaled_data_as_miller_array(
        [reflections], experiments, best_unit_cell
    )
    # Note, merge_equivalents does not raise an error if data is unique.
    merged = scaled_array.merge_equivalents(use_internal_variance=use_internal_variance)
    merged_anom = None

    if anomalous:
        anomalous_scaled = scaled_array.as_anomalous_array()
        merged_anom = anomalous_scaled.merge_equivalents(
            use_internal_variance=use_internal_variance
        )

    # Before merge, do assessment of the space_group
    if assess_space_group:
        merged_reflections = flex.reflection_table()
        merged_reflections["intensity"] = merged.array().data()
        merged_reflections["variance"] = flex.pow2(merged.array().sigmas())
        merged_reflections["miller_index"] = merged.array().indices()
        logger.info("Running systematic absences check")
        run_systematic_absences_checks(experiments, merged_reflections)

    try:
        stats, anom_stats = merging_stats_from_scaled_array(
            scaled_array,
            n_bins,
            use_internal_variance,
        )
    except DialsMergingStatisticsError as e:
        logger.error(e, exc_info=True)
        stats_summary = None
    else:
        stats_summary = make_merging_statistics_summary(stats)
        stats_summary += table_1_summary(stats, anom_stats)

    return merged, merged_anom, stats_summary


def show_wilson_scaling_analysis(merged_intensities, n_residues=200):
    """
    Report the wilson statistics for a merged intensity array

    Args:
        merged_intensities: A merged miller intensity array.
        n_residues: The number of residues to use for the wilson analysis.
    """
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
    """
    Perform French-Wilson truncation procedure on merged intensities.

    Args:
        merged_intensities: A merged miller intensity array (normal or anomalous).

    Returns:
        (tuple): tuple containing:
            amplitudes: A normal all-positive miller amplitude array
            anom_amplitudes: An anomalous all-positive amplitude array, if the
                input array has the anomalous_flag set, else None.
            dano: The array of anomalous differences, if the input array has the
                anomalous_flag set, else None.
    """
    logger.info("\nPerforming French-Wilson treatment of scaled intensities")
    out = StringIO()
    if merged_intensities.anomalous_flag():
        anom_amplitudes = merged_intensities.french_wilson(log=out)
        n_removed = merged_intensities.size() - anom_amplitudes.size()
        assert anom_amplitudes.is_xray_amplitude_array()
        amplitudes = anom_amplitudes.as_non_anomalous_array()
        amplitudes = amplitudes.merge_equivalents(use_internal_variance=False).array()
        dano = anom_amplitudes.anomalous_differences()
    else:
        anom_amplitudes = None
        dano = None
        amplitudes = merged_intensities.french_wilson(log=out)
        n_removed = merged_intensities.size() - amplitudes.size()
    logger.info("Total number of rejected intensities %s", n_removed)
    logger.debug(out.getvalue())
    return amplitudes, anom_amplitudes, dano


def dano_over_sigdano_stats(anomalous_amplitudes, n_bins=20):
    """Calculate the statistic for resolution bins and overall."""
    vals = flex.double()
    # First calculate dF/s(dF) per resolution bin
    anomalous_amplitudes.setup_binner(n_bins=n_bins)
    resolution_bin_edges = flex.double()
    for i_bin in anomalous_amplitudes.binner().range_used():
        sel = anomalous_amplitudes.binner().selection(i_bin)
        arr = anomalous_amplitudes.select(sel)
        vals.append(dano_over_sigdano(arr))
        resolution_bin_edges.append(anomalous_amplitudes.binner().bin_d_min(i_bin))
    resolution_bin_edges.append(anomalous_amplitudes.binner().d_min())
    return vals, resolution_bin_edges


def dano_over_sigdano(anomalous_amplitudes):
    """Calculate < |F(+) - F(-)| / sigma(F(+) - F(-))> i.e. <DANO/SIGDANO>."""
    diff = anomalous_amplitudes.anomalous_differences()
    if not diff.data() or not diff.sigmas():
        return 0.0
    return flex.mean(flex.abs(diff.data()) / diff.sigmas())


def make_dano_table(anomalous_amplitudes):
    """Calculate <dano/sigdano> in resolution bins and tabulate."""
    dFsdF, resolution_bin_edges = dano_over_sigdano_stats(anomalous_amplitudes)

    logger.info("Size of anomalous differences")
    header = ["d_max", "d_min", "<|ΔF|/σ(ΔF)>"]
    rows = []
    for i, dF in enumerate(dFsdF):
        rows.append(
            [
                f"{resolution_bin_edges[i]:6.2f}",
                f"{resolution_bin_edges[i+1]:6.2f}",
                f"{dF:6.3f}",
            ]
        )
    return tabulate(rows, header)


def make_dano_plots(anomalous_data):
    """
    Make dicts of data for plotting e.g. for plotly.

    Args:
        anomalous_data (dict) : A dict of (wavelength, anomalous array) data.

    Returns:
        dict: A dictionary containing the plotting data.
    """

    data = {
        "dF": {
            "dano": {
                "data": [],
                "help": """\
This plot shows the size of the anomalous differences of F relative to the uncertainties,
(<|F(+)-F(-)|/σ(F(+)-F(-))>). A value of 0.8 is indicative of pure noise, and
a suggested cutoff is when the value falls below 1.2, although these guides require
reliable sigma estimates. For further information see
https://strucbio.biologie.uni-konstanz.de/ccp4wiki/index.php?title=SHELX_C/D/E
""",
            },
        },
    }

    for i, (wave, anom) in enumerate(anomalous_data.items()):
        dFsdF, resolution_bin_edges = dano_over_sigdano_stats(anom)
        d_star_sq_bins = [
            0.5
            * (
                uctbx.d_as_d_star_sq(resolution_bin_edges[i])
                + uctbx.d_as_d_star_sq(resolution_bin_edges[i + 1])
            )
            for i in range(0, len(resolution_bin_edges[:-1]))
        ]
        d_star_sq_tickvals, d_star_sq_ticktext = d_star_sq_to_d_ticks(
            d_star_sq_bins, nticks=5
        )
        data["dF"]["dano"]["data"].append(
            {
                "x": d_star_sq_bins,
                "y": list(dFsdF),
                "type": "scatter",
                "name": "\u03BB" + f"={wave:.4f}",
            }
        )

    data["dF"]["dano"]["data"].append(
        {
            "x": d_star_sq_bins,
            "y": [0.8] * len(d_star_sq_bins),
            "type": "scatter",
            "mode": "lines",
            "name": "random noise level",
        }
    )
    data["dF"]["dano"]["data"].append(
        {
            "x": d_star_sq_bins,
            "y": [1.2] * len(d_star_sq_bins),
            "type": "scatter",
            "mode": "lines",
            "name": "an approximate <br>threshold for a<br>resolution cutoff",
        }
    )

    data["dF"]["dano"]["layout"] = {
        "title": "<|ΔF|/σ(ΔF)> vs resolution",
        "xaxis": {
            "title": "Resolution (Å)",
            "tickvals": d_star_sq_tickvals,
            "ticktext": d_star_sq_ticktext,
        },
        "yaxis": {"title": "<|ΔF|/σ(ΔF)>", "rangemode": "tozero"},
    }
    return data


def generate_html_report(mtz_file, filename):
    """Make a html report to plot dano/sigdano."""
    anom_data = {}
    # Collect F+, F-, SigF+, SigF- data, allowing for multiple wavelengths.
    for array in mtz_file.as_miller_arrays():
        if (
            len(array.info().labels) == 4
            and array.info().type_hints_from_file == "amplitude"
        ):
            anom_data[array.info().wavelength] = array
    data = {"dF": {}}
    if not anom_data:
        return
    data = make_dano_plots(anom_data)

    loader = ChoiceLoader(
        [
            PackageLoader("dials", "templates"),
            PackageLoader("dials", "static", encoding="utf-8"),
        ]
    )
    env = Environment(loader=loader)
    template = env.get_template("simple_report.html")
    html = template.render(
        page_title="DIALS merge report",
        panel_title="Anomalous signal",
        panel_id="Anomalous signal",
        graphs=data["dF"],
    )
    logger.info("Writing html report to %s", filename)
    with open(filename, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))
