"""
Observers for the scaling algorithm.
"""
from __future__ import absolute_import, division, print_function

import json
import logging
from collections import OrderedDict
from dials.util import tabulate

import six
from cctbx import uctbx
from dials.algorithms.scaling.plots import (
    plot_outliers,
    normal_probability_plot,
    error_model_variance_plot,
    error_regression_plot,
)
from dials.algorithms.scaling.model.model import (
    plot_scaling_models,
    make_combined_plots,
)
from dials.report.analysis import (
    reflection_tables_to_batch_dependent_properties,
    make_merging_statistics_summary,
    table_1_summary,
)
from dials.algorithms.scaling.error_model.error_model import (
    calc_sigmaprime,
    calc_deltahl,
)
from dials.algorithms.scaling.error_model.error_model_target import (
    calculate_regression_x_y,
)
from dials.report.plots import (
    scale_rmerge_vs_batch_plot,
    i_over_sig_i_vs_batch_plot,
    i_over_sig_i_vs_i_plot,
    ResolutionPlotsAndStats,
    IntensityStatisticsPlots,
    AnomalousPlotter,
    make_image_range_table,
)
from dials.algorithms.scaling.scale_and_filter import make_scaling_filtering_plots
from dials.algorithms.scaling.scaling_library import (
    merging_stats_from_scaled_array,
    DialsMergingStatisticsError,
)
from dials.util.batch_handling import batch_manager, get_image_ranges
from dials.util.exclude_images import get_valid_image_ranges
from dials.util.resolution_analysis import resolution_cc_half
from jinja2 import Environment, ChoiceLoader, PackageLoader
from scitbx.array_family import flex

logger = logging.getLogger("dials")


def assert_is_json_serialisable(thing, name, path=None):
    path = path or []
    if isinstance(thing, list):
        for n, element in enumerate(thing):
            assert_is_json_serialisable(element, name, path + [n])
    elif isinstance(thing, dict):
        for key, value in thing.items():
            assert_is_json_serialisable(value, name, path + [repr(key)])
    else:
        try:
            json.dumps(thing)
        except TypeError as e:
            raise TypeError(
                "JSON serialisation error '%s' for value '%s' type %s in %s%s"
                % (
                    e,
                    str(thing),
                    type(thing),
                    name,
                    "".join("[%s]" % step for step in path),
                )
            )


class ScalingSummaryContextManager(object):
    def __init__(self, script):
        self.script = script
        self.data = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.print_scaling_summary(self.script)

    def print_scaling_summary(self, script):
        """Log summary information after scaling."""
        logger.info(print_scaling_model_error_summary(script.experiments))
        valid_ranges = get_valid_image_ranges(script.experiments)
        image_ranges = get_image_ranges(script.experiments)
        msg = []
        for (img, valid, refl) in zip(image_ranges, valid_ranges, script.reflections):
            if valid:
                if len(valid) > 1 or valid[0][0] != img[0] or valid[-1][1] != img[1]:
                    msg.append(
                        "Excluded images for experiment id: %s, image range: %s, limited range: %s"
                        % (
                            refl.experiment_identifiers().keys()[0],
                            list(img),
                            list(valid),
                        )
                    )
        if msg:
            msg = ["Summary of image ranges removed:"] + msg
            logger.info("\n".join(msg))

        # report on partiality of dataset
        partials = flex.double()
        for r in script.reflections:
            if "partiality" in r:
                partials.extend(r["partiality"])
        not_full_sel = partials < 0.99
        not_zero_sel = partials > 0.01
        gt_half = partials > 0.5
        lt_half = partials < 0.5
        partial_gt_half_sel = not_full_sel & gt_half
        partial_lt_half_sel = not_zero_sel & lt_half
        logger.info("Summary of dataset partialities")
        header = ["Partiality (p)", "n_refl"]
        rows = [
            ["all reflections", str(partials.size())],
            ["p > 0.99", str(not_full_sel.count(False))],
            ["0.5 < p < 0.99", str(partial_gt_half_sel.count(True))],
            ["0.01 < p < 0.5", str(partial_lt_half_sel.count(True))],
            ["p < 0.01", str(not_zero_sel.count(False))],
        ]
        logger.info(tabulate(rows, header))
        logger.info(
            """
Reflections below a partiality_cutoff of %s are not considered for any
part of the scaling analysis or for the reporting of merging statistics.
Additionally, if applicable, only reflections with a min_partiality > %s
were considered for use when refining the scaling model.
""",
            script.params.cut_data.partiality_cutoff,
            script.params.reflection_selection.min_partiality,
        )
        stats = script.merging_statistics_result
        if stats:
            anom_stats, cut_stats, cut_anom_stats = (None, None, None)
            if not script.scaled_miller_array.space_group().is_centric():
                anom_stats = script.anom_merging_statistics_result
            logger.info(make_merging_statistics_summary(stats))
            try:
                d_min = resolution_cc_half(stats, limit=0.3).d_min
            except RuntimeError as e:
                logger.debug(f"Resolution fit failed: {e}")
            else:
                max_current_res = stats.bins[-1].d_min
                if d_min and d_min - max_current_res > 0.005:
                    logger.info(
                        "Resolution limit suggested from CC"
                        + "\u00BD"
                        + " fit (limit CC"
                        + "\u00BD"
                        + "=0.3): %.2f",
                        d_min,
                    )
                    try:
                        cut_stats, cut_anom_stats = merging_stats_from_scaled_array(
                            script.scaled_miller_array.resolution_filter(d_min=d_min),
                            script.params.output.merging.nbins,
                            script.params.output.use_internal_variance,
                        )
                    except DialsMergingStatisticsError:
                        pass
                    else:
                        if script.scaled_miller_array.space_group().is_centric():
                            cut_anom_stats = None
            logger.info(table_1_summary(stats, anom_stats, cut_stats, cut_anom_stats))


class ScalingHTMLContextManager(object):
    def __init__(self, script):
        self.script = script
        self.data = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.make_scaling_html(self.script)

    def make_scaling_html(self, scaling_script):
        """Collect data from the individual observers and write the html."""
        html_file = scaling_script.params.output.html
        json_file = scaling_script.params.output.json
        if not (html_file or json_file):
            return
        self.data.update(make_scaling_model_plots(scaling_script.experiments))
        self.data.update(
            make_outlier_plots(scaling_script.reflections, scaling_script.experiments)
        )
        self.data.update(
            make_error_model_plots(scaling_script.params, scaling_script.experiments)
        )
        self.data.update(make_merging_stats_plots(scaling_script))
        self.data.update(make_filtering_plots(scaling_script))
        if html_file:
            logger.info("Writing html report to %s", html_file)
            loader = ChoiceLoader(
                [
                    PackageLoader("dials", "templates"),
                    PackageLoader("dials", "static", encoding="utf-8"),
                ]
            )
            env = Environment(loader=loader)
            template = env.get_template("scaling_report.html")
            assert_is_json_serialisable(self.data, "self.data")
            html = template.render(
                page_title="DIALS scaling report",
                scaling_model_graphs=self.data["scaling_model"],
                scaling_tables=self.data["scaling_tables"],
                error_model_summary=self.data["error_model_summary"],
                resolution_plots=self.data["resolution_plots"],
                scaling_outlier_graphs=self.data["outlier_plots"],
                error_model_plots=self.data["error_model_plots"],
                anom_plots=self.data["anom_plots"],
                batch_plots=self.data["batch_plots"],
                image_range_tables=self.data["image_range_tables"],
                misc_plots=self.data["misc_plots"],
                filter_plots=self.data["filter_plots"],
            )
            with open(html_file, "wb") as f:
                f.write(html.encode("utf-8", "xmlcharrefreplace"))
        if json_file:
            logger.info("Writing html report data to %s", json_file)
            with open(json_file, "w") as outfile:
                json.dump(self.data, outfile)


def make_scaling_model_plots(experiments):
    data = {i: e.scaling_model.to_dict() for i, e in enumerate(experiments)}
    """Generate scaling model component plot data."""
    d = OrderedDict()
    combined_plots = make_combined_plots(data)
    if combined_plots:
        d.update(combined_plots)
    for key in sorted(data.keys()):
        scaling_model_plots = plot_scaling_models(data[key])
        for plot in scaling_model_plots.values():
            plot["layout"]["title"] += " (dataset %s)" % key
        for name, plot in six.iteritems(scaling_model_plots):
            d[name + "_" + str(key)] = plot
    graphs = {"scaling_model": d}
    return graphs


def print_scaling_model_error_summary(experiments):
    """Get a summary of the error distribution of the models."""
    models = [e.scaling_model.to_dict() for e in experiments]
    first_model = models[0]
    component = first_model["configuration_parameters"]["corrections"][0]
    msg = ""
    if "est_standard_devs" in first_model[component]:
        p_sigmas = flex.double()
        for model in models:
            for component in model["configuration_parameters"]["corrections"]:
                if "est_standard_devs" in model[component]:
                    params = flex.double(model[component]["parameters"])
                    sigmas = flex.double(model[component]["est_standard_devs"])
                    null_value = flex.double(
                        len(params), model[component]["null_parameter_value"]
                    )
                    p_sigmas.extend(flex.abs(params - null_value) / sigmas)
        log_p_sigmas = flex.log(p_sigmas)
        frac_high_uncertainty = (log_p_sigmas < 0.69315).count(True) / len(log_p_sigmas)
        if frac_high_uncertainty > 0.5:
            msg = (
                "Warning: Over half ({0:.2f}%) of model parameters have signficant\n"
                "uncertainty (sigma/abs(parameter) > 0.5), which could indicate a\n"
                "poorly-determined scaling problem or overparameterisation.\n"
            ).format(frac_high_uncertainty * 100)
        else:
            msg = (
                "{0:.2f}% of model parameters have signficant uncertainty\n"
                "(sigma/abs(parameter) > 0.5)\n"
            ).format(frac_high_uncertainty * 100)
    return msg


def make_outlier_plots(reflection_tables, experiments):
    """Make outlier plots for the HTML report."""
    data = {}
    for j, (table, expt) in enumerate(zip(reflection_tables, experiments)):
        outliers = table.get_flags(table.flags.outlier_in_scaling)
        x, y, z = table["xyzobs.px.value"].select(outliers).parts()
        if expt.scan:
            zrange = [
                i / expt.scan.get_oscillation()[1]
                for i in expt.scan.get_oscillation_range()
            ]
        else:
            zrange = [0, 0]
        data[j] = {
            "x": list(x),
            "y": list(y),
            "z": list(z),
            "image_size": expt.detector[0].get_image_size(),
            "z_range": zrange,
        }
    d = OrderedDict()
    for key in sorted(data):
        outlier_plots = plot_outliers(data[key])
        for plot in outlier_plots.values():
            if plot:  # may be null if no outliers
                plot["layout"]["title"] += " (dataset %s)" % key
        d["outlier_plot_" + str(key)] = outlier_plots["outlier_xy_positions"]
        d["outlier_plot_z" + str(key)] = outlier_plots["outliers_vs_z"]
    graphs = {"outlier_plots": d}
    return graphs


def make_error_model_plots(params, experiments):
    """Generate normal probability plot data."""
    data = {}
    if experiments[0].scaling_model.error_model:
        em = experiments[0].scaling_model.error_model
        if em.filtered_Ih_table:
            table = em.filtered_Ih_table
            data["intensity"] = table.intensities
            sigmaprime = calc_sigmaprime(em.parameters, table)
            data["delta_hl"] = calc_deltahl(table, table.calc_nh(), sigmaprime)
            data["inv_scale"] = table.inverse_scale_factors
            data["sigma"] = sigmaprime * data["inv_scale"]
            data["binning_info"] = em.binner.binning_info
            em.clear_Ih_table()
            if params.weighting.error_model.basic.minimisation == "regression":
                x, y = calculate_regression_x_y(em.filtered_Ih_table)
                data["regression_x"] = x
                data["regression_y"] = y
                data["model_a"] = em.parameters[0]
                data["model_b"] = em.parameters[1]
        data["summary"] = str(em)

    d = {"error_model_plots": {}, "error_model_summary": "No error model applied"}
    if "delta_hl" in data:
        d["error_model_plots"].update(normal_probability_plot(data))
        d["error_model_plots"].update(
            i_over_sig_i_vs_i_plot(data["intensity"], data["sigma"])
        )
        d["error_model_plots"].update(error_model_variance_plot(data))
        if "regression_x" in data:
            d["error_model_plots"].update(error_regression_plot(data))
    if "summary" in data:
        d["error_model_summary"] = data["summary"]
    return d


def make_filtering_plots(script):
    """Make filtering plots for HTML report"""
    if script.filtering_results:
        data = {
            "merging_stats": script.filtering_results.get_merging_stats(),
            "initial_expids_and_image_ranges": script.filtering_results.initial_expids_and_image_ranges,
            "cycle_results": script.filtering_results.get_cycle_results(),
            "expids_and_image_ranges": script.filtering_results.expids_and_image_ranges,
            "mode": script.params.filtering.deltacchalf.mode,
        }
        d = make_scaling_filtering_plots(data)
        return {"filter_plots": d}
    return {"filter_plots": {}}


def make_merging_stats_plots(script):
    """Make merging stats plots for HTML report"""
    d = {
        "scaling_tables": ([], []),
        "resolution_plots": OrderedDict(),
        "batch_plots": OrderedDict(),
        "misc_plots": OrderedDict(),
        "anom_plots": OrderedDict(),
        "image_range_tables": [],
    }
    if script.merging_statistics_result:
        stats = script.merging_statistics_result
        anom_stats = script.anom_merging_statistics_result
        is_centric = script.scaled_miller_array.space_group().is_centric()
        # Now calculate batch data
        (
            batches,
            rvb,
            isigivb,
            svb,
            batch_data,
        ) = reflection_tables_to_batch_dependent_properties(  # pylint: disable=unbalanced-tuple-unpacking
            script.reflections, script.experiments, script.scaled_miller_array,
        )
        bm = batch_manager(batches, batch_data)
        image_range_tables = make_image_range_table(script.experiments, bm)

        plotter = ResolutionPlotsAndStats(stats, anom_stats, is_centric)
        d["resolution_plots"].update(plotter.make_all_plots())
        d["scaling_tables"] = plotter.statistics_tables()
        d["batch_plots"].update(scale_rmerge_vs_batch_plot(bm, rvb, svb))
        d["batch_plots"].update(i_over_sig_i_vs_batch_plot(bm, isigivb))
        plotter = IntensityStatisticsPlots(
            script.scaled_miller_array, run_xtriage_analysis=False
        )
        d["resolution_plots"].update(plotter.generate_resolution_dependent_plots())
        if d["resolution_plots"]["cc_one_half"]["data"][2]:
            cc_anom = d["resolution_plots"]["cc_one_half"]["data"][2]["y"]
            significance = d["resolution_plots"]["cc_one_half"]["data"][3]["y"]
            sig = flex.double(cc_anom) > flex.double(significance)
            max_anom = 0
            for i, v in enumerate(sig):
                if v:
                    max_anom = i
                else:
                    break
            d_min = uctbx.d_star_sq_as_d(plotter.binner.limits())[max_anom + 1]
        else:
            d_min = 0.0
        d["misc_plots"].update(plotter.generate_miscellanous_plots())
        intensities_anom = script.scaled_miller_array.as_anomalous_array()
        intensities_anom = intensities_anom.map_to_asu().customized_copy(
            info=script.scaled_miller_array.info()
        )
        anom_plotter = AnomalousPlotter(intensities_anom, strong_cutoff=d_min)
        d["anom_plots"].update(anom_plotter.make_plots())
        d["image_range_tables"] = image_range_tables
    return d
