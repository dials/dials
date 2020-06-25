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
from dials.util.observer import Observer, singleton
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
from dials.util.batch_handling import batch_manager, get_image_ranges
from dials.util.exclude_images import get_valid_image_ranges
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


def register_default_scaling_observers(script):
    """Register the standard observers to the scaling script."""
    script.register_observer(
        event="merging_statistics", observer=MergingStatisticsObserver()
    )
    script.register_observer(
        event="run_script",
        observer=ScalingSummaryGenerator(),
        callback="print_scaling_summary",
    )
    script.register_observer(
        event="run_script",
        observer=ScalingHTMLGenerator(),
        callback="make_scaling_html",
    )
    register_scaler_observers(script.scaler)


def register_scaler_observers(scaler):
    """Register observers on the scaler."""
    scaler.register_observer(
        event="performed_error_analysis", observer=ErrorModelObserver()
    )

    scaler.register_observer(event="performed_scaling", observer=ScalingModelObserver())
    scaler.register_observer(
        event="performed_outlier_rejection", observer=ScalingOutlierObserver()
    )


def register_merging_stats_observers(script):
    """Register only obsevers needed to record and print merging stats."""
    script.register_observer(
        event="merging_statistics", observer=MergingStatisticsObserver()
    )
    script.register_observer(
        event="run_script",
        observer=ScalingSummaryGenerator(),
        callback="print_scaling_summary",
    )


def register_scale_and_filter_observers(script):
    script.register_observer(event="run_scale_and_filter", observer=FilteringObserver())
    script.register_observer(
        event="run_scale_and_filter",
        observer=ScalingHTMLGenerator(),
        callback="make_scaling_html",
    )
    try:
        script.unregister_observer(event="run_script", observer=ScalingHTMLGenerator())
    except KeyError:
        pass


@singleton
class ScalingSummaryGenerator(Observer):
    """
    Observer to summarise data
    """

    def print_scaling_summary(self, scaling_script):
        """Log summary information after scaling."""
        if ScalingModelObserver().data:
            logger.info(ScalingModelObserver().return_model_error_summary())
        valid_ranges = get_valid_image_ranges(scaling_script.experiments)
        image_ranges = get_image_ranges(scaling_script.experiments)
        msg = []
        for (img, valid, refl) in zip(
            image_ranges, valid_ranges, scaling_script.reflections
        ):
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
        for r in scaling_script.reflections:
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
            scaling_script.params.cut_data.partiality_cutoff,
            scaling_script.params.reflection_selection.min_partiality,
        )
        if MergingStatisticsObserver().data:
            logger.info(
                make_merging_statistics_summary(
                    MergingStatisticsObserver().data["statistics"],
                    MergingStatisticsObserver().data["anomalous_statistics"],
                )
            )


@singleton
class ScalingHTMLGenerator(Observer):
    """
    Observer to make a html report
    """

    def make_scaling_html(self, scaling_script):
        """Collect data from the individual observers and write the html."""
        html_file = scaling_script.params.output.html
        json_file = scaling_script.params.output.json
        if not (html_file or json_file):
            return
        self.data.update(ScalingModelObserver().make_plots())
        self.data.update(ScalingOutlierObserver().make_plots())
        self.data.update(ErrorModelObserver().make_plots())
        self.data.update(MergingStatisticsObserver().make_plots())
        self.data.update(FilteringObserver().make_plots())
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


@singleton
class ScalingModelObserver(Observer):
    """
    Observer to record scaling model data and make model plots.
    """

    def update(self, scaler):
        """Update the data in the observer."""
        active_scalers = getattr(scaler, "active_scalers", False)
        if not active_scalers:
            active_scalers = [scaler]
        for i, s in enumerate(active_scalers):
            self.data[i] = s.experiment.scaling_model.to_dict()

    def make_plots(self):
        """Generate scaling model component plot data."""
        d = OrderedDict()
        combined_plots = make_combined_plots(self.data)
        if combined_plots:
            d.update(combined_plots)
        for key in sorted(self.data.keys()):
            scaling_model_plots = plot_scaling_models(self.data[key])
            for plot in scaling_model_plots.values():
                plot["layout"]["title"] += " (dataset %s)" % key
            for name, plot in six.iteritems(scaling_model_plots):
                d[name + "_" + str(key)] = plot
        graphs = {"scaling_model": d}
        return graphs

    def return_model_error_summary(self):
        """Get a summary of the error distribution of the models."""
        first_model = list(self.data.values())[0]
        component = first_model["configuration_parameters"]["corrections"][0]
        msg = ""
        if "est_standard_devs" in first_model[component]:
            p_sigmas = flex.double()
            for model in self.data.values():
                for component in model["configuration_parameters"]["corrections"]:
                    if "est_standard_devs" in model[component]:
                        params = flex.double(model[component]["parameters"])
                        sigmas = flex.double(model[component]["est_standard_devs"])
                        null_value = flex.double(
                            len(params), model[component]["null_parameter_value"]
                        )
                        p_sigmas.extend(flex.abs(params - null_value) / sigmas)
            log_p_sigmas = flex.log(p_sigmas)
            frac_high_uncertainty = (log_p_sigmas < 0.69315).count(True) / len(
                log_p_sigmas
            )
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


@singleton
class ScalingOutlierObserver(Observer):
    """
    Observer to record scaling outliers and make outlier plots.
    """

    def update(self, scaler):
        active_scalers = getattr(scaler, "active_scalers")
        if not active_scalers:
            active_scalers = [scaler]
        for j, scaler in enumerate(active_scalers):
            outlier_isel = scaler.suitable_refl_for_scaling_sel.iselection().select(
                scaler.outliers
            )
            x, y, z = (
                scaler.reflection_table["xyzobs.px.value"].select(outlier_isel).parts()
            )
            if scaler.experiment.scan:
                zrange = [
                    i / scaler.experiment.scan.get_oscillation()[1]
                    for i in scaler.experiment.scan.get_oscillation_range()
                ]
            else:
                zrange = [0, 0]
            self.data[j] = {
                "x": list(x),
                "y": list(y),
                "z": list(z),
                "image_size": scaler.experiment.detector[0].get_image_size(),
                "z_range": zrange,
            }

    def make_plots(self):
        """Generate plot data of outliers on the detector and vs z."""
        d = OrderedDict()
        for key in sorted(self.data):
            outlier_plots = plot_outliers(self.data[key])
            for plot in outlier_plots.values():
                if plot:  # may be null if no outliers
                    plot["layout"]["title"] += " (dataset %s)" % key
            d["outlier_plot_" + str(key)] = outlier_plots["outlier_xy_positions"]
            d["outlier_plot_z" + str(key)] = outlier_plots["outliers_vs_z"]
        graphs = {"outlier_plots": d}
        return graphs


@singleton
class ErrorModelObserver(Observer):
    """
    Observer to record scaling error model data and make a plot.
    """

    def update(self, scaler):
        if scaler.error_model:
            if scaler.error_model.filtered_Ih_table:
                table = scaler.error_model.filtered_Ih_table
                self.data["intensity"] = table.intensities
                sigmaprime = calc_sigmaprime(scaler.error_model.parameters, table)
                self.data["delta_hl"] = calc_deltahl(table, table.calc_nh(), sigmaprime)
                self.data["inv_scale"] = table.inverse_scale_factors
                self.data["sigma"] = sigmaprime * self.data["inv_scale"]
                self.data["binning_info"] = scaler.error_model.binner.binning_info
                scaler.error_model.clear_Ih_table()
            if scaler.params.weighting.error_model.basic.minimisation == "regression":
                x, y = calculate_regression_x_y(scaler.error_model.filtered_Ih_table)
                self.data["regression_x"] = x
                self.data["regression_y"] = y
                self.data["model_a"] = scaler.error_model.parameters[0]
                self.data["model_b"] = scaler.error_model.parameters[1]
            self.data["summary"] = str(scaler.error_model)

    def make_plots(self):
        """Generate normal probability plot data."""
        d = {"error_model_plots": {}, "error_model_summary": "No error model applied"}
        if "delta_hl" in self.data:
            d["error_model_plots"].update(normal_probability_plot(self.data))
            d["error_model_plots"].update(
                i_over_sig_i_vs_i_plot(self.data["intensity"], self.data["sigma"])
            )
            d["error_model_plots"].update(error_model_variance_plot(self.data))
            if "regression_x" in self.data:
                d["error_model_plots"].update(error_regression_plot(self.data))
        if "summary" in self.data:
            d["error_model_summary"] = self.data["summary"]
        return d


@singleton
class FilteringObserver(Observer):
    """
    Observer to record data from the scaling and filtering algorithm.
    """

    def update(self, scaling_script):
        if scaling_script.filtering_results:
            self.data = {
                "merging_stats": scaling_script.filtering_results.get_merging_stats(),
                "initial_expids_and_image_ranges": scaling_script.filtering_results.initial_expids_and_image_ranges,
                "cycle_results": scaling_script.filtering_results.get_cycle_results(),
                "expids_and_image_ranges": scaling_script.filtering_results.expids_and_image_ranges,
                "mode": scaling_script.params.filtering.deltacchalf.mode,
            }

    def make_plots(self):
        """Make plots for scale and filter."""
        if not self.data:
            return {"filter_plots": {}}
        # Make merging stats plots, histograms and image ranges.
        d = make_scaling_filtering_plots(self.data)
        return {"filter_plots": d}


@singleton
class MergingStatisticsObserver(Observer):
    """
    Observer to record merging statistics data and make tables.
    """

    def update(self, scaling_script):
        if scaling_script.merging_statistics_result:
            self.data = {
                "statistics": scaling_script.merging_statistics_result,
                "anomalous_statistics": scaling_script.anom_merging_statistics_result,
                "is_centric": scaling_script.scaled_miller_array.space_group().is_centric(),
            }
            # Now calculate batch data
            (
                batches,
                rvb,
                isigivb,
                svb,
                batch_data,
            ) = reflection_tables_to_batch_dependent_properties(  # pylint: disable=unbalanced-tuple-unpacking
                scaling_script.reflections,
                scaling_script.experiments,
                scaling_script.scaled_miller_array,
            )
            self.data["scaled_miller_array"] = scaling_script.scaled_miller_array
            self.data["bm"] = batch_manager(batches, batch_data)
            self.data["r_merge_vs_batch"] = rvb
            self.data["scale_vs_batch"] = svb
            self.data["isigivsbatch"] = isigivb
            self.data["image_range_tables"] = [
                make_image_range_table(scaling_script.experiments, self.data["bm"])
            ]

    def make_plots(self):
        """Generate tables of overall and resolution-binned merging statistics."""
        d = {
            "scaling_tables": ([], []),
            "resolution_plots": OrderedDict(),
            "batch_plots": OrderedDict(),
            "misc_plots": OrderedDict(),
            "anom_plots": OrderedDict(),
            "image_range_tables": [],
        }
        if "statistics" in self.data:
            plotter = ResolutionPlotsAndStats(
                self.data["statistics"],
                self.data["anomalous_statistics"],
                is_centric=self.data["is_centric"],
            )
            d["resolution_plots"].update(plotter.make_all_plots())
            d["scaling_tables"] = plotter.statistics_tables()
            d["batch_plots"].update(
                scale_rmerge_vs_batch_plot(
                    self.data["bm"],
                    self.data["r_merge_vs_batch"],
                    self.data["scale_vs_batch"],
                )
            )
            d["batch_plots"].update(
                i_over_sig_i_vs_batch_plot(self.data["bm"], self.data["isigivsbatch"])
            )
            plotter = IntensityStatisticsPlots(
                self.data["scaled_miller_array"], run_xtriage_analysis=False
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
            intensities_anom = self.data["scaled_miller_array"].as_anomalous_array()
            intensities_anom = intensities_anom.map_to_asu().customized_copy(
                info=self.data["scaled_miller_array"].info()
            )
            anom_plotter = AnomalousPlotter(intensities_anom, strong_cutoff=d_min)
            d["anom_plots"].update(anom_plotter.make_plots())
            d["image_range_tables"] = self.data["image_range_tables"]
        return d
