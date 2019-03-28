"""
Observers for the scaling algorithm.
"""
from collections import OrderedDict
import logging
from scitbx.array_family import flex
from dials.util.observer import Observer, singleton
from dials.algorithms.scaling.plots import (
    plot_scaling_models,
    plot_outliers,
    normal_probability_plot,
)
from dials.report.analysis import reflection_tables_to_batch_dependent_properties
from dials.report.plots import (
    scale_rmerge_vs_batch_plot,
    i_over_sig_i_vs_batch_plot,
    statistics_tables,
    cc_one_half_plot,
)
from dials.util.batch_handling import batch_manager

from jinja2 import Environment, ChoiceLoader, PackageLoader

logger = logging.getLogger("dials")


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

    script.scaler.register_observer(
        event="performed_scaling", observer=ScalingModelObserver()
    )
    script.scaler.register_observer(
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


@singleton
class ScalingSummaryGenerator(Observer):
    """
    Observer to summarise data
    """

    def print_scaling_summary(self, scaling_script):
        if ScalingModelObserver().data:
            logger.info(ScalingModelObserver().return_model_error_summary())
        if MergingStatisticsObserver().data:
            logger.info(
                "\n\t----------Overall merging statistics (non-anomalous)----------\t\n"
            )
            logger.info(MergingStatisticsObserver().make_statistics_summary())


@singleton
class ScalingHTMLGenerator(Observer):

    """
    Observer to make a html report
    """

    def make_scaling_html(self, scaling_script):
        """Collect data from the individual observers and write the html."""
        if not scaling_script.params.output.html:
            return
        self.data.update(ScalingModelObserver().make_plots())
        self.data.update(ScalingOutlierObserver().make_plots())
        self.data.update(ErrorModelObserver().make_plots())
        self.data.update(MergingStatisticsObserver().make_plots())
        filename = scaling_script.params.output.html
        print("Writing html report to: %s" % filename)
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        env = Environment(loader=loader)
        template = env.get_template("scaling_report.html")
        html = template.render(
            page_title="DIALS scaling report",
            scaling_model_graphs=self.data["scaling_model"],
            scaling_tables=self.data["scaling_tables"],
            cc_one_half_plot=self.data["cc_one_half_plot"],
            scaling_outlier_graphs=self.data["outlier_plots"],
            normal_prob_plot=self.data["normal_prob_plot"],
            batch_plots=self.data["batch_plots"],
        )
        with open(filename, "wb") as f:
            f.write(html.encode("ascii", "xmlcharrefreplace"))


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
        for s in active_scalers:
            id_ = s.experiment.identifier
            self.data[id_] = s.experiment.scaling_model.to_dict()

    def make_plots(self):
        """Generate scaling model component plot data."""
        d = OrderedDict()
        for key in sorted(self.data.keys()):
            scaling_model_plots = plot_scaling_models(self.data[key])
            for name, plot in scaling_model_plots.iteritems():
                d.update({name + "_" + str(key): plot})
        graphs = {"scaling_model": d}
        return graphs

    def return_model_error_summary(self):
        """Get a summary of the error distribution of the models."""
        first_model = self.data.values()[0]
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
        active_scalers = getattr(scaler, "active_scalers", False)
        if not active_scalers:
            active_scalers = [scaler]
        for scaler in active_scalers:
            id_ = scaler.experiment.identifier
            outlier_isel = scaler.suitable_refl_for_scaling_sel.iselection().select(
                scaler.outliers
            )
            x, y, z = (
                scaler.reflection_table["xyzobs.px.value"].select(outlier_isel).parts()
            )
            self.data[id_] = {
                "x": list(x),
                "y": list(y),
                "z": list(z),
                "image_size": scaler.experiment.detector[0].get_image_size(),
                "z_range": [
                    i / scaler.experiment.scan.get_oscillation()[1]
                    for i in scaler.experiment.scan.get_oscillation_range()
                ],
            }

    def make_plots(self):
        """Generate plot data of outliers on the detector and vs z."""
        d = OrderedDict()
        for key in sorted(self.data.keys()):
            outlier_plots = plot_outliers(self.data[key])
            d.update(
                {"outlier_plot_" + str(key): outlier_plots["outlier_xy_positions"]}
            )
            d.update({"outlier_plot_z" + str(key): outlier_plots["outliers_vs_z"]})
        graphs = {"outlier_plots": d}
        return graphs


@singleton
class ErrorModelObserver(Observer):

    """
    Observer to record scaling error model data and make a plot.
    """

    def update(self, scaler):
        self.data["delta_hl"] = list(
            scaler.experiment.scaling_model.error_model.delta_hl
        )

    def make_plots(self):
        """Generate normal probability plot data."""
        d = {"normal_prob_plot": {}}
        if "delta_hl" in self.data:
            d["normal_prob_plot"] = normal_probability_plot(self.data)
        return d


@singleton
class MergingStatisticsObserver(Observer):

    """
    Observer to record merging statistics data and make tables.
    """

    def update(self, scaling_script):
        if scaling_script.merging_statistics_result:
            self.data = {
                "statistics": scaling_script.merging_statistics_result,
                "is_centric": scaling_script.scaled_miller_array.space_group().is_centric(),
            }
            # Now calculate batch data
            batches, rvb, isigivb, svb, batch_data = reflection_tables_to_batch_dependent_properties(  # pylint: disable=unbalanced-tuple-unpacking
                scaling_script.reflections,
                scaling_script.experiments,
                scaling_script.scaled_miller_array,
            )

            self.data["bm"] = batch_manager(batches, batch_data)
            self.data["r_merge_vs_batch"] = rvb
            self.data["scale_vs_batch"] = svb
            self.data["isigivsbatch"] = isigivb

    def make_plots(self):
        """Generate tables of overall and resolution-binned merging statistics."""
        d = {
            "scaling_tables": [[], []],
            "cc_one_half_plot": {},
            "batch_plots": OrderedDict(),
        }
        if "statistics" in self.data:
            d["scaling_tables"] = statistics_tables(self.data["statistics"])
            d["cc_one_half_plot"] = cc_one_half_plot(
                self.data["statistics"], is_centric=self.data["is_centric"]
            )
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
        return d

    def make_statistics_summary(self):
        """Format merging statistics information into an output string."""
        result = self.data["statistics"]
        overall = result.overall
        # First make overall summary
        msg = (
            "Resolution: {0:.2f} - {1:.2f} {sep}Observations: {2:d} {sep}"
            "Unique reflections: {3:d} {sep}Redundancy: {4:.1f} {sep}"
            "Completeness: {5:.2f}% {sep}Mean intensity: {6:.1f} {sep}"
            "Mean I/sigma(I): {7:.1f}"
        ).format(
            overall.d_max,
            overall.d_min,
            overall.n_obs,
            overall.n_uniq,
            overall.mean_redundancy,
            overall.completeness * 100,
            overall.i_mean,
            overall.i_over_sigma_mean,
            sep="\n",
        )
        if overall.n_neg_sigmas > 0:
            msg += "SigI < 0 (rejected): {0} observations\n".format(
                overall.n_neg_sigmas
            )
        if overall.n_rejected_before_merge > 0:
            msg += "I < -3*SigI (rejected): {0} observations\n".format(
                overall.n_rejected_before_merge
            )
        if overall.n_rejected_after_merge > 0:
            msg += "I < -3*SigI (rejected): {0} reflections\n".format(
                overall.n_rejected_after_merge
            )
        msg += "{sep}R-merge: {0:5.3f}{sep}R-meas:  {1:5.3f}{sep}R-pim:   {2:5.3f}{sep}".format(
            overall.r_merge, overall.r_meas, overall.r_pim, sep="\n"
        )

        # Next make statistics by resolution bin
        msg += "\nStatistics by resolution bin:\n"
        msg += (
            " d_max  d_min   #obs  #uniq   mult.  %comp       <I>  <I/sI>"
            + "    r_mrg   r_meas    r_pim   cc1/2   cc_ano\n"
        )
        for bin_stats in self.data["statistics"].bins:
            msg += bin_stats.format() + "\n"
        msg += self.data["statistics"].overall.format() + "\n"

        # Now show estimated cutoffs, based on iotbx code
        def format_d_min(value):
            """Format result values"""
            if value is None:
                return "(use all data)"
            return "%7.3f" % value

        msg += "\nResolution cutoff estimates:"
        msg += """
resolution of all data          : {0:7.3f}
based on CC(1/2) >= 0.33        : {1}
based on mean(I/sigma) >= 2.0   : {2}
based on R-merge < 0.5          : {3}
based on R-meas < 0.5           : {4}
based on completeness >= 90%    : {5}
based on completeness >= 50%    : {6}\n""".format(
            result.overall.d_min,
            format_d_min(result.estimate_d_min(min_cc_one_half=0.33)),
            format_d_min(result.estimate_d_min(min_i_over_sigma=2.0)),
            format_d_min(result.estimate_d_min(max_r_merge=0.5)),
            format_d_min(result.estimate_d_min(max_r_meas=0.5)),
            format_d_min(result.estimate_d_min(min_completeness=0.9)),
            format_d_min(result.estimate_d_min(min_completeness=0.5)),
        )
        msg += "NOTE: we recommend using all data out to the CC(1/2) limit for refinement\n"
        return msg
