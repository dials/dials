"""
Observers for the scaling algorithm.
"""
from collections import OrderedDict
from dials.util.observer import Observer, singleton
from dials.algorithms.scaling.plots import (
    plot_scaling_models,
    statistics_tables,
    plot_outliers,
    normal_probability_plot,
)
from jinja2 import Environment, ChoiceLoader, PackageLoader


def register_default_scaling_observers(script):
    """Register the standard observers to the scaling script."""
    script.register_observer(
        event="merging_statistics", observer=MergingStatisticsObserver()
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
    script.scaler.register_observer(
        event="performed_error_analysis", observer=ErrorModelObserver()
    )

@singleton
class ScalingHTMLGenerator(Observer):

    """
    Observer to make a html report
    """

    def make_scaling_html(self, scaling_script):
        """Collect data from the individual observers and write the html."""
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
            scaling_outlier_graphs=self.data["outlier_plots"],
            normal_prob_plot=self.data["normal_prob_plot"],
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
            self.data = {"statistics": scaling_script.merging_statistics_result}

    def make_plots(self):
        """Generate tables of overall and resolution-binned merging statistics."""
        d = {"scaling_tables": [[], []]}
        if "statistics" in self.data:
            d["scaling_tables"] = statistics_tables(self.data["statistics"])
        return d
