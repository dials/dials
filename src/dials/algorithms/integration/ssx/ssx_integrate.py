from __future__ import annotations

import functools
from abc import ABC, abstractmethod

import numpy as np
from jinja2 import ChoiceLoader, Environment, PackageLoader

import dials.extensions
from dials.array_family import flex


def generate_html_report(plots_data, filename):
    loader = ChoiceLoader(
        [
            PackageLoader("dials", "templates"),
            PackageLoader("dials", "static", encoding="utf-8"),
        ]
    )
    env = Environment(loader=loader)
    template = env.get_template("simple_report.html")
    html = template.render(
        page_title="DIALS SSX integration report",
        panel_title="Integration plots",
        graphs=plots_data,
    )
    with open(filename, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))


class SimpleIntegrator(ABC):

    """Define an interface for ssx prediction/integration processing"""

    def __init__(self, params):
        self.params = params
        self.collector = NullCollector()
        BackgroundAlgorithm = dials.extensions.Background.load(
            params.integration.background.algorithm
        )
        flex.reflection_table.background_algorithm = functools.partial(
            BackgroundAlgorithm, params
        )
        CentroidAlgorithm = dials.extensions.Centroid.load(
            params.integration.centroid.algorithm
        )
        flex.reflection_table.centroid_algorithm = functools.partial(
            CentroidAlgorithm, params
        )

    @abstractmethod
    def run(self, experiment, table):
        # all gathering/collecting of output data should be optionally done at
        # the run level, so that the calls to individual processing steps are as
        # fast as possible.
        # In general, most output/statistics are calculable on the experiment
        # or reflection table.
        # However the refine step could be expected to return a json with
        # relevant history
        pass

    @abstractmethod
    def preprocess(experiment, reflection_table, *args, **kwargs):
        pass

    @abstractmethod
    def refine(experiment, reflection_table, *args, **kwargs):
        pass

    @abstractmethod
    def predict(experiment, reflection_table, *args, **kwargs):
        pass

    @abstractmethod
    def integrate(experiment, reflection_table, *args, **kwargs):
        pass


class NullCollector(object):

    """
    Defines a null data collector for cases where you don't want
    to record data during the process.
    """

    def __init__(self):
        self.data = {}

    def initial_collect(self, *args, **kwargs):
        pass

    def collect_after_preprocess(self, *args, **kwargs):
        pass

    def collect_after_refinement(self, *args, **kwargs):
        pass

    def collect_after_prediction(self, *args, **kwargs):
        pass

    def collect_after_integration(self, *args, **kwargs):
        pass


class OutputCollector(object):

    """
    Defines a data collector to log common quantities for all algorithm choices
    for an individual image.
    """

    def __init__(self):
        self.data = {}

    # collects general output for reporting, independent of underlying models,
    # for integration of a single image.

    def initial_collect(self, experiment, reflection_table):
        self.data["initial_n_refl"] = reflection_table.size()
        xobs, yobs, _ = reflection_table["xyzobs.px.value"].parts()
        xcal, ycal, _ = reflection_table["xyzcal.px"].parts()
        rmsd = flex.mean((xobs - xcal) ** 2 + (yobs - ycal) ** 2) ** 0.5
        self.data["strong_rmsd"] = rmsd

    def collect_after_preprocess(self, experiment, reflection_table):
        xobs, yobs, _ = reflection_table["xyzobs.px.value"].parts()
        xcal, ycal, _ = reflection_table["xyzcal.px"].parts()
        rmsd = flex.mean((xobs - xcal) ** 2 + (yobs - ycal) ** 2) ** 0.5
        self.data["strong_rmsd_preprocessed"] = rmsd

    def collect_after_prediction(self, predicted, reference):
        matched, _, unmatched = predicted.match_with_reference(reference)
        self.data["n_strong_predicted"] = matched.count(True)
        self.data["n_strong_unpredicted"] = unmatched.size()

    def collect_after_integration(self, experiment, reflection_table):
        sel = reflection_table.get_flags(reflection_table.flags.integrated_sum)
        n_sum = sel.count(True)
        self.data["n_integrated"] = n_sum
        self.data["n_failed"] = reflection_table.size() - n_sum
        I = reflection_table["intensity.sum.value"].select(sel)
        var = reflection_table["intensity.sum.variance"].select(sel)
        if not var.all_gt(0):
            sel2 = var > 0
            I = I.select(sel2)
            var = var.select(sel2)
        self.data["i_over_sigma_overall"] = flex.mean(I / flex.sqrt(var))


class OutputAggregator:

    """
    Simple aggregator class to aggregate data from all images and generate
    json data for output/plotting.
    """

    def __init__(self):
        self.data = {}

    def add_dataset(self, collector, index):
        if collector.data:
            self.data[index] = collector.data

    def make_history_json(self):
        if not self.data:
            return {}
        history = {}
        for i, d in enumerate(self.data.values()):
            history[i] = {}
            if "likelihood" in d:
                history[i]["likelihood_per_iteration"] = d["likelihood"]
            if "parameters" in d:
                history[i]["active_parameters_per_iteration"] = d["parameters"]
        return history

    def make_plots(self):
        # just make some simple plots for now as a test
        if not self.data:
            return {}
        I_over_sigma = [d["i_over_sigma_overall"] for d in self.data.values()]
        n = list(self.data.keys())
        n_integrated = [d["n_integrated"] for d in self.data.values()]
        n_predicted = [d["n_strong_predicted"] for d in self.data.values()]
        n_strong = [
            (d["n_strong_predicted"] + d["n_strong_unpredicted"])
            for d in self.data.values()
        ]
        overall_rmsd = [d["strong_rmsd"] for d in self.data.values()]
        overall_rmsd_preprocessed = [
            d["strong_rmsd_preprocessed"] for d in self.data.values()
        ]

        hist = np.zeros((10,))
        for d in self.data.values():
            hist += d["partiality"]
        bins = np.linspace(0.05, 0.95, 10)

        plots_dict = {
            "I_over_sigma_overall": {
                "data": [
                    (
                        {
                            "x": n,
                            "y": I_over_sigma,
                            "type": "scatter",
                            "mode": "markers",
                        }
                    )
                ],
                "layout": {
                    "title": "Overall I/sigma per image",
                    "xaxis": {"title": "image number"},
                    "yaxis": {"title": "I/sigma"},
                },
            },
            "n_predicted": {
                "data": [
                    {
                        "x": n,
                        "y": n_predicted,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "Number of indexed spots predicted",
                    },
                    {
                        "x": n,
                        "y": n_strong,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "Total number of indexed spots",
                    },
                ],
                "layout": {
                    "title": "Number of predicted reflections per image",
                    "xaxis": {"title": "image number"},
                    "yaxis": {"title": "N. reflections"},
                },
            },
            "n_integrated": {
                "data": [
                    (
                        {
                            "x": n,
                            "y": n_integrated,
                            "type": "scatter",
                            "mode": "markers",
                        }
                    )
                ],
                "layout": {
                    "title": "Number of integrated reflections per image",
                    "xaxis": {"title": "image number"},
                    "yaxis": {"title": "N. reflections"},
                },
            },
            "overall_rmsds": {
                "data": [
                    {
                        "x": n,
                        "y": overall_rmsd,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "Initial 2D rmsd of strong spots",
                    },
                    {
                        "x": n,
                        "y": overall_rmsd_preprocessed,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "2D rmsd after preprocess",
                    },
                ],
                "layout": {
                    "title": "Rmsds of strong (indexed) reflections per image",
                    "xaxis": {"title": "image number"},
                    "yaxis": {"title": "RMSD (px)"},
                },
            },
            "partiality": {
                "data": [
                    (
                        {
                            "x": list(bins),
                            "y": list(hist),
                            "type": "scatter",
                            "mode": "markers",
                        }
                    )
                ],
                "layout": {
                    "title": "Partiality distribution",
                    "xaxis": {"title": "Partiality bin centre"},
                    "yaxis": {"title": "Number of reflections"},
                },
            },
        }

        value = list(self.data.values())[0]
        mosaicities = {}
        for k in value["profile_model_mosaicity"].keys():
            mosaicities["M_" + k] = np.zeros(shape=(len(self.data),))
        for i, value in enumerate(self.data.values()):
            sigma = value["profile_model_mosaicity"]
            for k, v in sigma.items():
                mosaicities["M_" + k][i] = v
        data = []
        data_angular = []
        for k, v in mosaicities.items():
            if "angular" in k:
                data_angular.append(
                    {
                        "x": n,
                        "y": list(v),
                        "type": "scatter",
                        "mode": "markers",
                        "name": k,
                        "yaxis": "y2",
                    }
                )
            else:
                data.append(
                    {
                        "x": n,
                        "y": list(v),
                        "type": "scatter",
                        "mode": "markers",
                        "name": k,
                    }
                )
        mosaic_plots = {
            "mosaicities": {
                "data": data,
                "layout": {
                    "title": "Profile model mosaicities per image",
                    "xaxis": {"title": "image number"},
                    "yaxis": {"title": "Invariant crystal mosaicity (Å⁻¹)"},
                },
            },
        }
        if data_angular:
            mosaic_plots["mosaicities"]["layout"]["yaxis2"] = {
                "anchor": "x",
                "side": "right",
                "title": "Angular mosaicity (degrees)",
                "overlaying": "y",
            }
            mosaic_plots["mosaicities"]["data"].extend(data_angular)
        plots_dict.update(mosaic_plots)

        return plots_dict
