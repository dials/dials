from __future__ import annotations

import numpy as np

from dxtbx import flumpy
from dxtbx.model import ExperimentList
from scitbx.array_family import flex

from dials.algorithms.integration.ssx.ssx_integrate import (
    OutputAggregator,
    OutputCollector,
    SimpleIntegrator,
)
from dials.algorithms.profile_model.ellipsoid.algorithm import (
    final_integrator,
    initial_integrator,
    predict,
    run_ellipsoid_refinement,
)
from dials.algorithms.profile_model.ellipsoid.indexer import reindex


class ToFewReflections(Exception):
    pass


class EllipsoidOutputCollector(OutputCollector):
    def collect_after_refinement(self, experiment, reflection_table, refiner_output):
        if "initial_rmsd_x" not in self.data:
            self.data["initial_rmsd_x"] = refiner_output[0][0]["rmsd"][0]
            self.data["initial_rmsd_y"] = refiner_output[0][0]["rmsd"][1]
        self.data["final_rmsd_x"] = refiner_output[0][-1]["rmsd"][0]
        self.data["final_rmsd_y"] = refiner_output[0][-1]["rmsd"][1]
        # also record some model info:
        self.data["profile_model"] = experiment.profile.to_dict()
        self.data["profile_model_mosaicity"] = experiment.profile.mosaicity()

    def collect_after_preprocess(self, experiment, reflection_table):
        self.data["n_strong_after_preprocess"] = reflection_table.size()
        sel = reflection_table.get_flags(reflection_table.flags.integrated_sum)
        n_sum = sel.count(True)
        self.data["n_strong_sum_integrated"] = n_sum
        if n_sum:
            xobs, yobs, _ = reflection_table["xyzobs.px.value"].select(sel).parts()
            xcal, ycal, _ = reflection_table["xyzcal.px"].select(sel).parts()
            rmsd = flex.mean((xobs - xcal) ** 2 + (yobs - ycal) ** 2) ** 0.5
            self.data["strong_rmsd_preprocessed"] = rmsd
        else:
            self.data["strong_rmsd_preprocessed"] = 0.0

    def collect_after_integration(self, experiment, reflection_table):
        super().collect_after_integration(experiment, reflection_table)
        p = flumpy.to_numpy(reflection_table["partiality"])
        bins = np.linspace(0, 1, 11)
        hist = np.histogram(p, bins)
        self.data["partiality"] = hist[0]


class EllipsoidIntegrator(SimpleIntegrator):
    def __init__(self, params, collect_data=False):
        super().__init__(params)
        if collect_data:
            self.collector = EllipsoidOutputCollector()

    def run(self, experiment, table):

        # first set ids to zero so can integrate (this is how integration
        # finds the image in the imageset)
        ids_map = dict(table.experiment_identifiers())
        table["id"] = flex.int(table.size(), 0)
        del table.experiment_identifiers()[list(ids_map.keys())[0]]
        table.experiment_identifiers()[0] = list(ids_map.values())[0]

        fix_list = []
        if self.params.profile.ellipsoid.unit_cell.fixed:
            fix_list.append("unit_cell")
        if self.params.profile.ellipsoid.orientation.fixed:
            fix_list.append("orientation")

        self.collector.initial_collect(experiment, table)

        for _ in range(self.params.profile.ellipsoid.refinement.n_macro_cycles):
            try:
                table, sigma_d = self.preprocess(experiment, table, self.params)
                self.collector.collect_after_preprocess(experiment, table)
            except ToFewReflections as e:
                raise RuntimeError(e)
            else:
                experiment, table, refiner_output = self.refine(
                    experiment,
                    table,
                    sigma_d,
                    profile_model=self.params.profile.ellipsoid.rlp_mosaicity.model,
                    fix_list=fix_list,
                    n_cycles=self.params.profile.ellipsoid.refinement.n_cycles,
                    capture_progress=isinstance(
                        self.collector, EllipsoidOutputCollector
                    ),
                )
                self.collector.collect_after_refinement(
                    experiment, table, refiner_output["refiner_output"]["history"]
                )

        predicted = self.predict(
            experiment,
            table,
            d_min=self.params.prediction.d_min,
            prediction_probability=self.params.profile.ellipsoid.prediction.probability,
        )
        # do we want to add unmatched i.e. strong spots which weren't predicted?
        self.collector.collect_after_prediction(predicted, table)

        predicted = self.integrate(experiment, predicted, sigma_d)
        self.collector.collect_after_integration(experiment, predicted)

        return experiment, predicted, self.collector

    # methods below are staticmethods, so that one can build different
    # workflows from the run method above using individual algorithm components

    @staticmethod
    def preprocess(experiment, reflection_table, params):
        reference = reindex(
            reflection_table,
            experiment,
            outlier_probability=params.profile.ellipsoid.refinement.outlier_probability,
            max_separation=params.profile.ellipsoid.refinement.max_separation,
            fail_on_bad_index=params.profile.ellipsoid.indexing.fail_on_bad_index,
        )
        if reference.size() < params.profile.ellipsoid.refinement.min_n_reflections:
            raise ToFewReflections(
                "Too few reflections to perform refinement: got %d, expected %d"
                % (
                    reference.size(),
                    params.profile.ellipsoid.refinement.min_n_reflections,
                )
            )
        reference, sigma_d = initial_integrator(ExperimentList([experiment]), reference)

        return reference, sigma_d

    @staticmethod
    def refine(
        experiment,
        reflection_table,
        sigma_d,
        profile_model,
        fix_list=None,
        n_cycles=1,
        capture_progress=False,
    ):
        fix_unit_cell = False
        fix_orientation = False
        if fix_list:
            if "unit_cell" in fix_list:
                fix_unit_cell = True
            if "orientation" in fix_list:
                fix_orientation = True

        expts, refls, output_data = run_ellipsoid_refinement(
            ExperimentList([experiment]),
            reflection_table,
            sigma_d,
            profile_model=profile_model,
            fix_unit_cell=fix_unit_cell,
            fix_orientation=fix_orientation,
            n_cycles=n_cycles,
            capture_progress=capture_progress,
        )
        return expts[0], refls, output_data

    @staticmethod
    def predict(experiment, reference, d_min=None, prediction_probability=0.9973):
        id_map = dict(reference.experiment_identifiers())
        reflection_table = predict(
            ExperimentList([experiment]),
            d_min=d_min,
            prediction_probability=prediction_probability,
        )
        ids_ = set(reflection_table["id"])
        assert ids_ == set(id_map.keys()), f"{ids_}, {id_map.keys()}"
        for id_ in ids_:
            reflection_table.experiment_identifiers()[id_] = id_map[id_]
        return reflection_table

    @staticmethod
    def integrate(experiment, reflection_table, sigma_d):
        reflection_table = final_integrator(
            ExperimentList([experiment]),
            reflection_table,
            sigma_d,
        )
        return reflection_table


class EllipsoidOutputAggregator(OutputAggregator):
    def make_plots(self):
        plots = super().make_plots()
        initial_rmsds_x = [d["initial_rmsd_x"] for d in self.data.values()]
        final_rmsds_x = [d["final_rmsd_x"] for d in self.data.values()]
        initial_rmsds_y = [d["initial_rmsd_y"] for d in self.data.values()]
        final_rmsds_y = [d["final_rmsd_y"] for d in self.data.values()]
        n = list(self.data.keys())
        hist = np.zeros((10,))
        for d in self.data.values():
            hist += d["partiality"]
        bins = np.linspace(0.05, 0.95, 10)
        rmsd_plots = {
            "refinement_rmsds": {
                "data": [
                    {
                        "x": n,
                        "y": initial_rmsds_x,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "Initial rmsd_x",
                    },
                    {
                        "x": n,
                        "y": final_rmsds_x,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "Final rmsd_x",
                    },
                    {
                        "x": n,
                        "y": initial_rmsds_y,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "Initial rmsd_y",
                    },
                    {
                        "x": n,
                        "y": final_rmsds_y,
                        "type": "scatter",
                        "mode": "markers",
                        "name": "Final rmsd_y",
                    },
                ],
                "layout": {
                    "title": "Rmsds of integrated reflections per image",
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
        plots.update(rmsd_plots)
        return plots
