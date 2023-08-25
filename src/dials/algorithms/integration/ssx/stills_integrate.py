from __future__ import annotations

import numpy as np

from dxtbx import flumpy
from dxtbx.model import ExperimentList

from dials.algorithms.integration.integrator import (
    Parameters,
    _finalize_stills,
    _initialize_stills,
)
from dials.algorithms.integration.ssx.ssx_integrate import (
    OutputCollector,
    SimpleIntegrator,
)
from dials.algorithms.profile_model.factory import ProfileModelFactory
from dials.algorithms.shoebox import MaskCode
from dials.array_family import flex
from dials.command_line.integrate import process_reference


class StillsOutputCollector(OutputCollector):
    def collect_after_integration(self, experiment, reflection_table):
        super().collect_after_integration(experiment, reflection_table)
        self.data["profile_model"] = experiment.profile.to_dict()

        invD = 1.0 / experiment.crystal.get_domain_size_ang()
        e = experiment.crystal.get_half_mosaicity_deg()
        self.data["profile_model_mosaicity"] = {
            "inv_mosaic_domain_size": invD,
            "angular_half_mosaicity": e,
        }
        p = flumpy.to_numpy(reflection_table["partiality"])
        bins = np.linspace(0, 1, 11)
        hist = np.histogram(p, bins)
        self.data["partiality"] = hist[0]


class StillsIntegrator(SimpleIntegrator):
    def __init__(self, params, collect_data=False):
        super().__init__(params)
        if collect_data:
            self.collector = StillsOutputCollector()

    def run(self, experiment, table):

        # first set ids to zero so can integrate (this is how integration
        # finds the image in the imageset)
        ids_map = dict(table.experiment_identifiers())
        table["id"] = flex.int(table.size(), 0)
        del table.experiment_identifiers()[list(ids_map.keys())[0]]
        table.experiment_identifiers()[0] = list(ids_map.values())[0]

        self.collector.initial_collect(experiment, table)

        # preprocess
        table = self.preprocess(table)
        self.collector.collect_after_preprocess(experiment, table)
        # refine - null

        predicted, elist = self.predict(experiment, table, self.params)
        self.collector.collect_after_prediction(predicted, table)
        integrated, elist = self.integrate(elist, predicted, self.params)
        # NB what about things like kapton correction?
        self.collector.collect_after_integration(elist[0], integrated)

        return elist[0], integrated, self.collector

    @staticmethod
    def preprocess(table):
        reference, _ = process_reference(table)
        return reference

    @staticmethod
    def refine():
        pass

    @staticmethod
    def predict(experiment, table, params):
        elist = ExperimentList([experiment])
        predicted = flex.reflection_table.from_predictions_multi(
            elist,
            dmin=params.prediction.d_min,
            dmax=params.prediction.d_max,
            margin=params.prediction.margin,
            force_static=params.prediction.force_static,
            padding=params.prediction.padding,
        )
        # need to set imageset id?
        elist = ProfileModelFactory.create(params, elist, table)
        predicted.compute_bbox(elist)
        return predicted, elist

    @staticmethod
    def integrate(experiments, table, params):

        _params = Parameters.from_phil(params.integration)
        experiments[0].scan = None
        _initialize_stills(experiments, _params, table)

        table["shoebox"] = flex.shoebox(table["panel"], table["bbox"], allocate=True)
        table.extract_shoeboxes(experiments[0].imageset)

        # From integratorexecutor
        table.is_overloaded(experiments)
        table.compute_mask(experiments)
        # Check for invalid pixels in foreground/background
        table.contains_invalid_pixels()
        sbox = table["shoebox"]
        nvalfg = sbox.count_mask_values(MaskCode.Valid | MaskCode.Foreground)
        nforeg = sbox.count_mask_values(MaskCode.Foreground)
        fraction_valid = nvalfg.as_double() / nforeg.as_double()
        selection = (
            fraction_valid < params.integration.profile.valid_foreground_threshold
        )
        table.set_flags(selection, table.flags.dont_integrate)

        # Process the data
        table.compute_background(experiments)
        table.compute_centroid(experiments)
        table.compute_summed_intensity()
        table["num_pixels.valid"] = sbox.count_mask_values(MaskCode.Valid)
        table["num_pixels.background"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background
        )
        table["num_pixels.background_used"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
        )
        table["num_pixels.foreground"] = nvalfg
        table, experiments = _finalize_stills(table, experiments, _params)
        return table, experiments
