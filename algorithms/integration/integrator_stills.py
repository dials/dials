#
# integrator_stills.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function


class ReflectionBlockIntegratorStills(object):
    """ A class to perform the integration. """

    def __init__(self, params, experiments, reference, extractor=None):
        """ Initialise the integrator. """
        from dials.algorithms import shoebox

        # Ensure we have 1 experiment at the moment
        assert len(experiments) == 1
        assert extractor is not None

        # Save the parameters
        self.params = params
        self.experiments = experiments
        self.extractor = extractor

        # Create the shoebox masker
        n_sigma = params.integration.shoebox.n_sigma
        assert n_sigma > 0
        self._mask_profiles = shoebox.MaskerEmpirical(
            experiments[0], reference=reference
        )

    def integrate(self):
        """ Integrate all the reflections. """
        from dials.array_family import flex
        from dials.algorithms.shoebox import MaskCode

        result = flex.reflection_table()
        for indices, reflections in self.extractor:
            self._mask_profiles(reflections, None)
            reflections.integrate(self.experiments[0])
            bg_code = MaskCode.Valid | MaskCode.BackgroundUsed
            fg_code = MaskCode.Valid | MaskCode.Foreground
            n_bg = reflections["shoebox"].count_mask_values(bg_code)
            n_fg = reflections["shoebox"].count_mask_values(fg_code)
            reflections["n_background"] = n_bg
            reflections["n_foreground"] = n_fg
            del reflections["shoebox"]
            del reflections["rs_shoebox"]
            result.extend(reflections)
        assert len(result) > 0
        result.sort("miller_index")
        return result


class IntegratorStills(object):
    """ Integrate reflections """

    def __init__(self, params, exlist, reference=None, predicted=None, shoeboxes=None):
        """Initialise the script."""

        assert reference is not None

        # Load the extractor based on the input
        if shoeboxes is not None:
            extractor = self._load_extractor(shoeboxes, params, exlist)
        else:
            if predicted is None:
                predicted = self._predict_reflections(params, exlist)
                # predicted = self._filter_reflections(params, exlist, predicted) # FIXME

            predicted = self._match_with_reference(predicted, reference)

            from annlib_ext import AnnAdaptor
            from dials.array_family import flex
            import math

            matcheddata = predicted.select(
                predicted.get_flags(predicted.flags.reference_spot)
            )

            A = AnnAdaptor(matcheddata["xyzcal.mm"].as_double(), 3, 10)
            A.query(predicted["xyzcal.mm"].as_double())

            bboxes = flex.int6()
            for i, ref in enumerate(predicted):
                nn_pred = [matcheddata[A.nn[i * 10 + j]] for j in range(10)]
                nn_ref = [
                    reference[reference["miller_index"].first_index(r["miller_index"])]
                    for r in nn_pred
                ]

                max_x = max([r["bbox"][1] - r["bbox"][0] for r in nn_ref])
                max_y = max([r["bbox"][3] - r["bbox"][2] for r in nn_ref])

                panel = exlist[ref["id"]].detector[ref["panel"]]
                imgsize_x, imgsize_y = panel.get_image_size()

                x1 = int(math.floor(ref["xyzcal.px"][0] - (max_x / 2)))
                x2 = int(math.ceil(ref["xyzcal.px"][0] + (max_x / 2)))
                y1 = int(math.floor(ref["xyzcal.px"][1] - (max_y / 2)))
                y2 = int(math.ceil(ref["xyzcal.px"][1] + (max_y / 2)))

                if x1 < 0:
                    x1 = 0
                if y1 < 0:
                    y1 = 0

                if x2 > imgsize_x:
                    x2 = imgsize_x
                if y2 > imgsize_y:
                    y2 = imgsize_y

                bboxes.append((x1, x2, y1, y2, 0, 1))

            predicted["bbox"] = bboxes

            extractor = self._create_extractor(params, exlist, predicted)

        # Initialise the integrator
        self._integrator = ReflectionBlockIntegratorStills(
            params, exlist, reference, extractor
        )

    def integrate(self):
        """ Integrate the reflections. """
        return self._integrator.integrate()

    def _match_with_reference(self, predicted, reference):
        """ Match predictions with reference spots. """

        from dials.algorithms.spot_finding.spot_matcher import SpotMatcher
        from dials.util.command_line import Command

        Command.start("Matching reference spots with predicted reflections")
        match = SpotMatcher(max_separation=1)
        rind, pind = match(reference, predicted)
        h1 = predicted.select(pind)["miller_index"]
        h2 = reference.select(rind)["miller_index"]
        mask = h1 == h2
        predicted.set_flags(pind.select(mask), predicted.flags.reference_spot)
        Command.end(
            "Matched %d reference spots with predicted reflections" % mask.count(True)
        )
        return predicted

    def _load_extractor(self, filename, params, exlist):
        """ Load the shoebox extractor. """
        from dials.model.serialize.reflection_block import ReflectionBlockExtractor

        assert len(exlist) == 1
        imageset = exlist[0].imageset
        return ReflectionBlockExtractor(
            filename, params.integration.shoebox.block_size, imageset
        )

    def _create_extractor(self, params, exlist, predicted):
        """ Create the extractor. """
        from dials.model.serialize.reflection_block import ReflectionBlockExtractor

        assert len(exlist) == 1
        imageset = exlist[0].imageset
        return ReflectionBlockExtractor(
            "shoebox.dat", params.integration.shoebox.block_size, imageset, predicted
        )

    def _predict_reflections(self, params, experiments):
        """ Predict all the reflections. """
        from dials.array_family import flex

        result = flex.reflection_table()
        for i, experiment in enumerate(experiments):
            predicted = flex.reflection_table.from_predictions(experiment)
            predicted["id"] = flex.int(len(predicted), i)
            result.extend(predicted)
        return result

    def _filter_reflections(self, params, experiments, reflections):
        """ Filter the reflections to integrate. """
        from dials.util.command_line import Command
        from dials.algorithms import filtering
        from dials.array_family import flex

        # Set all reflections which overlap bad pixels to zero
        Command.start("Filtering reflections by detector mask")
        if experiments[0].scan is None:
            array_range = 1
        else:
            array_range = experiments[0].scan.get_array_range()
        mask = filtering.by_detector_mask(
            reflections["bbox"],
            experiments[0].imageset.get_raw_data(0)[0] >= 0,
            array_range,
        )
        reflections.del_selected(not mask)
        Command.end("Filtered %d reflections by detector mask" % len(reflections))

        # Filter the reflections by zeta
        min_zeta = params.integration.filter.by_zeta
        if min_zeta > 0:
            Command.start("Filtering reflections by zeta >= %f" % min_zeta)
            zeta = reflections.compute_zeta(experiments[0])
            reflections.del_selected(flex.abs(zeta) < min_zeta)
            n = len(reflections)
            Command.end("Filtered %d reflections by zeta >= %f" % (n, min_zeta))
            return reflections
