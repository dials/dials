#!/usr/bin/env python
#
# dials.create_profile_model.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import logging

from libtbx.phil import parse

logger = logging.getLogger("dials.command_line.create_profile_model")

help_message = """

This program computes the profile model from the input reflections. It then
saves a modified models.expt file with the additional profile model
information. Usually this is performed during integration; however, on some
occasions it may be desirable to compute the profile model independently.

Examples::

  dials.create_profile_model models.expt observations.refl

"""

phil_scope = parse(
    """
  subtract_background = False
    .type = bool
    .help = "Subtract background from pixel data before computing profile"
    .expert_level = 2
  output = models_with_profiles.expt
    .type = str
    .help = "The filename for the experiments"

  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
""",
    process_includes=True,
)


class Script(object):
    """ Encapsulate the script in a class. """

    def __init__(self):
        """ Initialise the script. """
        from dials.util.options import OptionParser

        usage = "usage: dials.create_profile_model [options] models.expt spots.refl"
        self.parser = OptionParser(
            usage=usage,
            epilog=help_message,
            phil=phil_scope,
            read_reflections=True,
            read_experiments=True,
            check_format=False,
        )

    def run(self):
        """ Run the script. """
        from dials.algorithms.profile_model.factory import ProfileModelFactory
        from dials.util.command_line import Command
        from dials.array_family import flex
        from dials.util.options import flatten_reflections, flatten_experiments
        from dxtbx.model.experiment_list import ExperimentListDumper
        from dials.util import Sorry
        from dials.util import log

        log.config()

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=True)
        reflections = flatten_reflections(params.input.reflections)
        experiments = flatten_experiments(params.input.experiments)
        if len(reflections) == 0 and len(experiments) == 0:
            self.parser.print_help()
            return
        if len(reflections) != 1:
            raise Sorry("exactly 1 reflection table must be specified")
        if len(experiments) == 0:
            raise Sorry("no experiments were specified")
        if (not "background.mean" in reflections[0]) and params.subtract_background:
            raise Sorry("for subtract_background need background.mean in reflections")

        reflections, _ = self.process_reference(reflections[0], params)

        # Check pixels don't belong to neighbours
        self.filter_reference_pixels(reflections, experiments)

        # Predict the reflections
        logger.info("")
        logger.info("=" * 80)
        logger.info("")
        logger.info("Predicting reflections")
        logger.info("")
        predicted = flex.reflection_table.from_predictions_multi(
            experiments,
            dmin=params.prediction.d_min,
            dmax=params.prediction.d_max,
            margin=params.prediction.margin,
            force_static=params.prediction.force_static,
            padding=params.prediction.padding,
        )

        # Match with predicted
        matched, reflections, unmatched = predicted.match_with_reference(reflections)
        assert len(matched) == len(predicted)
        assert matched.count(True) <= len(reflections)
        if matched.count(True) == 0:
            raise Sorry(
                """
        Invalid input for reference reflections.
        Zero reference spots were matched to predictions
      """
            )
        elif len(unmatched) != 0:
            logger.info("")
            logger.info("*" * 80)
            logger.info(
                "Warning: %d reference spots were not matched to predictions"
                % (len(unmatched))
            )
            logger.info("*" * 80)
            logger.info("")

        # Create the profile model
        experiments = ProfileModelFactory.create(params, experiments, reflections)
        for model in experiments:
            sigma_b = model.profile.sigma_b(deg=True)
            sigma_m = model.profile.sigma_m(deg=True)
            if model.profile.is_scan_varying():  # scan varying
                mean_sigma_b = sum(sigma_b) / len(sigma_b)
                mean_sigma_m = sum(sigma_m) / len(sigma_m)
                logger.info("Sigma B: %f", mean_sigma_b)
                logger.info("Sigma M: %f", mean_sigma_m)
            else:
                logger.info("Sigma B: %f", sigma_b)
                logger.info("Sigma M: %f", sigma_m)

        # Write the parameters
        Command.start("Writing experiments to %s" % params.output)
        dump = ExperimentListDumper(experiments)
        with open(params.output, "w") as outfile:
            outfile.write(dump.as_json())
        Command.end("Wrote experiments to %s" % params.output)

    def process_reference(self, reference, params):
        """ Load the reference spots. """
        from time import time
        from dials.util import Sorry

        if reference is None:
            return None, None
        st = time()
        assert "miller_index" in reference
        assert "id" in reference
        logger.info("Processing reference reflections")
        logger.info(" read %d strong spots" % len(reference))
        mask = reference.get_flags(reference.flags.indexed)
        rubbish = reference.select(mask == False)
        if mask.count(False) > 0:
            reference.del_selected(mask == False)
            logger.info(" removing %d unindexed reflections" % mask.count(False))
        if len(reference) == 0:
            raise Sorry(
                """
        Invalid input for reference reflections.
        Expected > %d indexed spots, got %d
      """
                % (0, len(reference))
            )
        mask = reference.get_flags(reference.flags.centroid_outlier)
        if mask.count(True) > 0:
            rubbish.extend(reference.select(mask))
            reference.del_selected(mask)
            logger.info(
                " removing %d reflections marked as centroid outliers"
                % mask.count(True)
            )
        mask = reference["miller_index"] == (0, 0, 0)
        if mask.count(True) > 0:
            rubbish.extend(reference.select(mask))
            reference.del_selected(mask)
            logger.info(" removing %d reflections with hkl (0,0,0)" % mask.count(True))
        mask = reference["id"] < 0
        if mask.count(True) > 0:
            raise Sorry(
                """
        Invalid input for reference reflections.
        %d reference spots have an invalid experiment id
      """
                % mask.count(True)
            )
        logger.info(" using %d indexed reflections" % len(reference))
        logger.info(" found %d junk reflections" % len(rubbish))

        if "background.mean" in reference and params.subtract_background:
            logger.info(
                " subtracting background from %d reference reflections" % len(reference)
            )
            for spot in reference:
                spot["shoebox"].data -= spot["background.mean"]
        logger.info(" time taken: %g" % (time() - st))
        return reference, rubbish

    def filter_reference_pixels(self, reference, experiments):
        """
        Set any pixel closer to other reflections to background

        """
        modified_count = 0
        for experiment, indices in reference.iterate_experiments_and_indices(
            experiments
        ):
            subset = reference.select(indices)
            modified = subset["shoebox"].mask_neighbouring(
                subset["miller_index"],
                experiment.beam,
                experiment.detector,
                experiment.goniometer,
                experiment.scan,
                experiment.crystal,
            )
            modified_count += modified.count(True)
            reference.set_selected(indices, subset)
        logger.info(" masked neighbouring pixels in %d shoeboxes" % modified_count)
        return reference


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
