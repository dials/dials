#!/usr/bin/env python
#
# dials.integrate.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import logging

from dials.array_family import flex
from dials.util import Sorry

logger = logging.getLogger("dials.command_line.integrate")
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

help_message = """

This program is used to integrate the reflections on the diffraction images. It
is called with an experiment list outputted from dials.index or dials.refine and
a corresponding set of strong spots from which a profile model is calculated.
The program will output a set of integrated reflections and an experiment list
with additional profile model data. The data can be reintegrated using the same
profile model by inputting this integrated.expt file back into
dials.integate.

Examples::

  dials.integrate models.expt refined.refl

  dials.integrate models.expt refined.refl output.reflections=integrated.refl

  dials.integrate models.expt refined.refl profile.fitting=False

  dials.integrate models.expt refined.refl background.algorithm=glm

"""

# Create the phil scope
from libtbx.phil import parse

phil_scope = parse(
    """

  output {
    experiments = 'integrated.expt'
      .type = str
      .help = "The experiments output filename"

    reflections = 'integrated.refl'
      .type = str
      .help = "The integrated output filename"

    phil = 'dials.integrate.phil'
      .type = str
      .help = "The output phil file"

    log = 'dials.integrate.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.integrate.debug.log'
      .type = str
      .help = "The debug log filename"

    report = None
      .type = str
      .help = "The integration report filename (*.xml or *.json)"

    include_bad_reference = False
      .type = bool
      .help = "Include bad reference data including unindexed spots,"
              "and reflections whose predictions are messed up in the"
              "reflection table output. Reflections will have the"
              "'bad_reference' flag set."
  }

  scan_range = None
    .type = ints(size=2)
    .help = "Explicitly specify the images to be processed. Only applicable"
            "when experiment list contains a single imageset."
    .multiple = True

  create_profile_model = True
    .type = bool
    .help = "Create the profile model"

  sampling
    .expert_level = 1
  {

    reflections_per_degree = 50
      .help = "The number of predicted reflections per degree of the sweep "
              "to integrate."
      .type = float(value_min=0.)

    minimum_sample_size = 1000
      .help = "cutoff that determines whether subsetting of the input "
              "prediction list is done"
      .type = int

    maximum_sample_size = None
      .help = "The maximum number of predictions to integrate."
              "Overrides reflections_per_degree if that produces a"
              "larger sample size."
      .type = int(value_min=1)

    integrate_all_reflections = True
      .help = "Override reflections_per_degree and integrate all predicted"
              "reflections."
      .type = bool

  }

  exclude_images = None
    .type = ints
    .help = "Exclude images from integration (e.g. 1,2,3,4,5 etc)"

  verbosity = 0
    .type = int(value_min=0)
    .help = "The verbosity level"

  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
  include scope dials.algorithms.integration.stills_significance_filter.phil_scope
  include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope

""",
    process_includes=True,
)


class Script(object):
    """ The integration program. """

    def __init__(self, phil=phil_scope):
        """Initialise the script."""
        from dials.util.options import OptionParser

        # The script usage
        usage = "usage: dials.integrate [options] models.expt"

        # Create the parser
        self.parser = OptionParser(
            usage=usage,
            phil=phil,
            epilog=help_message,
            read_experiments=True,
            read_reflections=True,
        )

    def run(self, args=None):
        """ Perform the integration. """
        from dials.util.command_line import heading
        from dials.util.options import flatten_reflections, flatten_experiments
        from dials.util import log
        from time import time
        from dials.util import Sorry

        # Check the number of arguments is correct
        start_time = time()

        # Parse the command line
        params, options = self.parser.parse_args(args=args, show_diff_phil=False)
        reference = flatten_reflections(params.input.reflections)
        experiments = flatten_experiments(params.input.experiments)
        if len(reference) == 0 and len(experiments) == 0:
            self.parser.print_help()
            return
        if len(reference) == 0:
            reference = None
        elif len(reference) != 1:
            raise Sorry("more than 1 reflection file was given")
        else:
            reference = reference[0]
        if len(experiments) == 0:
            raise Sorry("no experiment list was specified")

        # Save phil parameters
        if params.output.phil is not None:
            with open(params.output.phil, "w") as outfile:
                outfile.write(self.parser.diff_phil.as_str())

        if __name__ == "__main__":
            # Configure logging
            log.config(
                params.verbosity, info=params.output.log, debug=params.output.debug_log
            )

        from dials.util.version import dials_version

        logger.info(dials_version())

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        for abs_params in params.absorption_correction:
            if abs_params.apply:
                if not (
                    params.integration.debug.output
                    and not params.integration.debug.separate_files
                ):
                    raise Sorry(
                        "Shoeboxes must be saved to integration intermediates to apply an absorption correction. "
                        + "Set integration.debug.output=True, integration.debug.separate_files=False and "
                        + "integration.debug.delete_shoeboxes=True to temporarily store shoeboxes."
                    )

        # Print if we're using a mask
        for i, exp in enumerate(experiments):
            mask = exp.imageset.external_lookup.mask
            if mask.filename is not None:
                if mask.data:
                    logger.info("Using external mask: %s" % mask.filename)
                    for tile in mask.data:
                        logger.info(
                            " Mask has %d pixels masked" % tile.data().count(False)
                        )

        # Print the experimental models
        for i, exp in enumerate(experiments):
            logger.info("=" * 80)
            logger.info("")
            logger.info("Experiments")
            logger.info("")
            logger.info("Models for experiment %d" % i)
            logger.info("")
            logger.info(str(exp.beam))
            logger.info(str(exp.detector))
            if exp.goniometer:
                logger.info(str(exp.goniometer))
            if exp.scan:
                logger.info(str(exp.scan))
            logger.info(str(exp.crystal))

        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Initialising"))
        logger.info("")

        # Load the data
        reference, rubbish = self.process_reference(reference)

        # Check pixels don't belong to neighbours
        if reference is not None:
            if exp.goniometer is not None and exp.scan is not None:
                self.filter_reference_pixels(reference, experiments)
        logger.info("")

        # Initialise the integrator
        from dials.algorithms.profile_model.factory import ProfileModelFactory
        from dials.algorithms.integration.integrator import IntegratorFactory

        # Modify experiment list if scan range is set.
        experiments, reference = self.split_for_scan_range(
            experiments, reference, params.scan_range
        )

        # Modify experiment list if exclude images is set
        experiments = self.exclude_images(experiments, params.exclude_images)

        # Predict the reflections
        logger.info("")
        logger.info("=" * 80)
        logger.info("")
        logger.info(heading("Predicting reflections"))
        logger.info("")
        predicted = flex.reflection_table.from_predictions_multi(
            experiments,
            dmin=params.prediction.d_min,
            dmax=params.prediction.d_max,
            margin=params.prediction.margin,
            force_static=params.prediction.force_static,
            padding=params.prediction.padding,
        )

        # Match reference with predicted
        if reference:
            matched, reference, unmatched = predicted.match_with_reference(reference)
            assert len(matched) == len(predicted)
            assert matched.count(True) <= len(reference)
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
            rubbish.extend(unmatched)

            if len(experiments) > 1:
                # filter out any experiments without matched reference reflections
                # f_: filtered
                from dxtbx.model.experiment_list import ExperimentList

                f_reference = flex.reflection_table()
                f_predicted = flex.reflection_table()
                f_rubbish = flex.reflection_table()
                f_experiments = ExperimentList()
                good_expt_count = 0

                def refl_extend(src, dest, eid):
                    tmp = src.select(src["id"] == eid)
                    tmp["id"] = flex.int(len(tmp), good_expt_count)
                    dest.extend(tmp)

                for expt_id, experiment in enumerate(experiments):
                    if len(reference.select(reference["id"] == expt_id)) != 0:
                        refl_extend(reference, f_reference, expt_id)
                        refl_extend(predicted, f_predicted, expt_id)
                        refl_extend(rubbish, f_rubbish, expt_id)
                        f_experiments.append(experiment)
                        good_expt_count += 1
                    else:
                        logger.info(
                            "Removing experiment %d: no reference reflections matched to predictions"
                            % expt_id
                        )

                reference = f_reference
                predicted = f_predicted
                experiments = f_experiments
                rubbish = f_rubbish

        # Select a random sample of the predicted reflections
        if not params.sampling.integrate_all_reflections:
            predicted = self.sample_predictions(experiments, predicted, params)

        # Compute the profile model
        if (
            params.create_profile_model
            and reference is not None
            and "shoebox" in reference
        ):
            experiments = ProfileModelFactory.create(params, experiments, reference)
        else:
            experiments = ProfileModelFactory.create(params, experiments)
            for expr in experiments:
                if expr.profile is None:
                    raise Sorry("No profile information in experiment list")
        del reference

        # Compute the bounding box
        predicted.compute_bbox(experiments)

        # Create the integrator
        logger.info("")
        integrator = IntegratorFactory.create(params, experiments, predicted)

        # Integrate the reflections
        reflections = integrator.integrate()

        # Append rubbish data onto the end
        if rubbish is not None and params.output.include_bad_reference:
            mask = flex.bool(len(rubbish), True)
            rubbish.unset_flags(mask, rubbish.flags.integrated_sum)
            rubbish.unset_flags(mask, rubbish.flags.integrated_prf)
            rubbish.set_flags(mask, rubbish.flags.bad_reference)
            reflections.extend(rubbish)

        # Correct integrated intensities for absorption correction, if necessary
        for abs_params in params.absorption_correction:
            if abs_params.apply and abs_params.algorithm == "fuller_kapton":
                from dials.algorithms.integration.kapton_correction import (
                    multi_kapton_correction,
                )

                experiments, reflections = multi_kapton_correction(
                    experiments, reflections, abs_params.fuller_kapton, logger=logger
                )()

        if params.significance_filter.enable:
            from dials.algorithms.integration.stills_significance_filter import (
                SignificanceFilter,
            )
            from dxtbx.model.experiment_list import ExperimentList

            sig_filter = SignificanceFilter(params)
            filtered_refls = sig_filter(experiments, reflections)
            accepted_expts = ExperimentList()
            accepted_refls = flex.reflection_table()
            logger.info(
                "Removed %d reflections out of %d when applying significance filter"
                % (len(reflections) - len(filtered_refls), len(reflections))
            )
            for expt_id, expt in enumerate(experiments):
                refls = filtered_refls.select(filtered_refls["id"] == expt_id)
                if len(refls) > 0:
                    accepted_expts.append(expt)
                    refls["id"] = flex.int(len(refls), len(accepted_expts) - 1)
                    accepted_refls.extend(refls)
                else:
                    logger.info(
                        "Removed experiment %d which has no reflections left after applying significance filter"
                        % expt_id
                    )

            if len(accepted_refls) == 0:
                raise Sorry("No reflections left after applying significance filter")
            experiments = accepted_expts
            reflections = accepted_refls

        # Delete the shoeboxes used for intermediate calculations, if requested
        if params.integration.debug.delete_shoeboxes and "shoebox" in reflections:
            del reflections["shoebox"]

        # Save the reflections
        self.save_reflections(reflections, params.output.reflections)
        self.save_experiments(experiments, params.output.experiments)

        # Write a report if requested
        if params.output.report is not None:
            integrator.report().as_file(params.output.report)

        # Print the total time taken
        logger.info("\nTotal time taken: %f" % (time() - start_time))

        return experiments, reflections

    def process_reference(self, reference):
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
        mask = reference.get_flags(reference.flags.bad_for_refinement, all=False)
        if mask.count(True) > 0:
            rubbish.extend(reference.select(mask))
            reference.del_selected(mask)
            logger.info(
                " removing %d reflections marked as bad for refinement"
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
        if (reference["panel"] == reference["shoebox"].panels()).count(False) > 0:
            raise RuntimeError(
                'reflection table "panel" column does not match "shoebox" panel'
            )
        logger.info(" using %d indexed reflections" % len(reference))
        logger.info(" found %d junk reflections" % len(rubbish))
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

    def save_reflections(self, reflections, filename):
        """ Save the reflections to file. """
        from time import time

        st = time()
        logger.info("Saving %d reflections to %s" % (len(reflections), filename))
        reflections.as_pickle(filename)
        logger.info(" time taken: %g" % (time() - st))

    def save_experiments(self, experiments, filename):
        """ Save the profile model parameters. """
        from time import time
        from dxtbx.model.experiment_list import ExperimentListDumper

        st = time()
        logger.info("Saving the experiments to %s" % filename)
        dump = ExperimentListDumper(experiments)
        with open(filename, "w") as outfile:
            outfile.write(dump.as_json())
        logger.info(" time taken: %g" % (time() - st))

    def sample_predictions(self, experiments, predicted, params):
        """ Select a random sample of the predicted reflections to integrate. """

        nref_per_degree = params.sampling.reflections_per_degree
        min_sample_size = params.sampling.minimum_sample_size
        max_sample_size = params.sampling.maximum_sample_size

        # this code is very similar to David's code in algorithms/refinement/reflection_manager.py!

        # constants
        from math import pi

        RAD2DEG = 180.0 / pi

        working_isel = flex.size_t()
        for iexp, exp in enumerate(experiments):

            sel = predicted["id"] == iexp
            isel = sel.iselection()
            # refs = self._reflections.select(sel)
            nrefs = sample_size = len(isel)

            # set sample size according to nref_per_degree (per experiment)
            if exp.scan and nref_per_degree:
                sweep_range_rad = exp.scan.get_oscillation_range(deg=False)
                width = abs(sweep_range_rad[1] - sweep_range_rad[0]) * RAD2DEG
                sample_size = int(nref_per_degree * width)
            else:
                sweep_range_rad = None

            # adjust sample size if below the chosen limit
            sample_size = max(sample_size, min_sample_size)

            # set maximum sample size if requested
            if max_sample_size:
                sample_size = min(sample_size, max_sample_size)

            # determine subset and collect indices
            if sample_size < nrefs:
                isel = isel.select(flex.random_selection(nrefs, sample_size))
            working_isel.extend(isel)

        # create subset
        return predicted.select(working_isel)

    def exclude_images(self, experiments, exclude_images):

        if exclude_images is not None and len(exclude_images) > 0:
            for experiment in experiments:
                imageset = experiment.imageset
                for index in exclude_images:
                    imageset.mark_for_rejection(index, True)

        return experiments

    def split_for_scan_range(self, experiments, reference, scan_range):
        """ Update experiments when scan range is set. """
        from dxtbx.model.experiment_list import ExperimentList
        from dxtbx.model.experiment_list import Experiment

        # Only do anything is the scan range is set
        if scan_range is not None and len(scan_range) > 0:

            # Ensure that all experiments have the same imageset and scan
            iset = [e.imageset for e in experiments]
            scan = [e.scan for e in experiments]
            assert all(x == iset[0] for x in iset)
            assert all(x == scan[0] for x in scan)

            # Get the imageset and scan
            iset = experiments[0].imageset
            scan = experiments[0].scan

            # Get the array range
            if scan is not None:
                frames_start, frames_end = scan.get_array_range()
                assert scan.get_num_images() == len(iset)
            else:
                frames_start, frames_end = (0, len(iset))

            # Create the new lists
            new_experiments = ExperimentList()
            new_reference_all = reference.split_by_experiment_id()
            new_reference = flex.reflection_table()
            for i in range(len(new_reference_all) - len(experiments)):
                new_reference_all.append(flex.reflection_table())
            assert len(new_reference_all) == len(experiments)

            # Loop through all the scan ranges and create a new experiment list with
            # the requested scan ranges.
            for scan_start, scan_end in scan_range:
                # Validate the requested scan range
                if scan_end == scan_start:
                    raise Sorry(
                        "Scan range end must be higher than start; pass {},{} for single image".format(
                            scan_start, scan_start + 1
                        )
                    )
                if scan_end < scan_start:
                    raise Sorry("Scan range must be in ascending order")
                elif scan_start < frames_start or scan_end > frames_end:
                    raise Sorry(
                        "Scan range must be within image range {}..{}".format(
                            frames_start, frames_end
                        )
                    )

                assert scan_end > scan_start
                assert scan_start >= frames_start
                assert scan_end <= frames_end

                index_start = scan_start - frames_start
                index_end = index_start + (scan_end - scan_start)
                assert index_start < index_end
                assert index_start >= 0
                assert index_end <= len(iset)
                new_iset = iset[index_start:index_end]
                if scan is None:
                    new_scan = None
                else:
                    new_scan = scan[index_start:index_end]
                for i, e1 in enumerate(experiments):
                    e2 = Experiment()
                    e2.beam = e1.beam
                    e2.detector = e1.detector
                    e2.goniometer = e1.goniometer
                    e2.crystal = e1.crystal
                    e2.profile = e1.profile
                    e2.imageset = new_iset
                    e2.scan = new_scan
                    new_reference_all[i]["id"] = flex.int(
                        len(new_reference_all[i]), len(new_experiments)
                    )
                    new_reference.extend(new_reference_all[i])
                    new_experiments.append(e2)
            experiments = new_experiments
            reference = new_reference

            # Print some information
            logger.info(
                "Modified experiment list to integrate over requested scan range"
            )
            for scan_start, scan_end in scan_range:
                logger.info(" scan_range = %d -> %d" % (scan_start, scan_end))
            logger.info("")

        # Return the experiments
        return experiments, reference


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
