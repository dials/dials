#!/usr/bin/env python
#
# dials.refine.py
#
#  Copyright (C) 2013 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: James Parkhurst and David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# DIALS_ENABLE_COMMAND_LINE_COMPLETION

"""
Refine the diffraction geometry of input experiments against the input indexed
reflections. For rotation scans, the model may be either static (the same for
all reflections) or scan-varying (dependent on image number in the scan).
Other basic parameters include control over output filenames, fixing of
certain parameters of each model and options that control the number of
reflections used in refinement.

Examples::

  dials.refine indexed.expt indexed.refl

  dials.refine indexed.expt indexed.refl scan_varying=(False/True/Auto)

"""

from __future__ import absolute_import, division, print_function
import sys
import logging
from time import time
import dials.util
import libtbx.phil
from libtbx import Auto
from dials.array_family import flex
from dials.algorithms.refinement import RefinerFactory
from dials.algorithms.refinement import DialsRefineConfigError, DialsRefineRuntimeError
from dials.algorithms.refinement.corrgram import create_correlation_plots
from dials.util import Sorry

logger = logging.getLogger("dials.command_line.refine")

# The phil scope
phil_scope = libtbx.phil.parse(
    """

  output {
    experiments = refined.expt
      .type = str
      .help = "The filename for refined experimental models"

    reflections = refined.refl
      .type = str
      .help = "The filename for reflections with updated predictions"

    include_unused_reflections = True
      .type = bool
      .help = "If True, keep reflections unused in refinement in updated"
              "reflections file. Otherwise, remove them"
      .expert_level = 1

    matches = None
      .type = str
      .help = "The filename for output of the reflection table for reflections"
              "used in refinement, containing extra columns used internally."
              "Intended for debugging purposes only"
      .expert_level = 2

    centroids = None
      .type = str
      .help = "The filename for the table of centroids at the end of"
              "refinement"
      .expert_level = 1

    parameter_table = None
      .type = str
      .help = "The filename for the table of scan varying parameter values"
      .expert_level = 1

    log = dials.refine.log
      .type = str

    include scope dials.algorithms.refinement.corrgram.phil_scope

    history = None
      .type = str
      .help = "The filename for output of the refinement history json"
      .expert_level = 1
  }

  n_static_macrocycles = 1
    .type = int(value_min=1)
    .help = "Number of macro-cycles of static refinement to perform"

  include scope dials.algorithms.refinement.refiner.phil_scope
""",
    process_includes=True,
)

# local overrides for refiner.phil_scope
phil_overrides = libtbx.phil.parse(
    """
  refinement
  {
    parameterisation.scan_varying = Auto
  }
"""
)

working_phil = phil_scope.fetch(sources=[phil_overrides])


def write_centroids_table(refiner, filename):
    """Write a table of centroids from refinement for debugging.

    Args:
        refiner: A dials.algorithms.refinement.refiner.Refiner.
        filename (str): The name of the file to which to write.

    Returns:
        None
    """

    matches = refiner.get_matches()

    header = "H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t" "Y_calc\tPhi_calc"
    msg_temp = "%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%5.3f\t%9.6f"
    has_del_psi = "delpsical.rad" in matches
    if has_del_psi:
        header += "\tDelta_Psi"
        msg_temp += "\t%9.6f"
    header += "\n"
    msg_temp += "\n"

    with open(filename, "w") as f:
        f.write(header)

        for m in matches:
            (h, k, l) = m["miller_index"]
            frame = m["xyzobs.px.value"][2]
            x_obs, y_obs, phi_obs = m["xyzobs.mm.value"]
            x_calc, y_calc, phi_calc = m["xyzcal.mm"]
            if has_del_psi:
                del_psi = m["delpsical.rad"]
                msg = msg_temp % (
                    h,
                    k,
                    l,
                    frame,
                    x_obs,
                    y_obs,
                    phi_obs,
                    x_calc,
                    y_calc,
                    phi_calc,
                    del_psi,
                )
            else:
                msg = msg_temp % (
                    h,
                    k,
                    l,
                    frame,
                    x_obs,
                    y_obs,
                    phi_obs,
                    x_calc,
                    y_calc,
                    phi_calc,
                )
            f.write(msg)


def run_macrocycle(params, reflections, experiments):
    """Run one macrocycle of refinement.

    One macrocycle of refinement is run, as specified by the PHIL
    parameters, using the centroids from the supplied reflections
    and the initial experimental geometry taken from experiments.


    Args:
        params: The working PHIL parameters.
        reflections: A reflection table containing observed centroids
        experiments: The initial dxtbx experimental geometry models

    Returns:
        tuple: The Refiner, the reflection table with updated predictions
            and flags, and the refinement history object.
    """
    # Get the refiner
    logger.info("Configuring refiner")
    try:
        refiner = RefinerFactory.from_parameters_data_experiments(
            params, reflections, experiments
        )
    except DialsRefineConfigError as e:
        raise Sorry(e.message)

    # Refine the geometry
    nexp = len(experiments)
    if nexp == 1:
        logger.info("Performing refinement of a single Experiment...")
    else:
        logger.info("Performing refinement of {} Experiments...".format(nexp))

    # Refine and get the refinement history
    try:
        history = refiner.run()
    except DialsRefineRuntimeError as e:
        raise Sorry(e.message)

    # Update predictions for all indexed reflections
    logger.info("Updating predictions for indexed reflections")
    preds = refiner.predict_for_indexed()

    # just copy over the columns of interest or columns that may have been
    # updated, leaving behind things added by e.g. scan-varying refinement
    # such as 'block', the U, B and UB matrices and gradients.
    for key in preds:
        if key in reflections.keys() or key in [
            "s1",
            "xyzcal.mm",
            "xyzcal.px",
            "entering",
            "delpsical.rad",
        ]:
            reflections[key] = preds[key]

    # set refinement flags
    assert len(preds) == len(reflections)
    reflections.unset_flags(
        flex.size_t_range(len(reflections)),
        reflections.flags.excluded_for_refinement
        | reflections.flags.used_in_refinement
        | reflections.flags.centroid_outlier
        | reflections.flags.predicted,
    )
    reflections.set_flags(
        preds.get_flags(preds.flags.excluded_for_refinement),
        reflections.flags.excluded_for_refinement,
    )
    reflections.set_flags(
        preds.get_flags(preds.flags.centroid_outlier),
        reflections.flags.centroid_outlier,
    )
    reflections.set_flags(
        preds.get_flags(preds.flags.used_in_refinement),
        reflections.flags.used_in_refinement,
    )
    reflections.set_flags(
        preds.get_flags(preds.flags.predicted), reflections.flags.predicted
    )

    return refiner, reflections, history


def run_dials_refine(experiments, reflections, params):
    """Functional interface to tasks performed by the program dials.refine.

    This runs refinement according to the PHIL parameters in params, using
    taking the initial experimental geometry from experiments, and the
    observed reflection centroids from reflections.

    Static refinement can be done in a number of macrocycles. By default,
    for scans, this will be followed by a macrocycle of scan-varying
    refinement.


    Args:
        experiments: The initial dxtbx experimental geometry models
        reflections: A reflection table containing observed centroids
        params: The working PHIL parameters

    Returns:
        tuple: The refined experiments, the updated reflection table, the
            Refiner object and the refinement history object.
    """

    # Modify options if necessary
    if params.output.correlation_plot.filename is not None:
        params.refinement.refinery.journal.track_parameter_correlation = True

    # scan_varying=Auto is special to dials.refine. Refiner expects it to be
    # True or False, so catch and set here
    scan_varying = params.refinement.parameterisation.scan_varying
    if scan_varying is Auto:
        params.refinement.parameterisation.scan_varying = False

    if params.n_static_macrocycles == 1:
        refiner, reflections, history = run_macrocycle(params, reflections, experiments)
        experiments = refiner.get_experiments()
    else:
        for i in range(params.n_static_macrocycles):
            logger.info("\nStatic refinement macrocycle {}".format(i + 1))
            refiner, reflections, history = run_macrocycle(
                params, reflections, experiments
            )
            experiments = refiner.get_experiments()

    # Scan-varying macrocycle, if appropriate
    if scan_varying is Auto and refiner.experiment_type == "scans":
        logger.info("\nScan-varying refinement")
        params.refinement.parameterisation.scan_varying = True
        refiner, reflections, history = run_macrocycle(params, reflections, experiments)
        experiments = refiner.get_experiments()

    return experiments, reflections, refiner, history


def run(args=None, phil=working_phil):
    """
    Set up refinement from command line options, files and PHIL parameters.
    Run refinement and save output files as specified.

    Called when running dials.refine as a command-line program

    Args:
        args (list): Additional command-line arguments
        phil: The working PHIL parameters

    Returns:
        None
    """

    import dials.util.log
    from dials.util.options import OptionParser
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_experiments
    import libtbx.load_env

    start_time = time()

    # The script usage
    usage = (
        "usage: %s [options] [param.phil] "
        "models.expt observations.refl" % libtbx.env.dispatcher_name
    )

    # Create the parser
    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    # Parse the command line
    params, options = parser.parse_args(args=args, show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    # Configure the logging
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    # Try to load the models and data
    nexp = len(experiments)
    if nexp == 0:
        sys.exit("No Experiments found in the input")
    if len(reflections) == 0:
        sys.exit("No reflection data found in the input")
    if len(reflections) > 1:
        sys.exit("Only one reflections list can be imported at present")
    reflections = reflections[0]

    # check input is suitable
    msg = (
        "The supplied reflection table does not have the required data " + "column: {0}"
    )
    for key in ["xyzobs.mm.value", "xyzobs.mm.variance"]:
        if key not in reflections:
            msg = msg.format(key)
            raise dials.util.Sorry(msg)

    from dials.util.version import dials_version

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Warn about potentially unhelpful options
    if params.refinement.mp.nproc > 1:
        logger.warning(
            "WARNING: setting nproc > 1 is only helpful in rare "
            "circumstances. It is not recommended for typical data processing "
            "tasks.\n"
        )

    # Run refinement
    experiments, reflections, refiner, history = run_dials_refine(
        experiments, reflections, params
    )

    # For the usual case of refinement of one crystal, print that model for information
    crystals = experiments.crystals()
    if len(crystals) == 1:
        logger.info("")
        logger.info("Final refined crystal model:")
        logger.info(crystals[0])

    # Write table of centroids to file, if requested
    if params.output.centroids:
        logger.info(
            "Writing table of centroids to '{}'".format(params.output.centroids)
        )
        write_centroids_table(refiner, params.output.centroids)

    # Write scan-varying parameters to file, if there were any
    if params.output.parameter_table:
        scans = experiments.scans()
        if len(scans) > 1:
            logger.info(
                "Writing a scan-varying parameter table is only supported "
                "for refinement of a single scan"
            )
        else:
            scan = scans[0]
            text = refiner.get_param_reporter().varying_params_vs_image_number(
                scan.get_array_range()
            )
            if text:
                logger.info(
                    "Writing scan-varying parameter table to {}".format(
                        params.output.parameter_table
                    )
                )
                f = open(params.output.parameter_table, "w")
                f.write(text)
                f.close()
            else:
                logger.info("No scan-varying parameter table to write")

    # Save the refined experiments to file
    output_experiments_filename = params.output.experiments
    logger.info("Saving refined experiments to {}".format(output_experiments_filename))
    from dxtbx.model.experiment_list import ExperimentListDumper

    dump = ExperimentListDumper(experiments)
    dump.as_json(output_experiments_filename)

    # Save reflections with updated predictions if requested (allow to switch
    # this off if it is a time-consuming step)
    if params.output.reflections:
        logger.info(
            "Saving reflections with updated predictions to {}".format(
                params.output.reflections
            )
        )
        if params.output.include_unused_reflections:
            reflections.as_pickle(params.output.reflections)
        else:
            sel = reflections.get_flags(reflections.flags.used_in_refinement)
            reflections.select(sel).as_pickle(params.output.reflections)

    # Save matches to file for debugging
    if params.output.matches:
        matches = refiner.get_matches()
        logger.info(
            "Saving matches (use for debugging purposes) to {}".format(
                params.output.matches
            )
        )
        matches.as_pickle(params.output.matches)

    # Create correlation plots
    if params.output.correlation_plot.filename is not None:
        create_correlation_plots(refiner, params.output)

    # Save refinement history
    if params.output.history:
        logger.info(
            "Saving refinement step history to {}".format(params.output.history)
        )
        history.to_json_file(params.output.history)

    # Log the total time taken
    logger.info("\nTotal time taken: {:.2f}s".format(time() - start_time))


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()
