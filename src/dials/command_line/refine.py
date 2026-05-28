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

from __future__ import annotations

import copy
import logging
import sys
from collections import namedtuple

import libtbx.phil
from dxtbx.model.experiment_list import ExperimentList
from libtbx import Auto

import dials.util
import dials.util.log
from dials.algorithms.refinement import (
    DialsRefineConfigError,
    DialsRefineRuntimeError,
    RefinerFactory,
)
from dials.algorithms.refinement.corrgram import create_correlation_plots
from dials.array_family import flex
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

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

  separate_independent_sets = True
    .type = bool
    .help = "If true, the experiment list will be separated into independent groups"
            "that do not share models, and these groups will be refined separately."

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

    reflections.outlier.nproc = Auto
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

    header = "H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\tY_calc\tPhi_calc"
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
    refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments
    )

    # Refine the geometry
    nexp = len(experiments)
    if nexp == 1:
        logger.info("Performing refinement of a single Experiment...")
    else:
        logger.info(f"Performing refinement of {nexp} Experiments...")

    # Refine and get the refinement history
    history = refiner.run()

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


def _find_disjoint_sets(experiments):
    # Extract parameterisable models from the experiments
    models = []
    for experiment in experiments:
        models.append(
            [
                m
                for m in [
                    experiment.beam,
                    experiment.crystal,
                    experiment.detector,
                    experiment.goniometer,
                ]
                if m is not None
            ]
        )

    # Record first set of models
    sets = [set(models[0])]
    ids = [
        [0],
    ]

    # Go through all other models, matching to previous sets
    for i, m in enumerate(models[1:]):
        new_set = set(m)
        disj = [new_set.isdisjoint(s) for s in sets]
        if all(disj):
            # no shared models, so form a new set
            sets.append(new_set)
            ids.append([i + 1])
        else:
            # models shared with at least one existing set
            for j, d in enumerate(disj):
                if d:
                    continue
                sets[j].update(new_set)
                ids[j].append(i + 1)

    # Now combine lists in ids if any are not unique (https://stackoverflow.com/a/4842897)
    accepted = []
    while len(ids) > 0:
        first, *rest = ids
        first = set(first)

        lf = -1
        while len(first) > lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r))) > 0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        accepted.append(first)
        ids = rest
    accepted = [sorted(s) for s in accepted]

    return accepted


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

    # Warn about potentially unhelpful options
    if params.refinement.mp.nproc > 1:
        logger.warning(
            "Setting nproc > 1 is only helpful in rare "
            "circumstances. It is not recommended for typical data processing "
            "tasks."
        )

    if params.refinement.parameterisation.scan_varying is not False:
        # duplicate crystal if necessary for scan varying - will need
        # to compare the scans with crystals - if not 1:1 will need to
        # split the crystals

        crystal_has_scan = {}
        for j, e in enumerate(experiments):
            if e.crystal in crystal_has_scan:
                if e.scan is not crystal_has_scan[e.crystal]:
                    logger.info(
                        "Duplicating crystal model for scan-varying refinement of experiment %d",
                        j,
                    )
                    e.crystal = copy.deepcopy(e.crystal)
            else:
                crystal_has_scan[e.crystal] = e.scan

    # Modify options if necessary
    if params.output.correlation_plot.filename is not None:
        params.refinement.refinery.journal.track_parameter_correlation = True

    # scan_varying=Auto is special to dials.refine. Refiner expects it to be
    # True or False, so catch and set here
    scan_varying = params.refinement.parameterisation.scan_varying
    if scan_varying is Auto:
        params.refinement.parameterisation.scan_varying = False

    # Similarly, keep track of sparse to reset that for scan-varying macrocycle
    sparse = params.refinement.parameterisation.sparse

    # Look for disjoint sets of experiments
    disjoint_sets = _find_disjoint_sets(experiments)

    # Look for "hidden" links between experiments from constraints or restraints
    bp = params.refinement.parameterisation.beam
    ucp = params.refinement.parameterisation.crystal.unit_cell
    op = params.refinement.parameterisation.crystal.orientation
    dp = params.refinement.parameterisation.detector
    gp = params.refinement.parameterisation.goniometer
    crosslinks = any(
        (
            bp.constraints,
            ucp.constraints,
            ucp.restraints.tie_to_group,
            op.constraints,
            dp.constraints,
            gp.constraints,
        )
    )

    RefinementSet = namedtuple(
        "RefinementSet", ["experiments", "reflections", "params", "original_ids"]
    )
    if len(disjoint_sets) == 1 or not params.separate_independent_sets:
        # No splitting required, this is one interdependent refinement job
        refinement_sets = [
            RefinementSet(experiments, reflections, params, disjoint_sets[0])
        ]
    elif crosslinks:
        # Refuse to split because there are restraints or constraints present
        logger.warning(
            "The experiments contain disjoint subsets that do not share models. "
            "However, these will not be refined independently, because restraints "
            "or constraints may link models between the subsets."
        )
        refinement_sets = [
            RefinementSet(experiments, reflections, params, disjoint_sets[0])
        ]
    else:
        # If outlier rejection is meant to be done across all experiments then
        # do that once here
        if not params.refinement.reflections.outlier.separate_experiments:
            reflections = RefinerFactory.reflections_after_outlier_rejection(
                params, reflections, experiments
            )
            params.refinement.reflections.outlier.algorithm = "null"

        # Set large objects in params to None for copying
        params.input.reflections = None
        params.input.experiments = None

        # Split into independent refinement jobs
        refinement_sets = []
        for ids in disjoint_sets:
            el = ExperimentList()
            for i in ids:
                el.append(experiments[i])

            refl = flex.reflection_table()
            new_id = 0
            for i in ids:
                refl_one_experiment = reflections.select(reflections["id"] == i)
                refl_one_experiment["id"] = flex.int(len(refl_one_experiment), new_id)
                new_id += 1
                refl.extend(refl_one_experiment)
            refinement_sets.append(RefinementSet(el, refl, copy.deepcopy(params), ids))

        # Report on independent refinement sets
        logger.info(
            "The experiments have been separated into independent groups that "
            "do not share models.\nRefinement will occur separately for each group:"
        )
        header = ["Group", "Experiment ids"]
        rows = []
        for i, ids in enumerate(disjoint_sets):
            rows.append([str(i), " ".join(str(e) for e in ids)])
        logger.info(dials.util.tabulate(rows, header))

    refinement_results = []
    for rs in refinement_sets:
        experiments = rs.experiments
        reflections = rs.reflections
        params = rs.params
        if len(refinement_sets) > 1:
            logger.info(
                "\nSelected group of experiments to refine with original ids: "
                + " ".join([str(i) for i in rs.original_ids])
            )
        if params.n_static_macrocycles == 1:
            refiner, reflections, history = run_macrocycle(
                params, reflections, experiments
            )
            experiments = refiner.get_experiments()
        else:
            for i in range(params.n_static_macrocycles):
                logger.info("\nStatic refinement macrocycle %s", i + 1)
                refiner, reflections, history = run_macrocycle(
                    params, reflections, experiments
                )
                experiments = refiner.get_experiments()

        # Scan-varying macrocycle, if appropriate
        if scan_varying is Auto and refiner.experiment_type == "scans":
            logger.info("\nScan-varying refinement")
            params.refinement.parameterisation.scan_varying = True
            params.refinement.parameterisation.sparse = sparse
            refiner, reflections, history = run_macrocycle(
                params, reflections, experiments
            )
            experiments = refiner.get_experiments()

        refinement_results.append((experiments, reflections, refiner, history))

    if len(refinement_results) == 1:
        experiments, reflections, refiner, history = refinement_results[0]
    else:
        # Rejoin results in the expected order of experiments
        experiments = {}
        reflections = flex.reflection_table()
        table_headers = []
        table_rows = {}
        for (el, refl, refiner, _), ids in zip(refinement_results, disjoint_sets):
            header, rows = refiner.calc_exp_rmsd_table()
            id_col = flex.int(len(refl))
            for new_id, orig_id in enumerate(ids):
                experiments[orig_id] = el[new_id]
                id_col.set_selected(refl["id"] == new_id, orig_id)
                rows[new_id][0] = str(orig_id)
                table_rows[orig_id] = rows[new_id]
            refl["id"] = id_col
            reflections.extend(refl)
            table_headers.append(header)

        experiments = ExperimentList([experiments[i] for i in range(len(experiments))])

        logger.info(
            "\nIndependently-refined groups of experiments have been recombined"
        )

        # There are multiple refiners and history objects. We don't have a way
        # to combine these usefully, so return None to avoid misleading the
        # caller that these are relevant for the full refinement
        refiner = None
        history = None

        # Report on RMSD by experiments after joining
        rows = [table_rows[i] for i in range(len(table_rows))]
        header = table_headers[0]
        if not (h == header for h in table_headers[1:]):
            logger.warning(
                "Cannot calculate RMSDs by experiment as units are inconsistent"
            )
        else:
            logger.info("\nRMSDs by experiment:")
            logger.info(dials.util.tabulate(rows, header))

    return experiments, reflections, refiner, history


@dials.util.show_mail_handle_errors()
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

    # The script usage
    usage = "usage: dials.refine [options] [param.phil] models.expt observations.refl"

    # Create the parser
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    # Parse the command line
    params, options = parser.parse_args(args=args, show_diff_phil=False)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Configure the logging
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    # Try to load the models and data
    nexp = len(experiments)
    if nexp == 0 or len(reflections) == 0:
        parser.print_help()
        return
    if len(reflections) > 1:
        sys.exit("Only one reflections list can be imported at present")
    reflections = reflections[0]

    # check input is suitable
    msg = (
        "The supplied reflection table does not have the required data " + "column: {0}"
    )
    for key in ["xyzobs.mm.value", "xyzobs.mm.variance"]:
        if key not in reflections:
            sys.exit(msg.format(key))

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Run refinement
    try:
        experiments, reflections, refiner, history = run_dials_refine(
            experiments, reflections, params
        )
    except (DialsRefineConfigError, DialsRefineRuntimeError) as e:
        sys.exit(str(e))

    # For the usual case of refinement of one crystal, print that model for information
    crystals = experiments.crystals()
    if len(crystals) == 1:
        logger.info("")
        logger.info("Final refined crystal model:")
        logger.info(crystals[0])

    # Write table of centroids to file, if requested
    if params.output.centroids:
        if not refiner:
            logger.warning(
                "Cannot write table of centroids as a single refiner object is not available"
            )
        else:
            logger.info(f"Writing table of centroids to '{params.output.centroids}'")
            write_centroids_table(refiner, params.output.centroids)

    # Write scan-varying parameters to file, if there were any
    if params.output.parameter_table:
        scans = experiments.scans()
        if len(scans) > 1:
            logger.info(
                "Writing a scan-varying parameter table is only supported "
                "for refinement of a single scan"
            )
        elif not refiner:
            logger.warning(
                "Cannot write scan-varying parameter table as a single refiner object is not available"
            )
        else:
            scan = scans[0]
            text = refiner.get_param_reporter().varying_params_vs_image_number(
                scan.get_array_range()
            )
            if text:
                logger.info(
                    "Writing scan-varying parameter table to %s",
                    params.output.parameter_table,
                )
                f = open(params.output.parameter_table, "w")
                f.write(text)
                f.close()
            else:
                logger.info("No scan-varying parameter table to write")

    # Save matches to file for debugging
    if params.output.matches:
        if not refiner:
            logger.warning(
                "Cannot write matches to file as a single refiner object is not available"
            )
        else:
            matches = refiner.get_matches()
            logger.info(
                "Saving matches (use for debugging purposes) to %s",
                params.output.matches,
            )
            matches.as_file(params.output.matches)

    # Create correlation plots
    if params.output.correlation_plot.filename is not None:
        if not refiner:
            logger.warning(
                "Cannot create correlation plots as a single refiner object is not available"
            )
        else:
            create_correlation_plots(refiner, params.output)

    # Save refinement history
    if params.output.history:
        if not history:
            logger.warning(
                "Cannot write refinement step history as a single history object is not available"
            )
        else:
            logger.info(f"Saving refinement step history to {params.output.history}")
            history.to_json_file(params.output.history)

    # Save the refined experiments to file
    output_experiments_filename = params.output.experiments
    logger.info(f"Saving refined experiments to {output_experiments_filename}")
    experiments.as_file(output_experiments_filename)

    # Save reflections with updated predictions. This causes a spike in memory
    # usage (https://github.com/dials/dials/issues/2024), so delete big objects
    # we no longer need first
    if params.output.reflections:
        del experiments
        del refiner
        del history
        logger.info(
            "Saving reflections with updated predictions to %s",
            params.output.reflections,
        )
        if params.output.include_unused_reflections:
            reflections.as_file(params.output.reflections)
        else:
            sel = reflections.get_flags(reflections.flags.used_in_refinement)
            reflections.select(sel).as_file(params.output.reflections)


if __name__ == "__main__":
    run()
