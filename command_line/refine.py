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
from __future__ import absolute_import, division, print_function

import logging
from libtbx.utils import Sorry
from libtbx import Auto
from dials.array_family import flex
from dials.algorithms.refinement import RefinerFactory

logger = logging.getLogger('dials.command_line.refine')

help_message = '''

Refine the diffraction geometry of input experiments against the input indexed
reflections. For rotation scans, the model may be either static (the same for
all reflections) or scan-varying (dependent on image number in the scan).
Other basic parameters include control over output filenames, fixing of
certain parameters of each model and options that control the number of
reflections used in refinement.

Examples::

  dials.refine experiments.json indexed.pickle

  dials.refine experiments.json indexed.pickle scan_varying=True

'''

# The phil scope
from libtbx.phil import parse
phil_scope = parse('''

  output {
    experiments = refined_experiments.json
      .type = str
      .help = "The filename for refined experimental models"

    reflections = refined.pickle
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

    debug_log = dials.refine.debug.log
      .type = str

    correlation_plot
      .expert_level = 1
    {
      filename = None
        .type = str
        .help = "The base filename for output of plots of parameter"
                "correlations. A file extension may be added to control"
                "the type of output file, if it is one of matplotlib's"
                "supported types. A pickle file with the same base filename"
                "will also be created, containing the correlation matrix and"
                "column labels for later inspection, replotting etc."

      col_select = None
        .type = strings
        .help = "Specific columns to include in the plots of parameter"
                "correlations, either specifed by parameter name or 0-based"
                "column index. Defaults to all columns."
                "This option is useful when there is a large number of"
                "parameters"

      steps = None
        .type = ints(value_min=0)
        .help = "Steps for which to make correlation plots. By default only"
                "the final step is plotted. Uses zero-based numbering, so"
                "the first step is numbered 0."
    }

    history = None
      .type = str
      .help = "The filename for output of the refinement history pickle"
      .expert_level = 1
  }

  do_scan_varying_macrocycle = False
    .type = bool
    .help = "Override whatever scan_varying is set to and force one round of"
            "static refinement followed by one round of scan-varying refinement"
    .expert_level = 3

  include scope dials.algorithms.refinement.refiner.phil_scope
''', process_includes=True)

# local overrides for refiner.phil_scope
phil_overrides = parse('''
  refinement
  {
    verbosity = 2
  }
''')

working_phil = phil_scope.fetch(sources=[phil_overrides])

class Script(object):
  '''A class for running the script.'''

  def __init__(self, phil=working_phil):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env


    # The script usage
    usage  = "usage: %s [options] [param.phil] " \
             "experiments.json reflections.pickle" \
               % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil,
      read_reflections=True,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

  def write_centroids_table(self, refiner, filename):

    matches = refiner.get_matches()

    header = ("H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t"
        "Y_calc\tPhi_calc")
    msg_temp = ("%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%5.3f\t%9.6f")
    has_del_psi = 'delpsical.rad' in matches
    if has_del_psi:
      header += "\tDelta_Psi"
      msg_temp += "\t%9.6f"
    header += "\n"
    msg_temp += "\n"

    with open(filename,"w") as f:
      f.write(header)

      for m in matches:
        (h, k, l) = m['miller_index']
        frame = m['xyzobs.px.value'][2]
        x_obs, y_obs, phi_obs = m['xyzobs.mm.value']
        x_calc, y_calc, phi_calc = m['xyzcal.mm']
        if has_del_psi:
          del_psi = m['delpsical.rad']
          msg = msg_temp % (h, k, l,
                          frame, x_obs, y_obs, phi_obs,
                          x_calc, y_calc, phi_calc, del_psi)
        else:
          msg = msg_temp % (h, k, l,
                          frame, x_obs, y_obs, phi_obs,
                          x_calc, y_calc, phi_calc)
        f.write(msg)

  @staticmethod
  def check_input(reflections):
    '''Check the input is suitable for refinement. So far just check keys in
    the reflection table. Maybe later check experiments have overlapping models
    etc.'''

    msg = "The supplied reflection table does not have the required data " + \
      "column: {0}"
    for key in ["xyzobs.mm.value", "xyzobs.mm.variance"]:
      if key not in reflections:
        msg = msg.format(key)
        raise Sorry(msg)

    # FIXME add other things to be checked here
    return

  @staticmethod
  def run_macrocycle(params, reflections, experiments):

    # Get the refiner
    logger.info('Configuring refiner')
    refiner = RefinerFactory.from_parameters_data_experiments(params,
        reflections, experiments)

    # Refine the geometry
    nexp = len(experiments)
    if nexp == 1:
      logger.info('Performing refinement of a single Experiment...')
    else:
      logger.info('Performing refinement of {0} Experiments...'.format(nexp))

    # Refine and get the refinement history
    history = refiner.run()
    return refiner, history

  def run(self, args=None):
    '''Execute the script.'''
    from time import time
    import six.moves.cPickle as pickle
    from dials.util import log
    from dials.util.options import flatten_reflections, flatten_experiments

    start_time = time()

    # Parse the command line
    params, options = self.parser.parse_args(args=args, show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    # Try to load the models and data
    nexp = len(experiments)
    if nexp == 0:
      print("No Experiments found in the input")
      self.parser.print_help()
      return
    if len(reflections) == 0:
      print("No reflection data found in the input")
      self.parser.print_help()
      return
    if len(reflections) > 1:
      raise Sorry("Only one reflections list can be imported at present")
    reflections = reflections[0]

    self.check_input(reflections)

    if __name__ == '__main__':
      # Configure the logging
      log.config(info=params.output.log,
        debug=params.output.debug_log)

    from dials.util.version import dials_version
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    # Warn about potentially unhelpful options
    if params.refinement.mp.nproc > 1:
      logger.warning("WARNING: setting nproc > 1 is only helpful in rare "
        "circumstances. It is not recommended for typical data processing "
        "tasks.\n")

    # Modify options if necessary
    if params.output.correlation_plot.filename is not None:
      params.refinement.refinery.journal.track_parameter_correlation = True
    do_sv_macrocycle = False
    if params.do_scan_varying_macrocycle:
      params.refinement.parameterisation.scan_varying = False
      do_sv_macrocycle = True

    # Do refinement
    if do_sv_macrocycle:
      logger.info("\nScan-static refinement")
      refiner, history = self.run_macrocycle(params, reflections, experiments)
      logger.info("\nScan-varying refinement")
      params.refinement.parameterisation.scan_varying = True
      refiner, history = self.run_macrocycle(params, reflections, experiments)
    else:
      refiner, history = self.run_macrocycle(params, reflections, experiments)

    if params.output.centroids:
      logger.info("Writing table of centroids to '{0}'".format(
        params.output.centroids))
      self.write_centroids_table(refiner, params.output.centroids)

    # Get the refined experiments
    experiments = refiner.get_experiments()

    # Write scan-varying parameters to file, if there were any
    if params.output.parameter_table:
      scans = experiments.scans()
      if len(scans) > 1:
        logger.info("Writing a scan-varying parameter table is only supported "
             "for refinement of a single scan")
      else:
        scan = scans[0]
        text = refiner.get_param_reporter().varying_params_vs_image_number(
            scan.get_array_range())
        if text:
          logger.info("Writing scan-varying parameter table to {0}".format(
            params.output.parameter_table))
          f = open(params.output.parameter_table,"w")
          f.write(text)
          f.close()
        else:
          logger.info("No scan-varying parameter table to write")

    crystals = experiments.crystals()
    if len(crystals) == 1:
      # output the refined model for information
      logger.info('')
      logger.info('Final refined crystal model:')
      logger.info(crystals[0])

    # Save the refined experiments to file
    output_experiments_filename = params.output.experiments
    logger.info('Saving refined experiments to {0}'.format(output_experiments_filename))
    from dxtbx.model.experiment_list import ExperimentListDumper
    dump = ExperimentListDumper(experiments)
    dump.as_json(output_experiments_filename)

    # Save reflections with updated predictions if requested (allow to switch
    # this off if it is a time-consuming step)
    if params.output.reflections:
      # Update predictions for all indexed reflections
      logger.info('Updating predictions for indexed reflections')
      preds = refiner.predict_for_indexed()

      # just copy over the columns of interest or columns that may have been
      # updated, leaving behind things added by e.g. scan-varying refinement
      # such as 'block', the U, B and UB matrices and gradients.
      for key in preds.keys():
        if key in reflections.keys() or key in \
            ['s1', 'xyzcal.mm', 'xyzcal.px', 'entering', 'delpsical.rad']:
          reflections[key] = preds[key]

      # set refinement flags
      assert len(preds) == len(reflections)
      reflections.unset_flags(flex.size_t_range(len(reflections)),
          reflections.flags.excluded_for_refinement |
          reflections.flags.used_in_refinement |
          reflections.flags.centroid_outlier |
          reflections.flags.predicted)
      reflections.set_flags(preds.get_flags(preds.flags.excluded_for_refinement),
          reflections.flags.excluded_for_refinement)
      reflections.set_flags(preds.get_flags(preds.flags.centroid_outlier),
          reflections.flags.centroid_outlier)
      reflections.set_flags(preds.get_flags(preds.flags.used_in_refinement),
          reflections.flags.used_in_refinement)
      reflections.set_flags(preds.get_flags(preds.flags.predicted),
          reflections.flags.predicted)

      logger.info('Saving reflections with updated predictions to {0}'.format(
        params.output.reflections))
      if params.output.include_unused_reflections:
        reflections.as_pickle(params.output.reflections)
      else:
        sel = reflections.get_flags(reflections.flags.used_in_refinement)
        reflections.select(sel).as_pickle(params.output.reflections)

    # For debugging, if requested save matches to file
    if params.output.matches:
      matches = refiner.get_matches()
      logger.info('Saving matches (use for debugging purposes) to {0}'.format(
        params.output.matches))
      matches.as_pickle(params.output.matches)

    # Correlation plot
    if params.output.correlation_plot.filename is not None:
      from os.path import splitext
      root, ext = splitext(params.output.correlation_plot.filename)
      if not ext: ext = ".pdf"

      steps = params.output.correlation_plot.steps
      if steps is None: steps = [history.get_nrows()-1]

      # extract individual column names or indices
      col_select = params.output.correlation_plot.col_select

      num_plots = 0
      for step in steps:
        fname_base = root
        if len(steps) > 1: fname_base += "_step%02d" % step

        corrmats, labels = refiner.get_parameter_correlation_matrix(step, col_select)
        if [corrmats, labels].count(None) == 0:
          from dials.algorithms.refinement.refinement_helpers import corrgram
          for resid_name, corrmat in corrmats.items():
            plot_fname = fname_base + "_" + resid_name + ext
            plt = corrgram(corrmat, labels)
            if plt is not None:
              logger.info('Saving parameter correlation plot to {}'.format(plot_fname))
              plt.savefig(plot_fname)
              plt.close()
              num_plots += 1
          mat_fname = fname_base + ".pickle"
          with open(mat_fname, 'wb') as handle:
            for k, corrmat in corrmats.items():
              corrmats[k] = corrmat.as_scitbx_matrix()
            logger.info('Saving parameter correlation matrices to {0}'.format(mat_fname))
            pickle.dump({'corrmats':corrmats, 'labels':labels}, handle)

      if num_plots == 0:
        msg = "Sorry, no parameter correlation plots were produced. Please set " \
              "track_parameter_correlation=True to ensure correlations are " \
              "tracked, and make sure correlation_plot.col_select is valid."
        logger.info(msg)

    # Write out refinement history, if requested
    if params.output.history:
      with open(params.output.history, 'wb') as handle:
        logger.info('Saving refinement step history to {0}'.format(
          params.output.history))
        pickle.dump(history, handle)

    # Log the total time taken
    logger.info("\nTotal time taken: {0:.2f}s".format(time() - start_time))

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
