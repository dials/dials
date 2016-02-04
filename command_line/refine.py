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
from __future__ import division
from copy import deepcopy
from libtbx.utils import Sorry

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
                "supported types"

      save_matrix = False
        .type = bool
        .help = "Save the matrix and column labels in a pickle file for"
                "later inspection, replotting etc."

      col_select = None
        .type = str
        .help = "Specific columns to include in the plots of parameter"
                "correlations, either specifed by parameter name or column"
                "number. Defaults to all columns."
                "This option is useful when there is a large number of"
                "parameters"
        .multiple = True

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

  include scope dials.algorithms.refinement.refiner.phil_scope
''', process_includes=True)

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
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
      phil=phil_scope,
      read_reflections=True,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

  def write_centroids_table(self, refiner, filename):

    matches = refiner.get_matches()

    f = open(filename,"w")
    header = ("H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t"
        "Y_calc\tPhi_calc")
    msg_temp = ("%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%5.3f\t%9.6f")
    has_del_psi = matches.has_key('delpsical.rad')
    if has_del_psi:
      header += "\tDelta_Psi"
      msg_temp += "\t%9.6f"
    header += "\n"
    msg_temp += "\n"
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
    f.close()
    return

  @staticmethod
  def parameter_correlation_plot(corrmat, labels):
    """Create a correlation matrix plot between columns of the Jacobian at
    the specified refinement step. Inspired by R's corrplot and
    https://github.com/louridas/corrplot/blob/master/corrplot.py"""

    try: # is corrmat a scitbx matrix?
      corrmat = corrmat.as_flex_double_matrix()
    except AttributeError: # assume it is already a flex double matrix
      pass
    assert corrmat.is_square_matrix()

    nr = corrmat.all()[0]
    assert nr == len(labels)

    from math import pi, sqrt
    try:
      import matplotlib
      matplotlib.use('Agg')
      import matplotlib.pyplot as plt
      import matplotlib.cm as cm
    except ImportError as e:
      msg = "matplotlib modules not available " + str(e)
      info(msg)
      return None

    plt.figure(1)
    ax = plt.subplot(1, 1, 1, aspect='equal')
    clrmap = cm.get_cmap('bwr')

    for x in xrange(nr):
      for y in xrange(nr):
        d = corrmat[x, y]
        d_abs = abs(d)
        circ = plt.Circle((x, y),radius=0.9*sqrt(d_abs)/2)
        circ.set_edgecolor('white')
        # put data into range [0,1] and invert so that 1 == blue and 0 == red
        facecolor = 1 - (0.5*d + 0.5)
        circ.set_facecolor(clrmap(facecolor))
        ax.add_artist(circ)
    ax.set_xlim(-0.5, nr-0.5)
    ax.set_ylim(-0.5, nr-0.5)

    ax.xaxis.tick_top()
    xtickslocs = range(len(labels))
    ax.set_xticks(xtickslocs)
    ax.set_xticklabels(labels, rotation=30, fontsize='small', ha='left')

    ax.invert_yaxis()
    ytickslocs = range(len(labels))
    ax.set_yticks(ytickslocs)
    ax.set_yticklabels(labels, fontsize='small')

    xtickslocs = [e + 0.5 for e in range(len(labels))]
    ax.set_xticks(xtickslocs, minor=True)
    ytickslocs = [e + 0.5 for e in range(len(labels))]
    ax.set_yticks(ytickslocs, minor=True)
    plt.grid(color='0.8', which='minor', linestyle='-')

    # suppress major tick marks
    ax.tick_params(which='major', width=0)

    # need this otherwise text gets clipped
    plt.tight_layout()

    # FIXME should this also have a colorbar as legend?
    return plt

  @staticmethod
  def check_input(reflections):
    '''Check the input is suitable for refinement. So far just check keys in
    the reflection table. Maybe later check experiments have overlapping models
    etc.'''

    msg = "The supplied reflection table does not have the required data " + \
      "column: {0}"
    for key in ["xyzobs.mm.value", "xyzobs.mm.variance"]:
      if not reflections.has_key(key):
        msg = msg.format(key)
        raise Sorry(msg)

    # FIXME add other things to be checked here
    return

  def run(self):
    '''Execute the script.'''
    from time import time
    import cPickle as pickle
    from logging import info
    from dials.util import log
    from dials.algorithms.refinement import RefinerFactory
    from dials.util.options import flatten_reflections, flatten_experiments

    start_time = time()

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    # Try to load the models and data
    nexp = len(experiments)
    if nexp == 0:
      print "No Experiments found in the input"
      self.parser.print_help()
      return
    if len(reflections) == 0:
      print "No reflection data found in the input"
      self.parser.print_help()
      return
    if len(reflections) > 1:
      raise Sorry("Only one reflections list can be imported at present")
    reflections = reflections[0]

    self.check_input(reflections)

    # Configure the logging
    log.config(params.refinement.verbosity,
      info=params.output.log,
      debug=params.output.debug_log)
    from dials.util.version import dials_version
    info(dials_version())

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      info('The following parameters have been modified:\n')
      info(diff_phil)

    # Modify options if necessary
    if params.output.correlation_plot.filename is not None:
      params.refinement.refinery.track_parameter_correlation = True

    # Get the refiner
    info('Configuring refiner')
    refiner = RefinerFactory.from_parameters_data_experiments(params,
        reflections, experiments)

    # Refine the geometry
    if nexp == 1:
      info('Performing refinement of a single Experiment...')
    else:
      info('Performing refinement of {0} Experiments...'.format(nexp))

    # Refine and get the refinement history
    history = refiner.run()

    if params.output.centroids:
      info("Writing table of centroids to '{0}'".format(
        params.output.centroids))
      self.write_centroids_table(refiner, params.output.centroids)

    # Write scan-varying parameters to file, if there were any
    if params.output.parameter_table:
      scan = refiner.get_scan()
      if scan:
        text = refiner.get_param_reporter().varying_params_vs_image_number(
            scan.get_array_range())
        if text:
          info("Writing scan-varying parameter table to {0}".format(
            params.output.parameter_table))
          f = open(params.output.parameter_table,"w")
          f.write(text)
          f.close()
        else:
          info("No scan-varying parameter table to write")

    # get the refined experiments
    experiments = refiner.get_experiments()
    crystals = experiments.crystals()

    if len(crystals) == 1:
      # output the refined model for information
      info('')
      info('Final refined crystal model:')
      info(crystals[0])

    # Save the refined experiments to file
    output_experiments_filename = params.output.experiments
    info('Saving refined experiments to {0}'.format(output_experiments_filename))
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    dump = ExperimentListDumper(experiments)
    dump.as_json(output_experiments_filename)

    # Save reflections with updated predictions if requested (allow to switch
    # this off if it is a time-consuming step)
    if params.output.reflections:
      # Update predictions for all indexed reflections
      info('Updating predictions for indexed reflections')
      preds = refiner.predict_for_indexed()

      # just copy over the columns of interest, leaving behind things
      # added by e.g. scan-varying refinement such as 'block', the
      # U, B and UB matrices and gradients.
      reflections['s1'] = preds['s1']
      reflections['xyzcal.mm'] = preds['xyzcal.mm']
      reflections['xyzcal.px'] = preds['xyzcal.px']
      if preds.has_key('entering'):
        reflections['entering'] = preds['entering']

      # set centroid_outlier flag
      mask = preds.get_flags(preds.flags.centroid_outlier)
      reflections.set_flags(mask, reflections.flags.centroid_outlier)

      # FIXME redo outlier rejection on the new predictions with updated geometry
      info('Saving reflections with updated predictions to {0}'.format(
        params.output.reflections))
      reflections.as_pickle(params.output.reflections)

    # For debugging, if requested save matches to file
    if params.output.matches:
      matches = refiner.get_matches()
      info('Saving matches (use for debugging purposes) to {0}'.format(
        params.output.matches))
      matches.as_pickle(params.output.matches)

    if params.output.correlation_plot.filename is not None:
      from os.path import splitext
      root, ext = splitext(params.output.correlation_plot.filename)
      if not ext: ext = ".pdf"

      steps = params.output.correlation_plot.steps
      if steps is None: steps = [history.get_nrows()-1]

      # flatten list of column names
      col_select = params.output.correlation_plot.col_select
      if len(col_select) != 0:
        col_select = " ".join(params.output.correlation_plot.col_select).split()
      else: col_select = None
      save_matrix = params.output.correlation_plot.save_matrix

      num_plots = 0
      for step in steps:
        fname_base = root + "_step%02d" % step
        plot_fname = fname_base + ext

        corrmat, labels = refiner.get_parameter_correlation_matrix(step, col_select)
        plt = self.parameter_correlation_plot(corrmat, labels)
        if plt is not None:
          info('Saving parameter correlation plot to {}'.format(plot_fname))
          plt.savefig(plot_fname)
          num_plots += 1

          if save_matrix:
            mat_fname = fname_base + ".pickle"
            with open(mat_fname, 'wb') as handle:
              py_mat = corrmat.as_scitbx_matrix() #convert to pickle-friendly form
              info('Saving parameter correlation matrix to {0}'.format(mat_fname))
              pickle.dump({'corrmat':py_mat, 'labels':labels}, handle)

      if num_plots == 0:
        msg = "Sorry, no parameter correlation plots were produced. Please set " \
              "track_parameter_correlation=True to ensure correlations are " \
              "tracked, and make sure correlation_plot.col_select is valid."
        info(msg)

    # Write out refinement history, if requested
    if params.output.history:
      with open(params.output.history, 'wb') as handle:
        info('Saving refinement step history to {0}'.format(
          params.output.history))
        pickle.dump(history, handle)

    # Log the total time taken
    info("\nTotal time taken: {0:.2f}s".format(time() - start_time))

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
