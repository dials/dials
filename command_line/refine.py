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

from __future__ import division

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # The phil scope
    phil_scope = parse('''

      output {
        experiments_filename = refined_experiments.json
          .type = str
          .help = "The filename for refined experimental models"

        centroids_filename = None
          .type = str
          .help = "The filename for the table of centroids at the end of"
                  "refinement"

        parameters_filename = None
          .type = str
          .help = "The filename for the table of scan varying parameter values"

        correlation_plot {
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

        reflections_filename = None
          .type = str
          .help = "The filename for output of refined reflections"
      }

      include scope dials.algorithms.refinement.refiner.phil_scope
    ''', process_includes=True)

    # The script usage
    usage  = "usage: %prog [options] [param.phil] " \
             "experiments.json reflections.pickle"

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True,
      read_experiments=True,
      check_format=False)

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
      print "matplotlib modules not available", e
      return None

    plt.figure(1)
    ax = plt.subplot(1, 1, 1, aspect='equal')
    poscm = cm.get_cmap('Blues')
    negcm = cm.get_cmap('Reds')

    for x in xrange(nr):
      for y in xrange(nr):
        d = corrmat[x, y]
        clrmap = poscm if d >= 0 else negcm
        d_abs = abs(d)
        circ = plt.Circle((x, y),radius=0.9*sqrt(d_abs)/2)
        circ.set_edgecolor('white')
        circ.set_facecolor(clrmap(d_abs))
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

  def run(self):
    '''Execute the script.'''
    from dials.algorithms.refinement import RefinerFactory
    from dials.util.options import flatten_reflections, flatten_experiments

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    # Try to load the models and data
    if len(experiments) == 0:
      raise RuntimeError("No Experiments found in the input")
    if len(reflections) == 0:
      raise RuntimeError("No reflection data found in the input")
    if len(reflections) > 1:
      raise RuntimeError("Only one reflections list can be imported at present")
    reflections = reflections[0]

    # Get the refiner
    print 'Configuring refiner'
    refiner = RefinerFactory.from_parameters_data_experiments(params,
        reflections, experiments)

    # Refine the geometry
    print 'Performing refinement'

    # Refine and get the refinement history
    refined = refiner.run()

    if params.output.centroids_filename:
      print "Writing table of centroids to '{0}'".format(
        params.output.centroids_filename)
      self.write_centroids_table(refiner, params.output.centroids_filename)

    # Write scan-varying parameters to file, if there were any
    if params.output.parameters_filename:
      scan = refiner.get_scan()
      if scan:
        text = refiner.get_param_reporter().varying_params_vs_image_number(
            scan.get_array_range())
        if text:
          print "Writing scan-varying parameter table to '{0}'".format(
            params.output.parameters_filename)
          f = open(params.output.parameters_filename,"w")
          f.write(text)
          f.close()
        else:
          print "No scan-varying parameter table to write"

    # get the refined experiments
    experiments = refiner.get_experiments()

    # Save the refined experiments to file
    output_experiments_filename = params.output.experiments_filename
    print 'Saving refined experiments to {0}'.format(output_experiments_filename)
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    dump = ExperimentListDumper(experiments)
    dump.as_json(output_experiments_filename)

    # Write out refined reflections, if requested
    if params.output.reflections_filename:
      matches = refiner.get_matches()
      print 'Saving refined reflections to {0}'.format(
        params.output.reflections_filename)
      matches.as_pickle(params.output.reflections_filename)

    if params.output.correlation_plot.filename is not None:
      from os.path import splitext
      root, ext = splitext(params.output.correlation_plot.filename)
      if not ext: ext = ".pdf"

      steps = params.output.correlation_plot.steps
      if steps is None: steps = [refined.get_nrows()-1]

      # flatten list of column names
      col_select = params.output.correlation_plot.col_select
      if len(col_select) != 0:
        col_select = " ".join(params.output.correlation_plot.col_select).split()
      else: col_select = None
      save_matrix = params.output.correlation_plot.save_matrix
      if save_matrix: import cPickle as pickle

      num_plots = 0
      for step in steps:
        fname_base = root + "_step%02d" % step
        plot_fname = fname_base + ext

        corrmat, labels = refiner.get_parameter_correlation_matrix(step, col_select)
        plt = self.parameter_correlation_plot(corrmat, labels)
        if plt is not None:
          plt.savefig(plot_fname)
          num_plots += 1

          if save_matrix:
            mat_fname = fname_base + ".pickle"
            with open(mat_fname, 'wb') as handle:
              py_mat = corrmat.as_scitbx_matrix() #convert to pickle-friendly form
              pickle.dump({'corrmat':py_mat, 'labels':labels}, handle)

      if num_plots == 0:
        print "Sorry, no parameter correlation plots were produced. Please set " \
              "track_parameter_correlation=True to ensure correlations are " \
              "tracked, and make sure correlation_plot.col_select is valid."

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
