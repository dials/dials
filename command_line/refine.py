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
from dials.util.command_line import Importer
from dxtbx.model.experiment.experiment_list import \
  ExperimentList, ExperimentListDumper

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # The script usage
    usage  = "usage: %prog [options] [param.phil] " \
             "experiments.json reflections.pickle"

    # The phil scope
    phil_scope = parse('''
      include scope dials.data.refinement.phil_scope
    ''', process_includes=True)

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope)

    # Output experiments filename option
    self.parser.add_option(
        '--output-experiments-filename',
        dest = 'output_experiments_filename',
        type = 'string', default = 'refined_experiments.json',
        help = 'Set the filename for refined experimental models.')

    self.parser.add_option(
        '--output-centroids-filename',
        dest = 'output_centroids_filename',
        type = 'string',
        help = 'Set the filename for the table of centroids at the end of refinement.')

    self.parser.add_option(
        '--output-parameters-filename',
        dest = 'output_parameters_filename',
        type = 'string',
        help = 'Set the filename for the table of scan varying parameter values'
               ' at the end of refinement.')

    self.parser.add_option(
        '--output-correlation-plot-filename',
        dest = 'output_corrplot_filename',
        type = 'string',
        help = 'Set the filename for output of a plot of parameter correlations.')

    self.parser.add_option(
        '--output-reflections-filename',
        dest = 'output_reflections_filename',
        type = 'string',
        help = 'Set the filename for output of refined reflections.')


  def write_centroids_table(self, refiner, filename):

    matches = refiner.get_matches()

    f = open(filename,"w")
    header = ("H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t"
        "Y_calc\tPhi_calc\n")
    f.write(header)

    for m in matches:
      msg = ("%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%5.3f\t"
            "%9.6f\n")
      (h, k, l) = m['miller_index']
      frame = m['xyzobs.px.value'][2]
      x_obs, y_obs, phi_obs = m['xyzobs.mm.value']
      x_calc, y_calc, phi_calc = m['xyzcal.mm']
      msg = msg % (h, k, l,
                   frame, x_obs, y_obs, phi_obs,
                   x_calc, y_calc, phi_calc)
      f.write(msg)
    f.close()
    return

  def run(self):
    '''Execute the script.'''
    from dials.algorithms.refinement import RefinerFactory
    import cPickle as pickle

    # Parse the command line
    params, options, args = self.parser.parse_args(show_diff_phil=True)

    # Check the number of arguments
    if len(args) < 2:
      self.parser.print_help()
      return

    importer = Importer(args, check_format=False, verbose=False)

    # Try to load the models and data
    experiments = importer.experiments
    if experiments is None:
      raise RuntimeError("No Experiments found in the input")
    reflections = importer.reflections
    if reflections is None:
      raise RuntimeError("No reflection data found in the input")
    if len(reflections) > 1:
      raise RuntimeError("Only one reflections list can be imported at present")
    reflections = importer.reflections[0]

    # Get the refiner
    print 'Configuring refiner'
    refiner = RefinerFactory.from_parameters_data_experiments(params,
        reflections, experiments)

    # Refine the geometry
    print 'Performing refinement'

    # Refine and get the refinement history
    refined = refiner.run()

    if options.output_centroids_filename:
      print "Writing table of centroids to '{0}'".format(
        options.output_centroids_filename)
      self.write_centroids_table(refiner, options.output_centroids_filename)

    # Write scan-varying parameters to file, if there were any
    if options.output_parameters_filename:
      scan = refiner.get_scan()
      if scan:
        text = refiner.get_param_reporter().varying_params_vs_image_number(
            scan.get_array_range())
        if text:
          print "Writing scan-varying parameter table to '{0}'".format(
            options.output_parameters_filename)
          f = open(options.output_parameters_filename,"w")
          f.write(text)
          f.close()
        else:
          print "No scan-varying parameter table to write"

    # get the refined experiments
    experiments = refiner.get_experiments()

    # Save the refined experiments to file
    output_experiments_filename = options.output_experiments_filename
    print 'Saving refined experiments to {0}'.format(output_experiments_filename)
    dump = ExperimentListDumper(experiments)
    dump.as_json(output_experiments_filename)

    # Write out refined reflections, if requested
    if options.output_reflections_filename:
      matches = refiner.get_matches()
      print 'Saving refined reflections to {0}'.format(
        options.output_reflections_filename)
      matches.as_pickle(options.output_reflections_filename)

    if options.output_corrplot_filename:
      if refined.parameter_correlation:
        plt = refiner.parameter_correlation_plot(len(refined.parameter_correlation)-1)
        plt.tight_layout()
        plt.savefig(options.output_corrplot_filename)
      else:
        print "Sorry, no parameter correlations were tracked. Please set " \
              "track_parameter_correlation=True"

    return

if __name__ == '__main__':
  script = Script()
  script.run()
