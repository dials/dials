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
from dials.util.script import ScriptRunner
from dials.util.command_line import Importer
from dials.model.experiment.experiment_list import \
  ExperimentList, ExperimentListDumper

class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] [param.phil] " \
             "experiments.json reflections.pickle"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output experiments filename option
    self.config().add_option(
        '--output-experiments-filename',
        dest = 'output_experiments_filename',
        type = 'string', default = 'refined_experiments.json',
        help = 'Set the filename for refined experimental models.')

    self.config().add_option(
        '--output-centroids-filename',
        dest = 'output_centroids_filename',
        type = 'string',
        help = 'Set the filename for the table of centroids at the end of refinement.')

    self.config().add_option(
        '--output-parameters-filename',
        dest = 'output_parameters_filename',
        type = 'string',
        help = 'Set the filename for the table of scan varying parameter values'
               ' at the end of refinement.')

    # Add a verbosity option
    self.config().add_option(
        "-v", "--verbosity",
        action="count", default=0,
        help="set verbosity level; -vv gives verbosity level 2.")

  def write_residuals_table(self, refiner, filename):

    matches = refiner.get_matches()

    f = open(filename,"w")
    header = ("H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t"
        "Y_calc\tPhi_calc\n")
    f.write(header)

    for m in matches:
      msg = ("%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%9.6f\t"
            "%5.3f\n")
      msg = msg % (m.miller_index[0], m.miller_index[1], m.miller_index[2],
                   m.frame_obs, m.x_obs, m.y_obs,
                   m.phi_obs, m.x_calc, m.y_calc, m.phi_calc)
      f.write(msg)
    f.close()

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.algorithms.refinement import RefinerFactory
    import cPickle as pickle

    # Check the number of arguments
    if len(args) < 2:
      self.config().print_help()
      return

    importer = Importer(args, check_format=False, verbose=True)

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
        reflections, experiments, verbosity=options.verbosity)

    # Refine the geometry
    print 'Performing refinement'

    # Refine and get the refinement history
    refined = refiner.run()

    if options.output_centroids_filename:
      print "Writing table of centroids to '{0}'".format(
        options.output_centroids_filename)
      self.write_residuals_table(refiner, options.output_centroids_filename)

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

    return

if __name__ == '__main__':
  script = Script()
  script.run()
