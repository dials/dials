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
from dxtbx.model.experiment.experiment_list import ExperimentListDumper

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

        correlation_plot_filename = None
          .type = str
          .help = "The filename of output of a plot of parameter correlations"

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
    dump = ExperimentListDumper(experiments)
    dump.as_json(output_experiments_filename)

    # Write out refined reflections, if requested
    if params.output.reflections_filename:
      matches = refiner.get_matches()
      print 'Saving refined reflections to {0}'.format(
        params.output.reflections_filename)
      matches.as_pickle(params.output_reflections_filename)

    if params.output.correlation_plot_filename:
      if refined.parameter_correlation:
        plt = refiner.parameter_correlation_plot(len(refined.parameter_correlation)-1)
        plt.tight_layout()
        plt.savefig(params.output.correlation_plot_filename)
      else:
        print "Sorry, no parameter correlations were tracked. Please set " \
              "track_parameter_correlation=True"

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
