#!/usr/bin/env python
#
# reciprocal_space.py
#
#  Copyright (C) 2014 Diamond Light Source, James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Script(object):
  ''' The class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # Create the phil parameters
    phil_scope = parse('''

      output = simulated.pickle
        .type = str
        .help = "The output pickle file"

      num = 1000
        .type = int
        .help = "The number of reflections"

      intensity = 1000
        .type = int
        .help = "The intensity of reflections"

      background = 0
        .type = str
        .help = "The background level"

      random = False
        .type = bool
        .help = "Random intensities and background"

      include scope dials.algorithms.profile_model.profile_model.phil_scope

    ''', process_includes=True)

    # Create the option parser
    usage = "usage: %s [options] experiment.json" \
      % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope)

  def run(self):
    ''' Run the script. '''
    from dials.algorithms.simulation.reciprocal_space import Simulator
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from dials.util.command_line import Command
    from libtbx.utils import Sorry
    from math import pi

    # Parse the command line arguments
    params, options, args = self.parser.parse_args()

    # Ensure we have enough arguments
    if len(args) != 1:
      self.parser.print_help()
      exit(0)

    # Check we have some profile parameters
    if len(params.profile) != 1:
      raise Sorry('no profile parameters specified')

    print 'Simulating with the following parameters:'
    print ' # Reflections: %d' % params.num
    print ' Intensity: %d' % params.intensity
    print ' Background: %s' % params.background
    print ' Sigma B: %f degrees' % params.profile[0].sigma_b
    print ' Sigma M: %f degrees' % params.profile[0].sigma_m
    print ' N Sigma: %f degrees' % params.profile[0].n_sigma
    print ' Random: %s' % str(params.random)

    # Get the experiments
    experiments = ExperimentListFactory.from_json_file(args[0], check_format=False)
    if len(experiments) != 1:
      raise Sorry('experiment list must contain exactly 1 experiment')
    experiment = experiments[0]

    # Do the simulation
    sigma_b = params.profile[0].sigma_b * pi / 180
    sigma_m = params.profile[0].sigma_m * pi / 180
    n_sigma = params.profile[0].n_sigma
    N = params.num
    I = params.intensity
    B = map(int, params.background.split(","))
    B = B + [0] * (4 - len(B))
    simulate = Simulator(experiment, sigma_b, sigma_m, n_sigma)
    if params.random:
      refl = simulate.with_random_intensity(N, I, B[0], B[1], B[2], B[3])
    else:
      refl = simulate.with_given_intensity(N, I, B[0], B[1], B[2], B[3])

    # Save the reflections to file
    Command.start('Writing reflections to %s' % params.output)
    refl.as_pickle(params.output)
    Command.end('Write reflections to %s' % params.output)


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
