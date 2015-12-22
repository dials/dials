#!/usr/bin/env python
#
# dials.create_profile_model.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from libtbx.phil import parse

help_message = '''

This program computes the profile model from the input reflections. It then
saves a modified experiments.json file with the additional profile model
information. Usually this is performed during integration; however, on some
occasions it may be desirable to compute the profile model independantly.

Examples::

  dials.create_profile_model experiments.json reflections.pickle

'''

phil_scope = parse('''
  output = experiments_with_profile_model.json
    .type = str
    .help = "The filename for the experiments"

  include scope dials.algorithms.profile_model.factory.phil_scope
''', process_includes=True)

class Script(object):
  ''' Encapsulate the script in a class. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] experiments.json spots.pickle" \
      % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      epilog=help_message,
      phil=phil_scope,
      read_reflections=True,
      read_experiments=True,
      check_format=False)

  def run(self):
    ''' Run the script. '''
    from dials.algorithms.profile_model.factory import ProfileModelFactory
    from dials.util.command_line import Command
    from dials.array_family import flex
    from dials.util.options import flatten_reflections, flatten_experiments
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    from libtbx.utils import Sorry

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    if len(reflections) == 0 and len(experiments) == 0:
      self.parser.print_help()
      return
    if len(reflections) != 1:
      raise Sorry('exactly 1 reflection table must be specified')
    if len(experiments) == 0:
      raise Sorry('no experiments were specified')
    reflections = reflections[0]

    from dials.array_family import flex
    Command.start('Removing invalid coordinates')
    xyz = reflections['xyzcal.mm']
    mask = flex.bool([x == (0, 0, 0) for x in xyz])
    reflections.del_selected(mask)
    Command.end('Removed invalid coordinates, %d remaining' % len(reflections))

    # For some reason, some input files have id as type int rather than uint
    reflections['id'] = flex.size_t(list(reflections['id']))

    # Create the profile model
    experiments = ProfileModelFactory.create(params, experiments, reflections)
    for model in experiments:
      sigma_b = model.profile.sigma_b(deg=True)
      sigma_m = model.profile.sigma_m(deg=True)
      if type(sigma_b) == type(1.0):
        print 'Sigma B: %f' % sigma_b
        print 'Sigma M: %f' % sigma_m
      else: # scan varying
        mean_sigma_b = sum(sigma_b) / len(sigma_b)
        mean_sigma_m = sum(sigma_m) / len(sigma_m)
        print 'Sigma B: %f' % mean_sigma_b
        print 'Sigma M: %f' % mean_sigma_m

    # Wrtie the parameters
    Command.start("Writing experiments to %s" % params.output)
    dump = ExperimentListDumper(experiments)
    with open(params.output, "w") as outfile:
      outfile.write(dump.as_json())
    Command.end("Wrote experiments to %s" % params.output)


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
