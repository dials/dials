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


class Script(object):
  ''' Encapsulate the script in a class. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # The phil parameters
    phil_scope = parse('''
      output = profile.phil
        .type = str
        .help = "The filename for the profile parameters"

      include scope dials.algorithms.profile_model.factory.phil_scope
    ''', process_includes=True)

    # The script usage
    usage  = "usage: %prog [options] experiments.json spots.pickle"
    self.parser = OptionParser(
      usage=usage,
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

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    if len(reflections) == 0 and len(experiments) == 0:
      self.parser.print_help()
      return
    assert(len(reflections) == 1)
    reflections = reflections[0]
    assert(len(experiments) > 0)

    from dials.array_family import flex
    Command.start('Removing invalid coordinates')
    xyz = reflections['xyzcal.mm']
    mask = flex.bool([x == (0, 0, 0) for x in xyz])
    reflections.del_selected(mask)
    Command.end('Removed invalid coordinates, %d remaining' % len(reflections))

    # Create the profile model
    profile_model = ProfileModelFactory.create(params, experiments, reflections)
    for model in profile_model:
      sigma_b = model.sigma_b(deg=True)
      sigma_m = model.sigma_m(deg=True)
      print 'Sigma B: %f' % sigma_b
      print 'Sigma M: %f' % sigma_m

    # Wrtie the parameters
    Command.start("Writing profile model to %s" % params.output)
    with open(params.output, "w") as outfile:
      outfile.write(profile_model.dump().as_str())
    Command.end("Wrote profile model to %s" % params.output)


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
