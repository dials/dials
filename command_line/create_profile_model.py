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
    ''')

    # The script usage
    usage  = "usage: %prog [options] experiments.json spots.pickle"
    self.parser = OptionParser(usage=usage, phil=phil_scope)

  def run(self):
    ''' Run the script. '''
    from dials.util.command_line import Importer
    from dials.algorithms.profile_model.profile_model import ProfileModelList
    from dials.util.command_line import Command
    from dials.array_family import flex

    # Parse the command line
    params, options, args = self.parser.parse_args(show_diff_phil=True)

    # Import the items
    importer = Importer(args, check_format=False)
    experiments = importer.experiments
    if experiments is None or len(experiments) != 1:
      self.parser.print_help()
      exit(0)
    if importer.reflections is None or len(importer.reflections) != 1:
      self.parser.print_help()
      exit(0)
    reflections = importer.reflections[0]

    from dials.array_family import flex
    Command.start('Removing invalid coordinates')
    xyz = reflections['xyzcal.mm']
    mask = flex.bool([x == (0, 0, 0) for x in xyz])
    reflections.del_selected(mask)
    Command.end('Removed invalid coordinates, %d remaining' % len(reflections))

    # Create the profile model
    profile_model = ProfileModelList.compute(experiments, reflections)
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
