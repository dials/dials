#!/usr/bin/env python
#
# dials.integrate.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

help_message = '''

This program is used to integrate the reflections on the diffraction images. It
is called with an experiment list outputted from dials.index or dials.refine.
The extend of the shoeboxes is specified through the profile parameters
shoebox.sigma_b and shoebox.sigma_m (use the --show-config option for more
details). These parameters can be specified directly, otherwise a set of strong
indexed reflections are needed to form the profile model; these are specified
using the -r (for reference) option. The program can also be called with a
specific set of predictions using the -p option.

Once a profile model is given and the size of the measurement boxes have been
calculated, the program will extract the reflections to file. The reflections
will then be integrated. The reflections can be integrated with different
options using the same measurement boxes by giving the measurement box file
using the -s option. This will skip reading the measurement boxes and go
directly to integrating the reflections.

Examples:

  dials.integrate experiments.json -r indexed.pickle

  dials.integrate experiments.json -r indexed.pickle -o integrated.pickle

  dials.integrate experiments.json shoebox.sigma_b=0.024 shoebox.sigma_m=0.044

  dials.integrate experiments.json -p predicted.pickle -r indexed.pickle

  dials.integrate experiments.json -s shoeboxes.dat

'''

class Script(ScriptRunner):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] experiment.json"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage,
                          epilog=help_message,
                          home_scope="integration")

    # Output filename option
    self.config().add_option(
      '-o', '--output',
      dest = 'output',
      type = 'string', default = 'integrated.pickle',
      help = 'Set the filename for integrated reflections.')

    # The predicted reflections to integrate
    self.config().add_option(
      '-p', '--predicted',
      dest = 'predicted',
      type = 'string', default = None,
      help = 'Specify predicted reflections.')

    # The intermediate shoebox data
    self.config().add_option(
      '-r', '--reference',
      dest = 'reference',
      type = 'string', default = None,
      help = 'Specify reference reflections.')

    # The intermediate shoebox data
    self.config().add_option(
      '-s', '--shoeboxes',
      dest = 'shoeboxes',
      type = 'string', default = None,
      help = 'Specify shoeboxes to integrate.')

  def main(self, params, options, args):
    ''' Perform the integration. '''
    from time import time

    # Check the number of arguments is correct
    start_time = time()

    # Check the number of command line arguments
    if len(args) != 1:
      self.config().print_help()
      return

    # Print the diff phil
    diff_phil_str = self.config().diff_phil().as_str()
    print 'Integrating with the following user specified parameters:\n'
    if (diff_phil_str is not ''):
      print diff_phil_str
    else:
      print 'All parameters set to defaults'

    # Load the experiment list
    exlist = self.load_experiments(args[0])

    # Load the data

    shoeboxes = reference = predicted = None

    if options.shoeboxes:
      shoeboxes = options.shoeboxes
    if options.reference:
      reference = self.load_reference(options.reference)
    if options.predicted:
      predicted = self.load_predicted(options.predicted)

    # Initialise the integrator
    if None in exlist.goniometers():
      from dials.algorithms.integration import IntegratorStills
      integrator = IntegratorStills(params, exlist, reference, predicted, shoeboxes)
    else:
      from dials.algorithms.integration import Integrator
      integrator = Integrator(params, exlist, reference, predicted, shoeboxes)

    # Integrate the reflections
    reflections = integrator.integrate()

    # Save the reflections
    self.save_reflections(reflections, options.output)

    # Print the total time taken
    print "\nTotal time taken: ", time() - start_time

  def load_predicted(self, filename):
    ''' Load the predicted reflections. '''
    from dials.array_family import flex
    if filename is not None:
      return flex.reflection_table.from_pickle(filename)
    return None

  def load_experiments(self, filename):
    ''' Load the experiment list. '''
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from dials.util.command_line import Command
    Command.start('Loading experiments from %s' % filename)
    exlist = ExperimentListFactory.from_json_file(filename)
    Command.end('Loaded experiments from %s' % filename)
    if len(exlist) == 0:
      raise RuntimeError('experiment list is empty')
    elif len(exlist.imagesets()) > 1 or len(exlist.detectors()) > 1:
      raise RuntimeError('experiment list contains > 1 imageset or detector')
    return exlist

  def load_reference(self, filename):
    ''' Load the reference spots. '''
    from dials.util.command_line import Command
    from dials.array_family import flex

    if not filename:
      return None

    Command.start('Loading reference spots from %s' % filename)
    reference = flex.reflection_table.from_pickle(filename)
    assert("miller_index" in reference)
    Command.end('Loaded reference spots from %s' % filename)
    Command.start('Removing reference spots with invalid coordinates')
    mask = flex.bool([x == (0, 0, 0) for x in reference['xyzcal.mm']])
    reference.del_selected(mask)
    mask = flex.bool([h == (0, 0, 0) for h in reference['miller_index']])
    reference.del_selected(mask)
    Command.end('Removed reference spots with invalid coordinates, %d remaining' %
                len(reference))
    return reference

  def save_reflections(self, reflections, filename):
    ''' Save the reflections to file. '''
    from dials.util.command_line import Command
    Command.start('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)
    Command.end('Saved %d reflections to %s' % (len(reflections), filename))


if __name__ == '__main__':
  script = Script()
  script.run()
