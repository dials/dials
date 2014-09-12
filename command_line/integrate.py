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

  dials.integrate experiments.json reference=indexed.pickle

  dials.integrate experiments.json reference=indexed.pickle integrated=integrated.pickle

  dials.integrate experiments.json profile.phil

  dials.integrate experiments.json predicted=predicted.pickle reference=indexed.pickle

  dials.integrate experiments.json shoeboxes=shoeboxes.dat

'''

class Script(object):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # Set the phil scope
    phil_scope = parse('''
      
      output {
        profile_model = 'profile_model.phil'
          .type = str
          .help = "The profile parameters output filename"

        reflections = 'integrated.pickle'
          .type = str
          .help = "The integrated output filename"
      }

      input {
        shoeboxes = None
          .type = str
          .help = "The shoebox input filename"
      }

      include scope dials.algorithms.integration.integrator.phil_scope
      include scope dials.algorithms.profile_model.profile_model.phil_scope

    ''', process_includes=True)

    # The script usage
    usage  = "usage: %prog [options] experiment.json"

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_reflections=True,
      read_experiments=True)

  def run(self):
    ''' Perform the integration. '''
    from time import time

    # Check the number of arguments is correct
    start_time = time()

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    
    # Get the expected data
    assert(len(params.input.experiments) == 1)
    exlist = params.input.experiments[0].data
    assert(len(exlist) == 1)
    if len(params.input.reflections) > 0:
      assert(len(params.input.reflections) == 1)
      reference = self.process_reference(params.input.reflections[0].data)
    else:
      reference = None

    # Load the data
    shoeboxes = None

    # Initialise the integrator
    if None in exlist.goniometers():
      from dials.algorithms.integration import IntegratorStills
      integrator = IntegratorStills(params, exlist, reference, None, shoeboxes)
    else:
      from dials.algorithms.integration import Integrator
      integrator = Integrator(params, exlist, reference, None, shoeboxes)

    # Integrate the reflections
    reflections = integrator.integrate()

    # Save the reflections
    self.save_reflections(reflections, params.output.reflections)

    # Print the total time taken
    print "\nTotal time taken: ", time() - start_time

  def process_reference(self, reference):
    ''' Load the reference spots. '''
    from dials.util.command_line import Command
    from dials.array_family import flex
    assert("miller_index" in reference)
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
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
