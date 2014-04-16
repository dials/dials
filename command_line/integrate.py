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


class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] [param.phil] " \
             "sweep.json crystal.json [reference.pickle]"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output filename option
    self.config().add_option(
        '-o', '--output-filename',
        dest = 'output_filename',
        type = 'string', default = 'integrated.pickle',
        help = 'Set the filename for integrated reflections.')

    # Option to save profiles with reflection data
    self.config().add_option(
        '-p', '--save-profiles',
        dest = 'save_profiles',
        action = 'store_true', default = False,
        help = 'Output profiles with reflection data.')

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.algorithms.integration import Integrator
    from dials.algorithms import shoebox
    from dials.model.serialize import load, dump
    from dials.util.command_line import Command
    from dials.util.command_line import Importer
    from time import time

    # Check the number of arguments is correct
    start_time = time()

    # Process the command line options
    Command.start('Processing command line options')
    importer = Importer(args)
    Command.end('Processed command line options')

    # Check the unhandled arguments
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg

    # Check the number of experiments
    if importer.experiments is None or len(importer.experiments) == 0:
      print 'Error: no experiment list specified'
      return
    elif len(importer.experiments) != 1:
      print 'Error: only 1 experiment can be processed at a time'
      return

    if importer.reflections is not None:
      assert(len(importer.reflections) == 1)
      reference = importer.reflections[0]

      from dials.array_family import flex
      Command.start('Removing invalid coordinates')
      xyz = reference['xyzcal.mm']
      mask = flex.bool([x == (0, 0, 0) for x in xyz])
      reference.del_selected(mask)
      Command.end('Removed invalid coordinates, %d remaining' % len(reference))
    else:
      reference = None

    # Get the reference and extracted stuff
    extracted = importer.extracted

    # Get the integrator from the input parameters
    print 'Configurating integrator from input parameters'
    integrate = Integrator(params.integration.shoebox.n_sigma,
                           params.integration.shoebox.n_blocks,
                           params.integration.filter.by_zeta)

    # Intregate the sweep's reflections
    print 'Integrating reflections'
    reflections = integrate(importer.experiments,
        reference=reference, extracted=extracted,
        save_profiles=options.save_profiles)

    # Save the reflections to file
    nvalid = len(reflections)
    Command.start('Saving {0} reflections to {1}'.format(
        nvalid, options.output_filename))
    reflections.as_pickle(options.output_filename)
    Command.end('Saved {0} reflections to {1}'.format(
        nvalid, options.output_filename))

    # Print the total time taken
    print "\nTotal time taken: ", time() - start_time

if __name__ == '__main__':
  script = Script()
  script.run()
