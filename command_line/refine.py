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


class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] [param.phil] " \
             "sweep.json crystal.json reflections.pickle"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output filename option
    self.config().add_option(
        '--output-sweep-filename',
        dest = 'output_sweep_filename',
        type = 'string', default = 'refined_sweep.json',
        help = 'Set the filename for refined experimental models.')

    # Output filename option
    self.config().add_option(
        '--output-crystal-filename',
        dest = 'output_crystal_filename',
        type = 'string', default = 'refined_crystal.json',
        help = 'Set the filename for refined crystal model.')

    # Output filename option
    self.config().add_option(
        '--output-reflections-filename',
        dest = 'output_reflections_filename',
        type = 'string', default = None,
        help = 'Set the filename for reflections predicted by the refined model.')

    # Add a verbosity option
    self.config().add_option(
        "-v", "--verbosity",
        action="count", default=0,
        help="set verbosity level; -vv gives verbosity level 2.")

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.algorithms.refinement import RefinerFactory
    from dials.model.serialize import load, dump
    import cPickle as pickle

    # Check the number of arguments is correct
    if len(args) != 3:
      self.config().print_help()
      return

    # Get the refiner
    print 'Configuring refiner'

    # Try to load the models
    print 'Loading models from {0} and {1}'.format(args[0], args[1])
    sweep = load.sweep(args[0])
    crystal = load.crystal(args[1])
    reflections = pickle.load(open(args[2], 'rb'))

    refiner = RefinerFactory.from_parameters_data_models(params,
        reflections, sweep, crystal=crystal, verbosity=options.verbosity)

    # Refine the geometry
    print 'Performing refinement'

    # Refine and get the refinement history
    refined = refiner.run()

    # update the input sweep
    sweep.set_beam(refiner.get_beam())
    sweep.set_detector(refiner.get_detector())
    sweep.set_goniometer(refiner.get_goniometer())

    # Save the refined geometry to file
    output_sweep_filename = options.output_sweep_filename
    print 'Saving refined geometry to {0}'.format(output_sweep_filename)
    dump.sweep(sweep, open(output_sweep_filename, 'w'))

    # Save the refined crystal to file
    output_crystal_filename = options.output_crystal_filename
    print 'Saving refined geometry to {0}'.format(output_crystal_filename)
    dump.crystal(refiner.get_crystal(), open(output_crystal_filename, 'w'))

    # Predict reflections and save to file
    output_reflections_filename = options.output_reflections_filename
    if output_reflections_filename:
      print "Predicting reflections with the refined model"
      rlist = refiner.predict_reflections()
      print "Saving reflections to {0}".format(output_reflections_filename)
      pickle.dump(rlist, open(output_reflections_filename, 'wb'),
          pickle.HIGHEST_PROTOCOL)

    return

if __name__ == '__main__':
  script = Script()
  script.run()
