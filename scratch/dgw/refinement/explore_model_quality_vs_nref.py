#!/usr/bin/env python
#
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

'''This script has the same interface as dials.refine2. It performs refinement
jobs, altering the number of input reflections each time. This allows an
investigation into the final model quality and execution time vs the number
of reflections in refinement.'''

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

        # Add a verbosity option
        self.config().add_option(
            "-v", "--verbosity",
            action="count", default=0,
            help="set verbosity level; -vv gives verbosity level 2")

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.algorithms.refinement import RefinerFactory
        from dials.model.serialize import load, dump
        from dials.model.data import ReflectionList
        import cPickle as pickle
        from time import time
        import sys

        # Check the number of arguments is correct
        if len(args) != 3:
            self.config().print_help()
            return

        # Hack the phil parameters to control the number of reflections
        params.refinement.reflections.use_all_reflections = False
        params.refinement.reflections.random_seed = None

        # and also to make the target achieved criterion more stringent
        params.refinement.target.bin_size_fraction = 0.1

        # print header for output table
        print "Nref RMSDx RMSDy RMSDphi Steps Target_acheived Time"

        #for n in range(1,170,2):
        for n in [e * 0.1 for e in range(1,100)]:

            # Set nref to use
            params.refinement.reflections.reflections_per_degree = n

            # Get the refiner from the input parameters
            refine = RefinerFactory.from_parameters(params, options.verbosity)

            # Try to load the models
            sweep = load.sweep(args[0])
            crystal = load.crystal(args[1])
            reflections = pickle.load(open(args[2], 'rb'))

            # Refine the geometry
            start_time = time()
            try:
                refined = refine(sweep, crystal, reflections)
            except Exception as e:
                print "ERROR occurred"
                continue
            time_taken = time() - start_time

            print refined.history.num_reflections[-1],
            print "%.6f %.6f %.8f" % refine.rmsds(),
            print "%d" % refined.get_num_steps(),
            print refined.test_for_termination(),
            print "%.3f" % time_taken

            # flush the buffer so we can see some output
            sys.stdout.flush()

        return


if __name__ == '__main__':
    script = Script()
    script.run()
