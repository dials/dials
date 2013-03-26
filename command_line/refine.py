#!/usr/bin/env python
#
# refine.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

def read_reflection_file(filename):
    '''Read reflections from pickle file.'''
    from dials.model.data import ReflectionList
    import pickle
    return pickle.load(open(reflection_file, 'rb'))
    
def run(reflection_file):
    '''Do the refinement.'''

    # Load the reflections from the pickle file
    reflections = read_reflection_file(reflection_file)    

    # Print some reflections
    for r in reflections:
        print r

if __name__ == '__main__':
    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/data.p"
    
    # Parse the command line options
    parser = OptionParser(usage)
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) < 1:
        print parser.print_help()
    else:
 
        # Get stuff from args
        reflection_file = args[0]

        # Run the refinement
        run(reflection_file)
