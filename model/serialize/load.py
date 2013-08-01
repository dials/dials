from __future__ import division
#!/usr/bin/env python
#
# dials.model.serialize.load.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# Import to give access from here
from dxtbx.serialize.load import imageset as sweep
from dxtbx.serialize.load import imageset_from_string as sweep_from_string

def crystal_from_string(string):
    ''' Load the string and return the models.

    Params:
        string The JSON string to load

    Returns:
        The models

    '''
    import json
    from dials.model.serialize.crystal import crystal_from_dict
    return crystal_from_dict(json.loads(string))

def crystal(infile):
    ''' Load the given JSON file.

    Params:
        infile The input filename or file object

    Returns:
        The models

    '''
    # If the input is a string then open and read from that file
    if isinstance(infile, str):
        with open(infile, 'r') as infile:
            return crystal_from_string(infile.read())

    # Otherwise assume the input is a file and read from it
    else:
        return crystal_from_string(infile.read())

def reflections(infile):
    ''' Load the given reflection file.

    Params:
        infile The input filename or file object

    Returns:
        The reflection list

    '''
    import cPickle as pickle

    # If the input is a string then open and read from that file
    if isinstance(infile, str):
        with open(infile, 'rb') as infile:
            return pickle.load(infile)

    # Otherwise assume the input is a file and read from it
    else:
        return pickle.load(infile)
