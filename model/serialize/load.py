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

def load_experimental_models_from_string(string):
    ''' Load the experimental models from a JSON string

    Params:
        string The JSON string

    Returns:
        (beam, detector, goniometer, scan)

    '''
    import json
    from dials.model.serialize.beam import beam_from_dict
    from dials.model.serialize.detector import detector_from_dict
    from dials.model.serialize.goniometer import goniometer_from_dict
    from dials.model.serialize.scan import scan_from_dict

    # Convert all string to dictionaries
    models = json.loads(string)

    # Return a tuple of the models
    return (beam_from_dict(models['beam']),
            detector_from_dict(models['detector']),
            goniometer_from_dict(models['goniometer']),
            scan_from_dict(models['scan']))

def load_crystal_model_from_string(string):
    ''' Load the crystal model from a JSON string

    Params:
        string The JSON string

    Returns:
        The crystal model

    '''
    import json
    from dials.model.serialize.crystal import crystal_from_dict
    return crystal_from_dict(json.loads(string))

def loads(string, model):
    ''' Load the string and return the models.

    Params:
        string The JSON string to load

    Returns:
        The models

    '''
    from dials.model.experiment import Beam, Detector, Goniometer, Scan
    from dials.model.experiment.crystal_model.crystal import Crystal

    # If this is a tuple of 4 elements
    if type(obj) == tuple and len(obj) == 4:

        # If the types for experimental models
        if (isinstance(obj[0], Beam) and
            isinstance(obj[1], Detector) and
            isinstance(obj[2], Goniometer) and
            isinstance(obj[3], Scan):
            return dump_experimental_models_to_string(*obj)

    # If this is a crystal model then dump
    elif isinstance(obj, Crystal):
        return dump_crystal_model_to_string(obj)

    # Otherwise raise an exception
    else:
        raise TypeError("Unknown type to serialize")

def load(infile):
    ''' Load the given JSON file.

    Params:
        infile The input filename or file object

    Returns:
        The models

    '''
    # If the input is a string then open and read from that file
    if isinstance(infile, str):
        with infile = open(infile, 'r'):
            return loads(infile.read())

    # Otherwise assume the input is a file and read from it
    else:
        return loads(infile.read())
