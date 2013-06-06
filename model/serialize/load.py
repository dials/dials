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

def model_from_dict(obj):
    ''' Convert the dictionary into a set of models.

    Params:
        obj The input dictionary

    Returns:
        A set of models in a named tuple

    '''
    from collections import namedtuple
    from dials.model.serialize.beam import beam_from_dict
    from dials.model.serialize.detector import detector_from_dict
    from dials.model.serialize.goniometer import goniometer_from_dict
    from dials.model.serialize.scan import scan_from_dict
    from dials.model.serialize.crystal import crystal_from_dict

    # Create the named tuple type
    Model = namedtuple('Model', ['filenames', 'beam', 'detector',
                                 'goniometer', 'scan', 'crystal'])

    # Convert all dictionaries to models
    if 'filenames' in obj:
        filenames = str(obj['filenames'])
    else:
        filenames = None

    if 'beam' in obj:
        beam = beam_from_dict(obj['beam'])
    else:
        beam = None

    if 'detector' in obj:
        detector = detector_from_dict(obj['detector'])
    else:
        detector = None

    if 'goniometer' in obj:
        goniometer = goniometer_from_dict(obj['goniometer'])
    else:
        goniometer = None
    if 'scan' in obj:
        scan = scan_from_dict(obj['scan'])
    else:
        scan = None

    if 'crystal' in obj:
        crystal = crystal_from_dict(obj['crystal'])
    else:
        crystal = None

    # Return models in a named tuple
    return Model(filenames=filenames,
                 beam=beam,
                 detector=detector,
                 goniometer=goniometer,
                 scan=scan,
                 crystal=crystal)

def loads(string):
    ''' Load the string and return the models.

    Params:
        string The JSON string to load

    Returns:
        The models

    '''
    import json
    return model_from_dict(json.loads(string))

def load(infile):
    ''' Load the given JSON file.

    Params:
        infile The input filename or file object

    Returns:
        The models

    '''
    # If the input is a string then open and read from that file
    if isinstance(infile, str):
        with open(infile, 'r') as infile:
            return loads(infile.read())

    # Otherwise assume the input is a file and read from it
    else:
        return loads(infile.read())
