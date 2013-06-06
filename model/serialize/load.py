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
    from dials.model.serialize.imageset import imageset_from_dict
    from dials.model.serialize.beam import beam_from_dict
    from dials.model.serialize.detector import detector_from_dict
    from dials.model.serialize.goniometer import goniometer_from_dict
    from dials.model.serialize.scan import scan_from_dict
    from dials.model.serialize.crystal import crystal_from_dict

    # Create the named tuple type
    Model = namedtuple('Model', ['imageset', 'beam', 'detector',
                                 'goniometer', 'scan', 'crystal'])

    # Convert all dictionaries to models
    if 'imageset' in obj:
        imageset = imageset_from_dict(obj['imageset'])
    else:
        imageset = None

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
    return Model(imageset=imageset,
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


def load2(infile):

    from collections import namedtuple
    from dxtbx.imageset import ImageSetFactory

    # Get all the information from the model file
    imageset, beam, detector, goniometer, scan, crystal = load(infile)

    # If the filenames parameter has been set
    if imageset != None:

        # If the beam is not set, then get it from the sweep
        if beam == None:
            beam = imageset.get_beam()

        # If the detector is not set, then get it from the sweep
        if detector == None:
            detector = imageset.get_detector()

        # If the goniometer is not set, then get it from the sweep
        if goniometer == None:
            goniometer = imageset.get_goniometer()

        # If the scan is not set, then get it from the sweep
        if scan == None:
            scan = imageset.get_scan()

    # Create a named tuple with the results
    Model = namedtuple('Model', ['sweep', 'beam', 'detector',
                                 'goniometer', 'scan', 'crystal'])

    # Return the tuple with the results
    return Model(sweep=sweep, beam=beam, detector=detector,
                 goniometer=gonimeter, scan=scan, crystal=crystal)
