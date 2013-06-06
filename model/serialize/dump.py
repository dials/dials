#!/usr/bin/env python
#
# dials.model.serialize.dump.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def model_to_dict(filenames=None, beam=None, detector=None,
                  goniometer=None, scan=None, crystal=None):
    ''' Dump the experimental models to a dictionary

    Params:
        filenames The filename list
        beam The beam model
        detector The detector model
        gonio The goniometer model
        scan The scan model
        crystal The crystal model

    Returns:
        The JSON string

    '''
    from collections import OrderedDict
    from dials.model.serialize.beam import beam_to_dict
    from dials.model.serialize.detector import detector_to_dict
    from dials.model.serialize.goniometer import goniometer_to_dict
    from dials.model.serialize.scan import scan_to_dict
    from dials.model.serialize.crystal import crystal_to_dict

    # Convert all models to dictionaries
    model = OrderedDict()
    if filenames != None: model['filenames'] = filenames
    if beam != None: model['beam'] = beam_to_dict(beam)
    if detector != None: model['detector'] = detector_to_dict(detector)
    if goniometer != None: model['goniometer'] = goniometer_to_dict(goniometer)
    if scan != None: model['scan'] = scan_to_dict(scan)
    if crystal != None: model['crystal'] = crystal_to_dict(crystal)

    # Return the model
    return model

def dumps(**kwargs):
    ''' Dump the given object to string.

    Input must be of the form

        dumps((beam, detector, goniometer, scan))

    or
        dumps(crystal)

    Params:
        filenames The filename list
        beam The beam model
        detector The detector model
        gonio The goniometer model
        scan The scan model
        crystal The crystal model

    Returns:
        The JSON string

    '''
    import json
    import textwrap

    # Return as a JSON string
    return '\n'.join(textwrap.wrap(json.dumps(model_to_dict(**kwargs)), 80))

def dump(outfile, **kwargs):
    ''' Dump the given object to file.

    Input must be of the form

        dumps((beam, detector, goniometer, scan))

    or
        dumps(crystal)

    Params:
        outfile The output file name or file object
        filenames The filename list
        beam The beam model
        detector The detector model
        gonio The goniometer model
        scan The scan model
        crystal The crystal model

    '''
    # If the input is a string then open and write to that file
    if isinstance(outfile, str):
        with open(outfile, 'w') as outfile:
            outfile.write(dumps(**kwargs))

    # Otherwise assume the input is a file and write to it
    else:
        outfile.write(dumps(**kwargs))
