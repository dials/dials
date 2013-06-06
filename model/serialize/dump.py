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

def compact_simple_list(match):
    ''' Callback function. Given a simple list match, compact it and ensure
    that it wraps around by 80 characters.

    Params:
        match The regular expression match

    Returns:
        The string to replace the expression with

    '''
    import textwrap

    # Calculate the initial indent as the length of the first match group
    initial_indent = len(match.group(1))

    # Get the lines in the match
    lines = match.group(2).splitlines()

    # Set the indent by finding the indent of the first lines
    if len(lines) > 1:
        subsequent_indent = len(lines[1]) - len(lines[1].lstrip())
    else:
        subsequent_indent = 0

    # Strip whitespace from the lines
    lines = [l.strip() for l in lines]

    # Create and return the string wrapped about 80 chars
    list_string = '\n'.join(textwrap.wrap(' '.join(lines),
        80, initial_indent=' '*initial_indent,
        subsequent_indent=' '*subsequent_indent)).lstrip()

    # Return the string
    return match.group(1) + list_string

def compact_simple_lists(string):
    ''' Find simple lists in the string and compact.

    Params:
        string The input JSON string

    Returns:
        The output JSON string

    '''
    import re
    return re.sub(r'(.*"\w+".*:.*)(\[[^\{\}\[\]]*\])', compact_simple_list, string)

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
    string = json.dumps(model_to_dict(**kwargs), indent=2)

    # Hack to make more readable
    return compact_simple_lists(string)

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
