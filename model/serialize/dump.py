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

def dumps(obj, compact=False):
    ''' Dump the given object to string.

    Params:
        obj The crystal model
        compact Write in compact representation

    Returns:
        The JSON string

    '''
    import json
    import textwrap
    from dials.model.serialize import crystal
    from dxtbx.serialize.dump import compact_simple_lists

    # Return as a JSON string
    if compact == False:
        string = json.dumps(crystal.to_dict(obj), indent=2)

        # Hack to make more readable
        string = compact_simple_lists(string)

    else:
        string = json.dumps(crystal.to_dict(obj), separators=(',',':'))

    # Return the string
    return string

def dump(obj, outfile, compact=False):
    ''' Dump the given object to file.

    Params:
        obj The crystal to dump
        outfile The output file name or file object
        compact Write in compact representation

    '''
    # If the input is a string then open and write to that file
    if isinstance(outfile, str):
        with open(outfile, 'w') as outfile:
            outfile.write(dumps_crystal(obj, compact))

    # Otherwise assume the input is a file and write to it
    else:
        outfile.write(dumps_crystal(obj, compact))
