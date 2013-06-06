#!/usr/bin/env python
#
# dials.model.serialize.imageset.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def imageset_to_dict(imageset):
    ''' Convert the imageset to a dictionary

    Params:
        sweep The imageset

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    from dxtbx.imageset import ImageSet, ImageSweep

    # If this is an imageset then return a list of filenames
    if isinstance(imageset, ImageSet):
        d = OrderedDict([("filenames", imageset.paths())])

    # Otherwise return a template and the image range
    elif isinstance(imageset, ImageSweep):
        template = imageset.get_template()
        array_range = imageset.get_array_range()
        image_range = (array_range[0] + 1, array_range[1])
        d = OrderedDict([
            ("template", imageset.get_template()),
            ("image_range", imageset.get_array_range())])
    else:
        raise TypeError("Unknown ImageSet Type")

    # Return the dictionary
    return d

def imageset_from_dict(d):
    ''' Convert the dictionary to a sweep

    Params:
        d The dictionary of parameters

    Returns:
        The sweep

    '''
    from dxtbx.imageset import ImageSetFactory

    # If this is a generic imageset then construct as normal
    if "filenames" in d:
        filenames = map(str, d['filenames'])
        return ImageSetFactory.new(filenames)[0]

    # If this is specifically a sweep then construct through the template
    else:
        template = str(d['template'])
        image_range = d['image_range']
        return ImageSetFactory.from_template(template, image_range)
