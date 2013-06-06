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

def basic_imageset_to_dict(imageset):
    ''' Convert an imageset to a dictionary

    Params:
        imageset The imageset

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    from dials.model.serialize.beam import beam_to_dict
    from dials.model.serialize.detector import detector_to_dict

    # Return the dictionary representation
    return OrderedDict([
        ("filenames", imageset.paths()),
        ("beam", beam_to_dict(imageset.get_beam())),
        ("detector", detector_to_dict(imageset.get_detector()))])

def imagesweep_to_dict(sweep):
    ''' Convert a sweep to a dictionary

    Params:
        sweep The sweep

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    from dials.model.serialize.beam import beam_to_dict
    from dials.model.serialize.detector import detector_to_dict
    from dials.model.serialize.goniometer import goniometer_to_dict
    from dials.model.serialize.scan import scan_to_dict

    # Return the dictionary representation
    return OrderedDict([
        ("template", template_format_to_string(sweep.get_template())),
        ("beam", beam_to_dict(sweep.get_beam())),
        ("detector", detector_to_dict(sweep.get_detector())),
        ("goniometer", goniometer_to_dict(sweep.get_goniometer())),
        ("scan", scan_to_dict(sweep.get_scan()))])

def replace_template_format_with_hash(match):
    return '#'*len(match.group(0))

def template_format_to_string(template):
    import re
    return re.sub(r'%0[0-9]+d', replace_template_format_with_hash, template)


def imageset_to_dict(imageset):
    ''' Convert the imageset to a dictionary

    Params:
        imageset The imageset

    Returns:
        A dictionary of the parameters

    '''
    from dxtbx.imageset import ImageSet, ImageSweep

    # If this is an imageset then return a list of filenames
    if isinstance(imageset, ImageSweep):
        return imagesweep_to_dict(imageset)
    elif isinstance(imageset, ImageSet):
        return basic_imageset_to_dict(imageset)
    else:
        raise TypeError("Unknown ImageSet Type")

def template_string_to_glob_expr(template):
    '''Convert the template to a glob expression.'''
    pfx = template.split('#')[0]
    sfx = template.split('#')[-1]
    return '%s%s%s' % (pfx, '[0-9]'*template.count('#'), sfx)

def template_string_number_index(template):
    '''Get the number idex of the template.'''
    pfx = template.split('#')[0]
    sfx = template.split('#')[-1]
    return len(pfx), len(pfx) + template.count('#')

def locate_files_matching_template_string(template):
    '''Return all files matching template.'''
    from glob import glob
    return glob(template_string_to_glob_expr(template))

def template_image_range(template):
    '''Return the image range of files with this template.'''

    # Find the files matching the template
    filenames = locate_files_matching_template_string(template)
    filenames = sorted(filenames)

    # Check that the template matches some files
    if len(filenames) == 0:
        raise ValueError('Template doesn\'t match any files.')

    # Get the templete format
    index = slice(*template_string_number_index(template))

    # Get the first and last indices
    first = int(filenames[0][index])
    last  = int(filenames[-1][index])

    # Reutrn the image range
    return (first, last)

def basic_imageset_from_dict(d):
    ''' Construct an ImageSet class from the dictionary.'''
    from dxtbx.imageset import ImageSetFactory
    from dials.model.serialize.beam import beam_from_dict
    from dials.model.serialize.detector import detector_from_dict

    # Get the filename list and create the imageset
    filenames = map(str, d['filenames'])
    imageset = ImageSetFactory.new(filenames)[0]

    # Set models
    imageset.set_beam(beam_from_dict(d.get('beam')))
    imageset.set_detector(detector_from_dict(d.get('detector')))

    # Return the imageset
    return imageset

def imagesweep_from_dict(d):
    '''Construct and image sweep from the dictionary.'''
    from dxtbx.imageset import ImageSetFactory
    from dials.model.serialize.beam import beam_from_dict
    from dials.model.serialize.detector import detector_from_dict
    from dials.model.serialize.goniometer import goniometer_from_dict
    from dials.model.serialize.scan import scan_from_dict

    # Get the template (required)
    template = str(d['template'])

    # Get the scan
    scan = scan_from_dict(d.get('scan'))

    # If the scan isn't set, find all available files
    if scan is None:
        image_range = template_image_range(template)
    else:
        image_range = scan.get_image_range()

    # Construct the sweep
    sweep = ImageSetFactory.from_template(template, image_range)[0]
    sweep.set_beam(beam_from_dict(d.get('beam')))
    sweep.set_goniometer(goniometer_from_dict(d.get('goniometer')))
    sweep.set_detector(detector_from_dict(d.get('detector')))

    # Return the sweep
    return sweep

def imageset_from_dict(d):
    ''' Convert the dictionary to a sweep

    Params:
        d The dictionary of parameters

    Returns:
        The sweep

    '''
    if d == None:
        return None
    if "filenames" in d:
        return basic_imageset_from_dict(d)
    elif "template" in d:
        return imagesweep_from_dict(d)
    else:
        raise TypeError("Unable to deserialize given imageset")
