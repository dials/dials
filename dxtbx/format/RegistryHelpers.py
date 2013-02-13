#!/usr/bin/env python
# RegistryHelpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Things to help the ImageFormat registry to work.

import os
import sys
import imp
import exceptions
import traceback

def InheritsFromFormat(PutativeFormatClass):
    '''Check that the PutativeFormatClass inherits on some level from a class
    named Format. This is aimed at making sure it should play nice with the
    ImageFormat registry etc.'''

    if PutativeFormatClass.__bases__:
        for base in PutativeFormatClass.__bases__:
            if InheritsFromFormat(base):
                return True

    if PutativeFormatClass.__name__ == 'Format':
        return True

    return False

def LookForFormatClasses():
    '''Look for files named Format(something).py in the sensible places (i.e.
    in the xia2 distribution and in the users home area) and return a list of
    paths. N.B. the class names themselves must be unique (otherwise there
    is no hope of importing them!)'''

    import dxtbx

    format_classes = []
    file_names = []

    # FIXME in here - dxtbx should already be in os.path - look for it there,
    # also wouldn't it be tidy to refer to a Phil parameter?

    format_dir = os.path.join(os.path.split(dxtbx.__file__)[0], 'format')

    if os.name == 'nt':
        home = os.path.join(os.environ['HOMEDRIVE'],
                            os.environ['HOMEPATH'])
    else:
        home = os.environ['HOME']

    for f in os.listdir(format_dir):
        if 'Format' in f[:6] and '.py' in f[-3:]:
            assert(not f in file_names)
            file_names.append(f)
            format_classes.append(os.path.join(format_dir, f))

    return format_classes

def LoadFormatClass(FormatClass):
    '''Load a format class module, which will trigger the automated self
    registration. This module will not therefore need to publish anything
    as the module will self publish. The idea being that these format classes
    were found by the search procedure above.'''

    format_class_name = os.path.split(FormatClass)[-1][:-3]
    format_class_path = os.path.split(FormatClass)[0]

    module, path, description = imp.find_module(format_class_name,
                                                [format_class_path])

    try:
        imp.load_module(format_class_name, module, path, description)
    except exceptions.Exception, e:
        traceback.print_exc(sys.stderr)
    finally:
        module.close()

    return

if __name__ == '__main__':

    for f in LookForFormatClasses():
        LoadFormatClass(f)
