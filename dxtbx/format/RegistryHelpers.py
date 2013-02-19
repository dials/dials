#!/usr/bin/env python
# RegistryHelpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Things to help the ImageFormat registry to work.

from __future__ import division

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

def LoadFormatClasses():
    '''Look for files named Format(something).py in the sensible
    places (i.e. in the xia2 distribution and in the users home area)
    and import the corresponding modules using their fully qualified
    names.'''

    import dxtbx.format

    # FIXME in here - dxtbx should already be in os.path - look for it there,
    # also wouldn't it be tidy to refer to a Phil parameter?

    format_dir = os.path.split(dxtbx.format.__file__)[0]

    if os.name == 'nt':
        home = os.path.join(os.environ['HOMEDRIVE'],
                            os.environ['HOMEPATH'])
    else:
        home = os.environ['HOME']

    for f in os.listdir(format_dir):
        if 'Format' in f[:6] and '.py' in f[-3:]:
            name = f[:-3]
            fqname = dxtbx.format.__name__ + '.' + name
            _LoadFormatModule(name, fqname, format_dir)

    format_dir = os.path.join(home, '.dxtbx')
    if os.path.exists(format_dir):
        for f in os.listdir(format_dir):
            if 'Format' in f[:6] and '.py' in f[:-3]:
                name = f[:-3]
                _LoadFormatModule(name, name, format_dir)

def _LoadFormatModule(name, fqname, path):
    '''Load a format class module, which will trigger the automated self
    registration. This module will not therefore need to publish anything
    as the module will self publish. The idea being that these format classes
    were found by the search procedure above.'''

    # Early return if module already imported.
    try:
        sys.modules[fqname]
        return
    except KeyError:
        pass

    module, pathname, description = imp.find_module(name, [path])

    try:
        imp.load_module(fqname, module, pathname, description)
    except exceptions.Exception, e:
        traceback.print_exc(sys.stderr)
    finally:
        module.close()

    return

if __name__ == '__main__':

    LoadFormatClasses()
