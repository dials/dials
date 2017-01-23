#!/usr/bin/env python
#
# lookup.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division
from libtbx.phil import parse

phil_scope = parse('''

lookup
  .help = "Parameters specifying lookup file path"
{
  mask = None
    .help = "The path to the mask file."
    .type = str
}

''')
