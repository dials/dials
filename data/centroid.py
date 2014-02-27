#!/usr/bin/env python
#
# centroid.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from libtbx.phil import parse

phil_scope = parse('''

centroid
  .help = "Configure the centroid algorithm."
{
  algorithm = *simple
    .help = "The centroid algorithm."
    .type = choice
}

''')
