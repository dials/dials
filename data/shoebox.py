#!/usr/bin/env python
#
# shoebox.py
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

shoebox
  .help = "Parameters used in the integration shoebox"
{
  n_blocks = 1
    .help = "The number of blocks to integrate in."
    .type = int

  n_sigma = 3
    .help = "The number of standard deviations of the beam divergence and the"
            "mosaicity to use for the bounding box size."
    .type = float

  sigma_b = None
    .help = "The E.S.D. of the beam divergence"
    .type = float

  sigma_m = None
    .help = "The E.S.D. of the reflecting range"
    .type = float
}

''')
