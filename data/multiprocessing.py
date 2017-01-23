#!/usr/bin/env python
#
# multiprocessing.py
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

mp {
  method = *multiprocessing sge lsf pbs
    .type = choice
    .help = "The multiprocessing method to use"

  nproc = 1
    .type = int(value_min=1)
    .help = "The number of processes to use."

  nthreads = 1
    .type = int(value_min=1)
    .help = "The number of local threads to use for openmp."
}

''')
