#!/usr/bin/env python
#
# background.py
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

background
  .help = "Configure the background algorithm."
{
  algorithm = null *xds fable flat inclined curved
    .help = "The background algorithm."
    .type = choice

  xds
    .help = "Parameters for xds background subtraction."
  {
    min_pixels = 10
      .help = "The minimum number of pixels to use in calculating the"
              "background intensity."
      .type = int

    n_sigma = 3.0
      .help = "The number of standard deviations above the mean pixel intensity"
              "below which the pixel will be classified as background."
      .type = float
  }

  fable
    .help = "Parameters for fable background subtraction."
  {
    min_pixels = 10
      .help = "The minimum number of pixels to use in calculating the"
              "background intensity."
      .type = int

    n_sigma = 3.0
      .help = "The number of standard deviations above the mean pixel intensity"
              "below which the pixel will be classified as background."
      .type = float
  }
}
''')
