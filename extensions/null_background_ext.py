#!/usr/bin/env python
#
# null_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface


class NullBackgroundExt(BackgroundIface):
  ''' An extension class implementing Null background subtraction. '''

  name = 'null'

  def compute_background(self, reflections, image_volume=None):
    '''
    Compute the background.

    :param reflections: The list of reflections

    '''
    from dials.algorithms.background import set_shoebox_background_value
    from dials.array_family import flex
    set_shoebox_background_value(reflections['shoebox'], 0)
    return flex.bool(len(reflections), True)
