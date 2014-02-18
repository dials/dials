#!/usr/bin/env python
#
# simple_centroid_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import CentroidIface


class SimpleCentroidExt(CentroidIface):
  ''' An extension class implementing a simple centroid algorithm. '''

  name = 'simple'

  def compute_centroid(self, reflections):
    ''' Compute the centroid. '''
    from dials.util.command_line import Command

    # Compute the reflections
    Command.start('Calculating reflection centroids')
    reflections['shoebox'].centroid_valid()
    Command.end('Calculated %d reflection centroids' % len(reflections))
