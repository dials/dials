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

  default = True

  def compute_centroid(self, reflections):
    ''' Compute the centroid. '''
    from logging import info
    from dials.array_family import flex
    from time import time

    # Compute the reflections
    start_time = time()
    info('')
    info(' Beginning centroid calculation')
    info('  using %d reflections' % len(reflections))
    centroid = reflections['shoebox'].centroid_valid()
    value = flex.vec3_double(len(centroid))
    variance = flex.vec3_double(len(centroid))
    for i in range(len(centroid)):
      value[i] = centroid[i].px.position
      variance[i] = centroid[i].px.variance
    reflections['xyzobs.px.value'] = value
    reflections['xyzobs.px.variance'] = variance
    info('  successfully processed %d reflections' % len(reflections))
    info('  time taken: %g seconds' % (time() - start_time))
