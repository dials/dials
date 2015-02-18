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

  def __init__(self, params, experiments):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param experiments: The experiment list

    '''
    self.experiments = experiments

  def compute_centroid(self, reflections):
    '''
    Compute the centroid.

    :param reflections: The list of reflections

    '''
    from logging import info
    from dials.array_family import flex
    from time import time
    from operator import mul

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

    # Convert to mm for each reflection
    value = flex.vec3_double(len(centroid))
    variance = flex.vec3_double(len(centroid))
    for i in range(len(reflections)):
      exp_id = reflections[i]['id']
      pan_id = reflections[i]['panel']
      pos = reflections[i]['xyzobs.px.value']
      var = reflections[i]['xyzobs.px.variance']
      exp = self.experiments[exp_id]
      panel = exp.detector[pan_id]
      pixel_size = panel.get_pixel_size()
      if exp.scan is None:
        oscillation = (0,0)
      else:
        oscillation = exp.scan.get_oscillation(deg=False)
      scale = pixel_size + (oscillation[1],)
      scale2 = map(mul, scale, scale)
      x, y, z = pos
      xy_mm = panel.pixel_to_millimeter((x, y))
      if exp.scan is None:
        z_rad = 0
      else:
        z_rad = exp.scan.get_angle_from_array_index(z, deg=False)

      # Set the position, variance and squared width in mm/rad
      # N.B assuming locally flat pixel to millimeter transform
      # for variance calculation.
      pos_mm = xy_mm + (z_rad,)
      var_mm = map(mul, var, scale2)
      value[i] = pos_mm
      variance[i] = var_mm
    reflections['xyzobs.mm.value'] = value
    reflections['xyzobs.mm.variance'] = variance
    info('  successfully processed %d reflections' % len(reflections))
    info('  time taken: %g seconds' % (time() - start_time))
