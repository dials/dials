#
# centroid_px_to_mm.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

def centroid_px_to_mm(detector, scan, position, variance, sd_error):
  '''Convenience function to calculate centroid in mm/rad from px'''
  from operator import mul

  # Get the pixel to millimeter function
  assert(len(detector) == 1)
  return centroid_px_to_mm_panel(
      detector[0], scan, position, variance, sd_error)

def centroid_px_to_mm_panel(panel, scan, position, variance, sd_error):
  '''Convenience function to calculate centroid in mm/rad from px'''
  from operator import mul

  # Get the pixel to millimeter function
  pixel_size = panel.get_pixel_size()
  if scan is None:
    oscillation = (0,0)
  else:
    oscillation = scan.get_oscillation(deg=False)
  scale = pixel_size + (oscillation[1],)
  scale2 = map(mul, scale, scale)

  # Convert Pixel coordinate into mm/rad
  x, y, z = position
  xy_mm = panel.pixel_to_millimeter((x, y))

  if scan is None:
    z_rad = 0
  else:
    z_rad = scan.get_angle_from_array_index(z, deg=False)

  # Set the position, variance and squared width in mm/rad
  # N.B assuming locally flat pixel to millimeter transform
  # for variance calculation.
  position_mm = xy_mm + (z_rad,)
  variance_mm = map(mul, variance, scale2)
  sd_error_mm = map(mul, sd_error, scale2)

  # Return the stuff in mm/rad
  return position_mm, variance_mm, sd_error_mm
