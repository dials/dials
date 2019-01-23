#!/usr/bin/env python
#
# masking.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import logging

from iotbx.phil import parse
from scitbx.array_family import flex

logger = logging.getLogger(__name__)

phil_scope = parse("""
  border = 0
    .type = int
    .help = "The border around the edge of the image."

  use_trusted_range = True
    .type = bool
    .help = "Use the trusted range to mask bad pixels."

  d_min = None
    .help = "The high resolution limit in Angstrom for a pixel to be"
            "accepted by the filtering algorithm."
    .type = float(value_min=0)

  d_max = None
    .help = "The low resolution limit in Angstrom for a pixel to be"
            "accepted by the filtering algorithm."
    .type = float(value_min=0)

  resolution_range = None
    .multiple = true
    .type = floats(2)
    .help = "an untrusted resolution range"

  untrusted
    .multiple = True
  {

    panel = 0
      .type = int
      .help = "The panel number"

    circle = None
      .type = ints(3)
      .help = "An untrusted circle (xc, yc, r)"

    rectangle = None
      .type = ints(4)
      .help = "An untrusted rectangle (x0, x1, y0, y1)"

    polygon = None
      .type = ints(value_min=0)
      .help = "The pixel coordinates (fast, slow) that define the corners "
              "of the untrusted polygon. Spots whose centroids fall within "
              "the bounds of the untrusted polygon will be rejected."

    pixel = None
      .type = ints(2, value_min=0)
      .help = "An untrusted pixel (x, y)"

  }

  ice_rings {
    filter = False
      .type = bool
    unit_cell = 4.498,4.498,7.338,90,90,120
      .type = unit_cell
      .help = "The unit cell to generate d_spacings for powder rings."
      .expert_level = 1
    space_group = 194
      .type = space_group
      .help = "The space group used to generate d_spacings for powder rings."
      .expert_level = 1
    width = 0.002
      .type = float(value_min=0.0)
      .help = "The width of an ice ring (in 1/d^2)."
      .expert_level = 1
    d_min = None
      .type = float(value_min=0.0)
      .help = "The high resolution limit (otherwise use detector d_min)"
  }

""", process_includes=True)

def generate_ice_ring_resolution_ranges(beam, panel, params):
  '''
  Generate a set of resolution ranges from the ice ring parameters

  '''
  from cctbx import crystal
  from math import sqrt

  if params.filter is True:

    # Get the crystal symmetry
    crystal_symmetry = crystal.symmetry(
      unit_cell=params.unit_cell,
      space_group=params.space_group.group())

    # Get the half width
    half_width = params.width * 0.5

    # Set the high resolution
    if params.d_min is None:
      d_min = panel.get_max_resolution_at_corners(beam.get_s0())
    else:
      d_min = params.d_min

    # Build the miller set
    ms = crystal_symmetry.build_miller_set(
      anomalous_flag=False,
      d_min=d_min)
    ms = ms.sort(by_value="resolution")

    # Yield all the d ranges
    for j, d in enumerate(ms.d_spacings().data()):
      d_sq_inv = 1.0 / (d**2)
      d_sq_inv_min = d_sq_inv - half_width
      d_sq_inv_max = d_sq_inv + half_width
      d_min = sqrt(1.0 / d_sq_inv_min)
      d_max = sqrt(1.0 / d_sq_inv_max)
      yield (d_min, d_max)


class MaskGenerator(object):
  ''' Generate a mask. '''

  def __init__(self, params):
    ''' Set the parameters. '''
    self.params = params

  def generate(self, imageset):
    ''' Generate the mask. '''
    from dials.util.ext import ResolutionMaskGenerator
    from dials.util.ext import mask_untrusted_rectangle
    from dials.util.ext import mask_untrusted_circle
    from dials.util.ext import mask_untrusted_polygon
    from dials.array_family import flex

    # Get the detector and beam
    detector = imageset.get_detector()
    beam = imageset.get_beam()

    # Get the first image
    image = imageset.get_raw_data(0)
    assert(len(detector) == len(image))

    # Create the mask for each image
    masks = []
    for index, (im, panel) in enumerate(zip(image, detector)):

      # The image width height
      height, width = im.all()

      # Create the basic mask from the trusted range
      if self.params.use_trusted_range:
        low, high = panel.get_trusted_range()
        imd = im.as_double()
        mask = (imd > low) & (imd < high)
      else:
        mask = flex.bool(flex.grid(im.all()), True)

      # Add a border around the image
      if self.params.border > 0:
        logger.info("Generating border mask:")
        logger.info(" border = %d" % self.params.border)
        border = self.params.border
        height, width = mask.all()
        borderx = flex.bool(flex.grid(border, width), False)
        bordery = flex.bool(flex.grid(height, border), False)
        mask[0:border,:] = borderx
        mask[-border:,:] = borderx
        mask[:,0:border] = bordery
        mask[:,-border:] = bordery

      # Apply the untrusted regions
      for region in self.params.untrusted:
        if region.panel == index:
          if region.circle is not None:
            xc, yc, radius = region.circle
            logger.info("Generating circle mask:")
            logger.info(" panel = %d" % region.panel)
            logger.info(" xc = %d" % xc)
            logger.info(" yc = %d" % yc)
            logger.info(" radius = %d" % radius)
            mask_untrusted_circle(mask, xc, yc, radius)
          if region.rectangle is not None:
            x0, x1, y0, y1 = region.rectangle
            logger.info("Generating rectangle mask:")
            logger.info(" panel = %d" % region.panel)
            logger.info(" x0 = %d" % x0)
            logger.info(" y0 = %d" % y0)
            logger.info(" x1 = %d" % x1)
            logger.info(" y1 = %d" % y1)
            mask_untrusted_rectangle(mask, x0, x1, y0, y1)
          if region.polygon is not None:
            assert len(region.polygon) % 2 == 0, "Polygon must contain 2D coords"
            vertices = []
            for i in range(int(len(region.polygon)/2)):
              x = region.polygon[2*i]
              y = region.polygon[2*i+1]
              vertices.append((x,y))
            polygon = flex.vec2_double(vertices)
            logger.info("Generating polygon mask:")
            logger.info(" panel = %d" % region.panel)
            for vertex in vertices:
              logger.info(" coord = (%d, %d)" % (vertex))
            mask_untrusted_polygon(mask, polygon)
          if region.pixel is not None:
            mask[region.pixel] = False

      # Create the resolution mask generator
      class ResolutionMaskGeneratorGetter(object):
        def __init__(self, beam, panel):
          self.beam = beam
          self.panel = panel
          self.result = None
        def __call__(self):
          if self.result is None:
            self.result = ResolutionMaskGenerator(beam, panel)
          return self.result
      get_resolution_mask_generator = ResolutionMaskGeneratorGetter(beam, panel)

      # Generate high and low resolution masks
      if self.params.d_min is not None:
        logger.info("Generating high resolution mask:")
        logger.info(" d_min = %f" % self.params.d_min)
        get_resolution_mask_generator().apply(mask, 0, self.params.d_min)
      if self.params.d_max is not None:
        logger.info("Generating low resolution mask:")
        logger.info(" d_max = %f" % self.params.d_max)
        d_min = self.params.d_max
        d_max = max(d_min + 1, 1e9)
        get_resolution_mask_generator().apply(mask, d_min, d_max)

      # Mask out the resolution range
      for drange in self.params.resolution_range:
        d_min=min(drange)
        d_max=max(drange)
        assert d_min < d_max, "d_min must be < d_max"
        logger.info("Generating resolution range mask:")
        logger.info(" d_min = %f" % d_min)
        logger.info(" d_max = %f" % d_max)
        get_resolution_mask_generator().apply(mask, d_min, d_max)

      # Mask out the resolution ranges for the ice rings
      for drange in generate_ice_ring_resolution_ranges(
          beam, panel, self.params.ice_rings):
        d_min=min(drange)
        d_max=max(drange)
        assert d_min < d_max, "d_min must be < d_max"
        logger.info("Generating ice ring mask:")
        logger.info(" d_min = %f" % d_min)
        logger.info(" d_max = %f" % d_max)
        get_resolution_mask_generator().apply(mask, d_min, d_max)

      # Add to the list
      masks.append(mask)

    # Return the mask
    return tuple(masks)


class GoniometerShadowMaskGenerator(object):

  def __init__(self, goniometer, extrema_at_datum, axis):
    self.goniometer = goniometer
    self._extrema_at_datum = extrema_at_datum
    self.axis = axis

  def extrema_at_scan_angle(self, scan_angle):
    from scitbx import matrix

    axes = self.goniometer.get_axes()
    angles = self.goniometer.get_angles()
    scan_axis = self.goniometer.get_scan_axis()
    angles[scan_axis] = scan_angle
    extrema = self._extrema_at_datum.deep_copy()

    for i in range(len(axes)):
      sel = (self.axis <= i)
      rotation = matrix.col(
        axes[i]).axis_and_angle_as_r3_rotation_matrix(angles[i], deg=True)
      extrema.set_selected(sel, rotation.elems * extrema.select(sel))

    return extrema

  def project_extrema(self, detector, scan_angle):
    from dials.util.ext import is_inside_polygon
    coords = self.extrema_at_scan_angle(scan_angle)
    shadow_boundary = []

    for p_id, p in enumerate(detector):
      # project coordinates onto panel plane
      a = p.get_D_matrix() * coords
      x, y, z = a.parts()
      valid = z > 0
      x.set_selected(valid, x.select(valid)/z.select(valid))
      y.set_selected(valid, y.select(valid)/z.select(valid))

      if valid.count(True) < 3:
        # no shadow projected onto this panel
        shadow_boundary.append(flex.vec2_double())
        continue

      # Compute convex hull of shadow points
      points = flex.vec2_double(x.select(valid), y.select(valid))
      shadow = flex.vec2_double(convex_hull(points))
      shadow *= 1/p.get_pixel_size()[0]

      shadow_orig = shadow.deep_copy()

      for i in (0, p.get_image_size()[0]):
        points = flex.vec2_double(flex.double(p.get_image_size()[1], i),
                                  flex.double_range(0, p.get_image_size()[1]))
        inside = is_inside_polygon(shadow_orig, points)
        # only add those points needed to define vertices of shadow
        inside_isel = inside.iselection()
        outside_isel = (~inside).iselection()
        while inside_isel.size():
          j = inside_isel[0]
          shadow.append(points[j])
          outside_isel = outside_isel.select(outside_isel > j)
          if outside_isel.size() == 0:
            shadow.append(points[inside_isel[-1]])
            break
          sel = inside_isel >= outside_isel[0]
          if sel.count(True) == 0:
            shadow.append(points[inside_isel[-1]])
            break
          inside_isel = inside_isel.select(sel)

      for i in (0, p.get_image_size()[1]):
        points = flex.vec2_double(flex.double_range(0, p.get_image_size()[0]),
                                  flex.double(p.get_image_size()[0], i))
        inside = is_inside_polygon(shadow_orig, points)
        # only add those points needed to define vertices of shadow
        inside_isel = inside.iselection()
        outside_isel = (~inside).iselection()
        while inside_isel.size():
          j = inside_isel[0]
          shadow.append(points[j])
          outside_isel = outside_isel.select(outside_isel > j)
          if outside_isel.size() == 0:
            shadow.append(points[inside_isel[-1]])
            break
          sel = inside_isel >= outside_isel[0]
          if sel.count(True) == 0:
            shadow.append(points[inside_isel[-1]])
            break
          inside_isel = inside_isel.select(sel)

      # Select only those vertices that are within the panel dimensions
      n_px = p.get_image_size()
      x, y = shadow.parts()
      valid = (x >= 0) & (x <= n_px[0]) & (y >= 0) & (y <= n_px[1])
      shadow = shadow.select(valid)

      # sort vertices clockwise from centre of mass
      from scitbx.math import principal_axes_of_inertia_2d
      com = principal_axes_of_inertia_2d(shadow).center_of_mass()
      sx, sy = shadow.parts()
      shadow = shadow.select(
        flex.sort_permutation(flex.atan2(sy-com[1],sx-com[0])))

      shadow_boundary.append(shadow)

    return shadow_boundary

  def get_mask(self, detector, scan_angle):
    shadow_boundary = self.project_extrema(detector, scan_angle)
    from dials.util.ext import mask_untrusted_polygon
    mask = []
    for panel_id in range(len(detector)):
      m = None
      if shadow_boundary[panel_id].size() > 3:
        m = flex.bool(
          flex.grid(reversed(detector[panel_id].get_image_size())), True)
        mask_untrusted_polygon(m, shadow_boundary[panel_id])
      mask.append(m)
    return mask


#https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python
#https://github.com/thepracticaldev/orly-full-res/blob/master/copyingandpasting-big.png
def convex_hull(points):
  """Computes the convex hull of a set of 2D points.

  Input: an iterable sequence of (x, y) pairs representing the points.
  Output: a list of vertices of the convex hull in counter-clockwise order,
    starting from the vertex with the lexicographically smallest coordinates.
  Implements Andrew's monotone chain algorithm. O(n log n) complexity.
  """

  # Sort the points lexicographically (tuples are compared lexicographically).
  # Remove duplicates to detect the case we have just one unique point.
  points = sorted(set(points))

  # Boring case: no points or a single point, possibly repeated multiple times.
  if len(points) <= 1:
    return points

  # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
  # Returns a positive value, if OAB makes a counter-clockwise turn,
  # negative for clockwise turn, and zero if the points are collinear.
  def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

  # Build lower hull
  lower = []
  for p in points:
    while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
      lower.pop()
    lower.append(p)

  # Build upper hull
  upper = []
  for p in reversed(points):
    while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
      upper.pop()
    upper.append(p)

  # Concatenation of the lower and upper hulls gives the convex hull.
  # Last point of each list is omitted because it is repeated at the beginning of the other list.
  return lower[:-1] + upper[:-1]
