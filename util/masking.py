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

from __future__ import division
from iotbx.phil import parse
from scitbx.array_family import flex

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
    width = 0.06
      .type = float(value_min=0.0)
      .help = "The width of an ice ring (in d-spacing)."
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
      d_min = d - half_width
      d_max = d + half_width
      yield (d_min, d_max)


class MaskGenerator(object):
  ''' Generate a mask. '''

  def __init__(self, params):
    ''' Set the parameters. '''
    self.params = params

  def generate(self, imageset):
    ''' Generate the mask. '''
    from dials.util import ResolutionMaskGenerator
    from dials.util import mask_untrusted_rectangle
    from dials.util import mask_untrusted_circle
    from dials.util import mask_untrusted_polygon
    from dials.util import mask_untrusted_resolution_range
    from dials.array_family import flex
    from math import floor, ceil
    from logging import info

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
        info("Generating border mask:")
        info(" border = %d" % self.params.border)
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
            info("Generating circle mask:")
            info(" panel = %d" % region.panel)
            info(" xc = %d" % xc)
            info(" yc = %d" % yc)
            info(" radius = %d" % radius)
            mask_untrusted_circle(mask, xc, yc, radius)
          if region.rectangle is not None:
            x0, x1, y0, y1 = region.rectangle
            info("Generating rectangle mask:")
            info(" panel = %d" % region.panel)
            info(" x0 = %d" % x0)
            info(" y0 = %d" % y0)
            info(" x1 = %d" % x1)
            info(" y1 = %d" % y1)
            mask_untrusted_rectangle(mask, x0, x1, y0, y1)
          if region.polygon is not None:
            assert len(region.polygon) % 2 == 0, "Polygon must contain 2D coords"
            vertices = []
            for i in range(int(len(region.polygon)/2)):
              x = region.polygon[2*i]
              y = region.polygon[2*i+1]
              vertices.append((x,y))
            polygon = flex.vec2_double(vertices)
            info("Generating polygon mask:")
            info(" panel = %d" % region.panel)
            for vertex in vertices:
              info(" coord = (%d, %d)" % (vertex))
            mask_untrusted_polygon(mask, polygon)

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
        info("Generating high resolution mask:")
        info(" d_min = %f" % self.params.d_min)
        get_resolution_mask_generator().apply(mask, 0, self.params.d_min)
      if self.params.d_max is not None:
        info("Generating low resolution mask:")
        info(" d_max = %f" % self.params.d_max)
        d_min = self.params.d_max
        d_max = max(d_min + 1, 1e9)
        get_resolution_mask_generator().apply(mask, d_min, d_max)

      # Mask out the resolution range
      for drange in self.params.resolution_range:
        d_min=min(drange)
        d_max=max(drange)
        assert d_min < d_max, "d_min must be < d_max"
        info("Generating resolution range mask:")
        info(" d_min = %f" % d_min)
        info(" d_max = %f" % d_max)
        get_resolution_mask_generator().apply(mask, d_min, d_max)

      # Mask out the resolution ranges for the ice rings
      for drange in generate_ice_ring_resolution_ranges(
          beam, panel, self.params.ice_rings):
        d_min=min(drange)
        d_max=max(drange)
        assert d_min < d_max, "d_min must be < d_max"
        info("Generating ice ring mask:")
        info(" d_min = %f" % d_min)
        info(" d_max = %f" % d_max)
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
    coords = self.extrema_at_scan_angle(scan_angle)
    shadow_boundary = []

    for p_id, p in enumerate(detector):
      # project coordinates onto panel plane
      a = p.get_D_matrix() * coords
      x, y, z = a.parts()
      valid = z > 0
      x.set_selected(valid, x/z)
      y.set_selected(valid, y/z)

      if valid.count(True) == 0:
        # no shadow projected onto this panel
        shadow_boundary.append(flex.vec2_double())
        continue

      # Compute convex hull of shadow points
      from scipy.spatial import ConvexHull
      import numpy as np
      points = np.array([list(x.select(valid)), list(y.select(valid))]).transpose()
      hull = ConvexHull(points, incremental=False)
      vertices = hull.vertices
      shadow = flex.vec2_double(points[v] for v in vertices)
      shadow *= 1/p.get_pixel_size()[0]

      # Use delaunay triangulation to find out which pixels around the edge of
      # a panel are within the shadow region
      from scipy.spatial import Delaunay
      delaunay = Delaunay(np.array(list(shadow)))

      for i in (0, p.get_image_size()[0]):
        points = flex.vec2_double(flex.double(p.get_image_size()[1], i),
                                  flex.double_range(0, p.get_image_size()[1]))
        inside = flex.bool(delaunay.find_simplex(list(points)) >= 0)
        # only add those points needed to define vertices of shadow
        for j in range(len(points)):
          if j == 0 and inside[j]:
            shadow.append(points[j])
          elif inside[j] and not inside[j-1]:
            shadow.append(points[j])
          elif not inside[j] and inside[j-1]:
            shadow.append(points[j-1])

      for i in (0, p.get_image_size()[1]):
        points = flex.vec2_double(flex.double_range(1, p.get_image_size()[0]),
                                  flex.double(p.get_image_size()[0]-1, i))
        inside = flex.bool(delaunay.find_simplex(list(points)) >= 0)
        # only add those points needed to define vertices of shadow
        for j in range(len(points)):
          if j == 0 and inside[j]:
            shadow.append(points[j])
          elif inside[j] and not inside[j-1]:
            shadow.append(points[j])
          elif not inside[j] and inside[j-1]:
            shadow.append(points[j-1])

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
    from dials.util import mask_untrusted_polygon
    mask = [flex.bool(flex.grid(reversed(p.get_image_size())), True) for p in detector]
    for panel_id in range(len(detector)):
      if shadow_boundary[panel_id].size() > 3:
        mask_untrusted_polygon(mask[panel_id], shadow_boundary[panel_id])
    return mask
