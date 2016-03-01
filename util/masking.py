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
    image = imageset[0]
    if not isinstance(image, tuple):
      image = (image,)
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
      resolution_mask_generator = ResolutionMaskGenerator(beam, panel)

      # Generate high and low resolution masks
      if self.params.d_min is not None:
        info("Generating high resolution mask:")
        info(" d_min = %f" % self.params.d_min)
        resolution_mask_generator.apply(mask, 0, self.params.d_min)
      if self.params.d_max is not None:
        info("Generating low resolution mask:")
        info(" d_max = %f" % self.params.d_max)
        d_min = self.params.d_max
        d_max = max(d_min + 1, 1e9)
        resolution_mask_generator.apply(mask, d_min, d_max)

      # Mask out the resolution range
      for drange in self.params.resolution_range:
        d_min=min(drange)
        d_max=max(drange)
        assert d_min < d_max, "d_min must be < d_max"
        info("Generating resolution range mask:")
        info(" d_min = %f" % d_min)
        info(" d_max = %f" % d_max)
        resolution_mask_generator.apply(mask, d_min, d_max)

      # Mask out the resolution ranges for the ice rings
      for drange in generate_ice_ring_resolution_ranges(
          beam, panel, self.params.ice_rings):
        d_min=min(drange)
        d_max=max(drange)
        assert d_min < d_max, "d_min must be < d_max"
        info("Generating ice ring mask:")
        info(" d_min = %f" % d_min)
        info(" d_max = %f" % d_max)
        resolution_mask_generator.apply(mask, d_min, d_max)

      # Add to the list
      masks.append(mask)

    # Return the mask
    return tuple(masks)
