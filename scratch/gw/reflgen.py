from __future__ import division

def reflgen(num_refl):

  from dials.model.data import ReflectionList, Reflection
  from dials.algorithms import shoebox
  import random
  import math

  rlist = ReflectionList(num_refl)

  # FIXME in here generate random cell and random U matrix

  # FIXME decide shape of reflection in reciprocal space for transformation
  # reasons - will need this to assign shoebox size and also generate the
  # basic density

  for j in range(num_refl):
    rlist[j].miller_index = (random.randint(0, 20),
                             random.randint(0, 20),
                             random.randint(0, 20))
    rlist[j].rotation_angle = 2 * math.pi * random.random()

    # FIXME write some prediction code or use some to derive these things
    #
    # P* = RUBh; s1 (beam_vector) = s0 + p*

    rlist[j].beam_vector = (0, 0, 0)
    rlist[j].image_coord_px = (0, 0)
    rlist[j].image_coord_mm = (0, 0)
    rlist[j].frame_number = 0
    rlist[j].panel_number = 0
    rlist[j].bounding_box = (0, 1, 0, 1, 0, 1)
    rlist[j].centroid_position = (0, 0, 0)
    rlist[j].centroid_variance = (0, 0, 0)
    rlist[j].centroid_sq_width = (0, 0, 0)
    rlist[j].intensity = 0
    rlist[j].intensity_variance = 0
    rlist[j].corrected_intensity = 0
    rlist[j].corrected_intensity_variance = 0

  shoebox.allocate(rlist)

  for j in range(num_refl):

    # FIXME in here fake up background and peak and mask but don't assign
    # the "true" values here unless we want to test things? Maybe that
    # would be best ... make sure we can reproduce the observations?

    shoebox = rlist[j].shoebox
    mask = rlist[j].shoebox_mask
    background = rlist[j].shoebox_background
