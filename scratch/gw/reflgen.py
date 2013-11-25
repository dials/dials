from __future__ import division

from libtbx.phil import parse

master_phil = parse("""
nrefl = 0
  .type = int
nproc = 0
  .type = int
shoebox_size {
  x = 10
    .type = int
  y = 10
    .type = int
  z = 10
    .type = int
}
spot_size {
  x = 1.0
    .type = float
  y = 1.0
    .type = float
  z = 1.0
    .type = float
}
spot_offset {
  x = 0.0
    .type = float
  y = 0.0
    .type = float
  z = 0.0
    .type = float
}
mask_nsigma = 3.0
  .type = float
counts = 0
  .type = int
background = 0
  .type = int
pixel_mask = *all static precise
  .type = choice
background_method = *xds mosflm
  .type = choice
integration_methpd = *xds mosflm
  .type = choice
output {
  over = None
    .type = path
  under = None
    .type = path
  all = None
    .type = path
}
rotation {
  axis {
    x = 0.0
      .type = float
    y = 0.0
      .type = float
    z = 0.0
      .type = float
  }

  angle = 0.0
    .type = float
}
""")

def simple_gaussian_spots(params):
  from dials.model.data import ReflectionList
  from dials.algorithms import shoebox
  import random
  import math

  from scitbx import matrix
  r = params.rotation
  axis = matrix.col((r.axis.x, r.axis.y, r.axis.z))
  if axis.length() > 0:
    rotation = axis.axis_and_angle_as_r3_rotation_matrix(r.angle, deg=True)
  else:
    rotation = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  # generate mask and peak values

  from dials.algorithms.shoebox import MaskCode
  mask_peak = MaskCode.Valid|MaskCode.Foreground
  mask_back = MaskCode.Valid|MaskCode.Background

  rlist = ReflectionList(params.nrefl)

  from dials.util.command_line import ProgressBar
  p = ProgressBar(title = 'Generating reflections')

  for j in range(params.nrefl):
    p.update(j * 100.0 / params.nrefl)
    rlist[j].miller_index = (random.randint(0, 20),
                             random.randint(0, 20),
                             random.randint(0, 20))
    rlist[j].rotation_angle = 2 * math.pi * random.random()
    rlist[j].beam_vector = (0, 0, 0)
    rlist[j].image_coord_px = (0, 0)
    rlist[j].image_coord_mm = (0, 0)
    rlist[j].frame_number = 0
    rlist[j].panel_number = 0
    rlist[j].bounding_box = (0, params.shoebox_size.x, 0, params.shoebox_size.y,
                             0, params.shoebox_size.z)
    rlist[j].centroid_position = (0, 0, 0)
    rlist[j].centroid_variance = (0, 0, 0)
    rlist[j].centroid_sq_width = (0, 0, 0)
    rlist[j].intensity = 0
    rlist[j].intensity_variance = 0
    rlist[j].corrected_intensity = 0
    rlist[j].corrected_intensity_variance = 0

  p.finished('Generating %d reflections' % params.nrefl)
  shoebox.allocate(rlist)

  p = ProgressBar(title = 'Generating shoeboxes')

  for _, refl in enumerate(rlist):

    p.update(_ * 100.0 / params.nrefl)
    mask = refl.shoebox_mask

    if params.pixel_mask == 'precise':
      # flag everything as background: peak will me assigned later
      for j in range(len(mask)):
        mask[j] = mask_back
    elif params.pixel_mask == 'all':
      # flag we have no idea what anything is
      mask_none = MaskCode.Valid|MaskCode.Foreground|MaskCode.Background
      for j in range(len(mask)):
        mask[j] = mask_none
    elif params.pixel_mask == 'static':
      import itertools
      from scitbx.array_family import flex
      from math import sqrt
      x0 = params.spot_offset.x + params.shoebox_size.x / 2
      y0 = params.spot_offset.x + params.shoebox_size.y / 2
      z0 = params.spot_offset.x + params.shoebox_size.z / 2
      sx = params.mask_nsigma * params.spot_size.x
      sy = params.mask_nsigma * params.spot_size.y
      sz = params.mask_nsigma * params.spot_size.z

      # The x, y, z indices
      z, y, x = zip(*itertools.product(*(range(n) for n in mask.all())))
      xyz = flex.vec3_double(flex.double(x), flex.double(y), flex.double(z))

      # Calculate SUM(((xj - xj0) / sxj)**2) for each element
      xyz0 = (x0, y0, z0)
      isxyz = (1.0/sx, 1.0/sy, 1.0/sz)
      dxyz = sum([(x * isx)**2 for x, isx in
        zip(((xyz - xyz0) * rotation).parts(), isxyz)])

      # Set the mask values
      index = dxyz <= 1.0
      index.reshape(mask.accessor())
      mask.set_selected(index, MaskCode.Valid | MaskCode.Foreground)
      mask.set_selected(index != True, MaskCode.Valid | MaskCode.Background)

    sbox = refl.shoebox

    # reflection itself, including setting the peak region if we're doing that
    # FIXME use flex arrays to make the rotation bit more efficient as this is
    # now rather slow...

    counts_true = 0
    for j in range(params.counts):
      _x = random.gauss(0, params.spot_size.x)
      _y = random.gauss(0, params.spot_size.y)
      _z = random.gauss(0, params.spot_size.z)

      Rxyz = rotation * matrix.col((_x, _y, _z)).elems

      x = int(Rxyz[0] + params.spot_offset.x + params.shoebox_size.x / 2)
      y = int(Rxyz[1] + params.spot_offset.y + params.shoebox_size.y / 2)
      z = int(Rxyz[2] + params.spot_offset.z + params.shoebox_size.z / 2)

      if x < 0 or x >= params.shoebox_size.x:
        continue
      if y < 0 or y >= params.shoebox_size.y:
        continue
      if z < 0 or z >= params.shoebox_size.z:
        continue
      sbox[z, y, x] += 1
      counts_true += 1
      if params.pixel_mask == 'precise':
        mask[z, y, x] = mask_peak

    refl.intensity = counts_true

    # background:flat; FIXME can replace this with Poissonian added random
    # number, which will give the same answer... also would like to find some
    # reasonable way of gathering random numbers drawn from an inclined
    # distribution...

    for j in range(params.background * len(sbox)):
      x = random.randint(0, params.shoebox_size.x - 1)
      y = random.randint(0, params.shoebox_size.y - 1)
      z = random.randint(0, params.shoebox_size.z - 1)
      sbox[z, y, x] += 1

  p.finished('Generating %d shoeboxes' % params.nrefl)

  return rlist

def background_xds(rlist):
  from dials.algorithms.background import XdsSubtractor
  background = XdsSubtractor()
  background(None, None, rlist)
  return

def background_inclined(rlist):
  from dials.algorithms.background import InclinedSubtractor
  background = InclinedSubtractor()
  background(None, None, rlist)
  return

def integrate_3d_summation(rlist):
  from dials.algorithms.integration import Summation3d
  integration = Summation3d()
  integration(None, None, rlist)
  return

def main(params):
  rlist = simple_gaussian_spots(params)
  correct_intensities = [r.intensity for r in rlist]
  for r in rlist:
    r.intensity = 0

  if params.background_method == 'xds':
    background_xds(rlist)
  elif params.background_method == 'mosflm':
    assert(params.pixel_mask != 'all')
    background_inclined(rlist)

  integrate_3d_summation(rlist)

  integrated_intensities = [r.intensity for r in rlist]

  # now scan through the reflection list and find those where the integration
  # gave an apparently duff answer i.e. outside of 3 sigma from correct value

  from dials.model.data import ReflectionList
  import math

  overestimates = ReflectionList()
  underestimates = ReflectionList()

  for j, (c, i) in enumerate(zip(correct_intensities, integrated_intensities)):
    sigma = math.sqrt(c)
    if math.fabs(c - i) < 3 * sigma:
      continue
    if i > params.counts:
      overestimates.append(rlist[j])
    else:
      underestimates.append(rlist[j])

  print '%d overestimates, %d underestimates' % (len(overestimates),
                                                 len(underestimates))

  # now pickle these, perhaps

  import cPickle as pickle

  if params.output.under:
    pickle.dump(underestimates, open(params.output.under, 'w'))
  if params.output.over:
    pickle.dump(overestimates, open(params.output.over, 'w'))
  if params.output.all:
    pickle.dump(rlist, open(params.output.all, 'w'))

if __name__ == '__main__':
  import sys

  # FXIME use phil parameters for
  # nproc - done
  # FIXME add Phil parameters for
  # background plane method with B(x, y) = c + ax + by; don't know how to make
  #                                                   ; random dist. like this

  # FIXME parse Phil parameters better
  # FIXME pass in Phil file?
  working_phil = master_phil
  for arg in sys.argv[1:]:
    working_phil = working_phil.fetch(parse(arg))
  main(working_phil.extract())
