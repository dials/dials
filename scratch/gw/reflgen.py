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
counts = 0
  .type = int
background = 0
  .type = int
pixel_mask = *none static precise
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
""")

def simple_gaussian_spots(num_refl, signal, background):
  from dials.model.data import ReflectionList
  from dials.algorithms import shoebox
  import random
  import math

  rlist = ReflectionList(num_refl)

  from dials.util.command_line import ProgressBar
  p = ProgressBar(title = 'Generating reflections')

  for j in range(num_refl):
    p.update(j * 100.0 / num_refl)
    rlist[j].miller_index = (random.randint(0, 20),
                             random.randint(0, 20),
                             random.randint(0, 20))
    rlist[j].rotation_angle = 2 * math.pi * random.random()
    rlist[j].beam_vector = (0, 0, 0)
    rlist[j].image_coord_px = (0, 0)
    rlist[j].image_coord_mm = (0, 0)
    rlist[j].frame_number = 0
    rlist[j].panel_number = 0
    rlist[j].bounding_box = (0, 10, 0, 10, 0, 10)
    rlist[j].centroid_position = (0, 0, 0)
    rlist[j].centroid_variance = (0, 0, 0)
    rlist[j].centroid_sq_width = (0, 0, 0)
    rlist[j].intensity = 0
    rlist[j].intensity_variance = 0
    rlist[j].corrected_intensity = 0
    rlist[j].corrected_intensity_variance = 0

  p.finished('Generating %d reflections' % num_refl)
  shoebox.allocate(rlist)

  p = ProgressBar(title = 'Generating shoeboxes')

  for _, refl in enumerate(rlist):

    p.update(_ * 100.0 / num_refl)
    sbox = refl.shoebox

    x0 = 5
    y0 = 5
    z0 = 5

    # reflection itself, including setting the peak region

    counts_true = 0
    for j in range(signal):
      x = int(random.gauss(x0, 1))
      y = int(random.gauss(y0, 1))
      z = int(random.gauss(z0, 1))
      if x < 0 or x >= 10:
        continue
      if y < 0 or y >= 10:
        continue
      if z < 0 or z >= 10:
        continue
      sbox[z, y, x] += 1
      counts_true += 1

    refl.intensity = counts_true

    # background:flat

    for j in range(background * len(sbox)):
      x = random.randint(0, 9)
      y = random.randint(0, 9)
      z = random.randint(0, 9)
      sbox[z, y, x] += 1

  # set the mask - it's all good

  from dials.algorithms.shoebox import MaskCode
  mask_value = MaskCode.Valid|MaskCode.Foreground|MaskCode.Background
  for refl in rlist:
    mask = refl.shoebox_mask
    for j in range(len(mask)):
      mask[j] = mask_value
  p.finished('Generating %d shoeboxes' % num_refl)

  return rlist

def simple_gaussian_spots_proper_mask(num_refl, signal, background):
  from dials.model.data import ReflectionList
  from dials.algorithms import shoebox
  import random
  import math

  # generate mask and peak values

  from dials.algorithms.shoebox import MaskCode
  mask_peak = MaskCode.Valid|MaskCode.Foreground
  mask_back = MaskCode.Valid|MaskCode.Background

  rlist = ReflectionList(num_refl)

  from dials.util.command_line import ProgressBar
  p = ProgressBar(title = 'Generating reflections')

  for j in range(num_refl):
    p.update(j * 100.0 / num_refl)
    rlist[j].miller_index = (random.randint(0, 20),
                             random.randint(0, 20),
                             random.randint(0, 20))
    rlist[j].rotation_angle = 2 * math.pi * random.random()
    rlist[j].beam_vector = (0, 0, 0)
    rlist[j].image_coord_px = (0, 0)
    rlist[j].image_coord_mm = (0, 0)
    rlist[j].frame_number = 0
    rlist[j].panel_number = 0
    rlist[j].bounding_box = (0, 10, 0, 10, 0, 10)
    rlist[j].centroid_position = (0, 0, 0)
    rlist[j].centroid_variance = (0, 0, 0)
    rlist[j].centroid_sq_width = (0, 0, 0)
    rlist[j].intensity = 0
    rlist[j].intensity_variance = 0
    rlist[j].corrected_intensity = 0
    rlist[j].corrected_intensity_variance = 0

  p.finished('Generating %d reflections' % num_refl)
  shoebox.allocate(rlist)

  p = ProgressBar(title = 'Generating shoeboxes')

  for _, refl in enumerate(rlist):

    p.update(_ * 100.0 / num_refl)

    # flag everything as background
    mask = refl.shoebox_mask
    for j in range(len(mask)):
      mask[j] = mask_back

    sbox = refl.shoebox

    x0 = 5
    y0 = 5
    z0 = 5

    # reflection itself, including setting the peak region

    counts_true = 0
    for j in range(signal):
      x = int(random.gauss(x0, 1))
      y = int(random.gauss(y0, 1))
      z = int(random.gauss(z0, 1))
      if x < 0 or x >= 10:
        continue
      if y < 0 or y >= 10:
        continue
      if z < 0 or z >= 10:
        continue
      sbox[z, y, x] += 1
      mask[z, y, x] = mask_peak
      counts_true += 1

    refl.intensity = counts_true

    # background:flat

    for j in range(background * len(sbox)):
      x = random.randint(0, 9)
      y = random.randint(0, 9)
      z = random.randint(0, 9)
      sbox[z, y, x] += 1

  p.finished('Generating %d shoeboxes' % num_refl)

  return rlist

def integrate_xds_background_3d_summation(rlist):
  from dials.algorithms.background import XdsSubtractor
  from dials.algorithms.integration import Summation3d

  background = XdsSubtractor()
  background(None, None, rlist)
  integration = Summation3d()
  integration(None, None, rlist)
  return

def integrate_inclined_background_3d_summation(rlist):
  from dials.algorithms.background import InclinedSubtractor
  from dials.algorithms.integration import Summation3d

  background = InclinedSubtractor()
  background(None, None, rlist)
  integration = Summation3d()
  integration(None, None, rlist)
  return

def main(params):
  if params.pixel_mask == 'none':
    rlist = simple_gaussian_spots(params.nrefl, params.counts, params.background)
  elif params.pixel_mask == 'precise':
    rlist = simple_gaussian_spots_proper_mask(params.nrefl, params.counts, 
                                              params.background)
  correct_intensities = [r.intensity for r in rlist]
  for r in rlist:
    r.intensity = 0

  if params.background_method == 'xds':
    integrate_xds_background_3d_summation(rlist)
  elif params.background_method == 'mosflm':
    integrate_inclined_background_3d_summation(rlist)
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

  # FIXME add Phil parameters for
  # nrefl - done
  # nproc - done
  # size of box - done x y z
  # width of spot - done x y z
  # rotation of spot density - not done
  # background plane method with B(x, y) = c + ax + by; don't know how to make
  #                                                   ; random dist. like this
  # counts
  # background
  # background method
  # integration method
  # set proper background mask or no
  # set static mask

  working_phil = master_phil
  for arg in sys.argv[1:]:
    working_phil = working_phil.fetch(parse(arg))
  main(working_phil.extract())
