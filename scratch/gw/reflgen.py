from __future__ import division

def reflgen(num_refl):

  from dials.model.data import ReflectionList
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

def simple_gaussian_spots(num_refl, signal, background):
  from dials.model.data import ReflectionList, Reflection
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

    # reflection itself

    counts_true = 0
    for j in range(signal):
      x = int(round(random.gauss(x0, 1)))
      y = int(round(random.gauss(y0, 1)))
      z = int(round(random.gauss(z0, 1)))
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

def integrate_xds_background_3d_summation(rlist):
  from dials.algorithms.background import XdsSubtractor
  from dials.algorithms.integration import Summation3d

  background = XdsSubtractor()
  background(None, None, rlist)
  integration = Summation3d()
  integration(None, None, rlist)
  return

def meansd(values):
  import math
  mean = sum(values) / len(values)
  var = sum([(v - mean) ** 2 for v in values]) / (len(values) - 1)
  return mean, math.sqrt(var)

def main(argv):
  nrefl, counts, background = map(int, argv)
  rlist = simple_gaussian_spots(nrefl, counts, background)
  correct_intensities = [r.intensity for r in rlist]
  for r in rlist:
    r.intensity = 0
  integrate_xds_background_3d_summation(rlist)
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
    if i > counts:
      overestimates.append(rlist[j])
    else:
      underestimates.append(rlist[j])

  print '%d overestimates, %d underestimates' % (len(overestimates),
                                                 len(underestimates))

  # now pickle these

  import cPickle as pickle

  pickle.dump(underestimates,
              open('%d_%d_under.pickle' % (counts, background), 'w'))

  pickle.dump(overestimates,
              open('%d_%d_over.pickle' % (counts, background), 'w'))

if __name__ == '__main__':
  import sys
  main(sys.argv[1:])
