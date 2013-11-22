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

def start(num_refl):
  from dials.model.data import ReflectionList, Reflection
  from dials.algorithms import shoebox
  import random
  import math

  rlist = ReflectionList(num_refl)

  for j in range(num_refl):
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

  shoebox.allocate(rlist)

  for refl in rlist:

    counts = random.randint(0, 10000)
    background = random.randint(0, 10000)
    
    # FIXME in here fake up background and peak and mask but don't assign
    # the "true" values here unless we want to test things? Maybe that
    # would be best ... make sure we can reproduce the observations?

    sbox = refl.shoebox

    x0 = 5
    y0 = 5
    z0 = 5

    # reflection itself

    counts_true = 0
    for j in range(counts):
      x = int(round(random.gauss(x0, 1)))
      y = int(round(random.gauss(x0, 1)))
      z = int(round(random.gauss(x0, 1)))
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

    for j in range(background):
      x = random.randint(0, 9)
      y = random.randint(0, 9)
      z = random.randint(0, 9)
      sbox[z, y, x] += 1

  return rlist

def print_shoebox(sbox):
  ndim = sbox.nd()
  dims = sbox.focus()
  
  print '-' * 80
  
  for k in range(dims[0]):
    for j in range(dims[1]):
      for i in range(dims[0]):
        print '%4d' % int(sbox[k, j, i]),
      print
    print '-' * 80
      
  return

def print_summary(refl):
  print 'Claimed / total intensity: %9.2f %9.2f' % (refl.intensity, 
                                                    sum(refl.shoebox))

def main():
  profiles = start(100)
  for p in profiles:
    print_summary(p)
    
if __name__ == '__main__':
  main()
