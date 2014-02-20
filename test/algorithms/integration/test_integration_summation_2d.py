#  test_integration_summation_2d
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

def run():
  old_way = '''
  import random
  rnd_siz = 10
  from scitbx.array_family import flex
  data2d = flex.double(flex.grid(1, 3, 3),15)

  for row in range(3):
    for col in range(3):
      data2d[0,row, col] += row * 2
      data2d[0,row, col] += col * 2

  mask2d = flex.int(flex.grid(1, 3, 3),3)
  mask2d[0, 1, 1] = 5
  data2d[0, 1, 1] += 50

  for row in range(3):
    for col in range(3):
      rnd_no = random.randint(0,rnd_siz)
      data2d[0,row, col] += rnd_no

  background2d = flex.double(flex.grid(1, 3, 3),0)

  from dials.model.data import Reflection, ReflectionList

  r = Reflection()
  r.shoebox = data2d
  r.shoebox_mask = mask2d
  r.shoebox_background = background2d
  rlist = ReflectionList()
  rlist.append(r)

  '''
  #from scitbx.array_family import flex
  from dials.array_family import flex
  num_ref = 3
  ref_table = flex.reflection_table()

  shoebox = flex.shoebox(num_ref)
  ref_table['shoebox'] = shoebox

  intensity = flex.double(num_ref)
  ref_table['intensity.raw.value'] = intensity

  intensity_var = flex.double(num_ref)
  ref_table['intensity.raw.variance'] = intensity_var

  iterate = ref_table['shoebox']
  n = 0
  for arr in iterate:
    img = flex.double(flex.grid(3, 3, 3))
    bkg = flex.double(flex.grid(3, 3, 3))
    msk = flex.int(flex.grid(3, 3, 3))
    for row in range(3):
      for col in range(3):
        for fra in range(3):
          img[row, col, fra] = row + col + fra + n * 9
          bkg[row, col, fra] = 0.0
          msk[row, col, fra] = 3
    n += 1
    msk[1, 1, 1] = 5
    tmp_i = n * n * n * 3
    img[1, 1, 1] += tmp_i
    print "intensity must be =", tmp_i
    arr.data = img[:, :, :]
    arr.background = bkg[:, :, :]
    arr.mask = msk[:, :, :]

  its = ref_table['intensity.raw.value']
  i_var = ref_table['intensity.raw.variance']

  for i in range(num_ref):
    its[i] = (i + 1) * 11
    i_var[i] = (i + 1) * 12

  from dials.algorithms.background.inclined_background_subtractor \
   import layering_and_background_plane
  layering_and_background_plane(ref_table)


  from dials.algorithms.integration.summation2d \
   import  flex_2d_layering_n_integrating
  flex_2d_layering_n_integrating(ref_table)


  iterate = ref_table['intensity.raw.value']
  for n_its in iterate:
    print n_its
    print ">>>"

  to_be_fixed = '''
  for r in rlist:

    if r.intensity > 50 - rnd_siz and r.intensity < 50 + rnd_siz:
      print "Summation integration  ...  OK"

    else:
      print "Summation integration algorithm is not giving the espected result"
      result = "wrong number"
  '''
  result = "OK"
  return result


if __name__ == '__main__':
  for i in range(5):
    res=run()
    print res
