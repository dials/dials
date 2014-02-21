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
  i_to_compare = []

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
    i_to_compare.append(tmp_i)
    img[1, 1, 1] += tmp_i

    arr.data = img[:, :, :]
    arr.background = bkg[:, :, :]
    arr.mask = msk[:, :, :]

  #its = ref_table['intensity.raw.value']
  #i_var = ref_table['intensity.raw.variance']

  #for i in range(num_ref):
  #  its[i] = (i + 1) * 11
  #  i_var[i] = (i + 1) * 12

  from dials.algorithms.background.inclined_background_subtractor \
   import layering_and_background_plane
  layering_and_background_plane(ref_table)


  from dials.algorithms.integration.summation2d \
   import  flex_2d_layering_n_integrating
  flex_2d_layering_n_integrating(ref_table)

  result = "OK"
  resl_its = ref_table['intensity.raw.value']
  for n_its in range(len(resl_its)):
    if(resl_its[n_its] == i_to_compare[n_its]):
      print "Ok ", n_its
    else:
      print "Wrong num", n_its
      result = "wrong"

  to_be_fixed = '''
  for r in rlist:

    if r.intensity > 50 - rnd_siz and r.intensity < 50 + rnd_siz:
      print "Summation integration  ...  OK"

    else:
      print "Summation integration algorithm is not giving the espected result"
      result = "wrong number"
  '''

  return result


if __name__ == '__main__':
  for i in range(5):
    res=run()
    print res
