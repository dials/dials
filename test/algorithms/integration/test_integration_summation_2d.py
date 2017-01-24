#  test_integration_summation_2d
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

def run(i, imp):
  from random import randint
  from dials.array_family import flex

  #building a reflection table
  num_ref = 5
  ref_table = flex.reflection_table()

  shoebox = flex.shoebox(num_ref)
  ref_table['shoebox'] = shoebox

  intensity = flex.double(num_ref)
  ref_table['intensity.sum.value'] = intensity

  intensity_var = flex.double(num_ref)
  ref_table['intensity.sum.variance'] = intensity_var

  iterate = ref_table['shoebox']
  i_to_compare = []

  # bulding the shoebox with a desired content
  # which is a reflection with noise included

  n = 0
  for arr in iterate:
    img = flex.double(flex.grid(3, 3, 3))
    bkg = flex.double(flex.grid(3, 3, 3))
    msk = flex.int(flex.grid(3, 3, 3))
    for row in range(3):
      for col in range(3):
        for fra in range(3):
          img[row, col, fra] = row + col + fra + n * 9 + randint(0, i)
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

  # calling the functions that we need to test
  # first select the algorithm for background calculation

  if imp == "inclined":
    print "testing inclined_background_subtractor"
    from dials.algorithms.background.inclined_background_subtractor \
     import layering_and_background_plane
    layering_and_background_plane(ref_table)
  elif imp == "flat":
    print "testing flat_background_subtractor"
    from dials.algorithms.background.flat_background_subtractor \
     import layering_and_background_avg
    layering_and_background_avg(ref_table)
  elif imp == "curved":
    print "testing curved_background_subtractor"
    from dials.algorithms.background.curved_background_subtractor \
     import layering_and_background_modl
    layering_and_background_modl(ref_table)

  # no matter which algorithm was used for background calculation
  # the integration summation must remain compatible

  from dials.algorithms.integration.summation2d \
    import  flex_2d_layering_n_integrating
  flex_2d_layering_n_integrating(ref_table)

  # comparing results

  result = "OK"
  resl_its = ref_table['intensity.sum.value']
  resl_var = ref_table['intensity.sum.variance']
  for n_its in range(len(resl_its)):
    if resl_its[n_its] <= i_to_compare[n_its] + i and \
       resl_its[n_its] >= i_to_compare[n_its] - i and \
       resl_var[n_its] > resl_its[n_its]:
      print "Ok ", n_its
    else:
      print "Wrong num", n_its

      print "i =", i
      print "resl_its[n_its] =", resl_its[n_its]
      print "i_to_compare[n_its] =", i_to_compare[n_its]
      print "resl_var[n_its] =", resl_var[n_its]

      result = "wrong"
      raise RuntimeError('wrong result')
  return result


if __name__ == '__main__':
  for i in range(5):
    res1 = run(i, "flat")
    print res1
    res2 = run(i, "inclined")
    print res2
    res3 = run(i, "curved")
    print res3
