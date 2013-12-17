#  test_integration_summation_2d
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
import random
rnd_siz = 10
def run():
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
  from scitbx.array_family import flex

  r = Reflection()
  r.shoebox = data2d
  r.shoebox_mask = mask2d
  r.shoebox_background = background2d
  rlist = ReflectionList()
  rlist.append(r)

  from dials.algorithms.background.inclined_background_subtractor \
   import layering_and_background_plane
  layering_and_background_plane(rlist)


  from dials.algorithms.integration.summation2d \
   import  flex_2d_layering_n_integrating
  flex_2d_layering_n_integrating(rlist)

  for r in rlist:

    if r.intensity > 50 - rnd_siz and r.intensity < 50 + rnd_siz:
      print "Summation integration  ...  OK"
      result = "OK"
    else:
      print "Summation integration algorithm is not giving the espected result"
      result = "wrong number"

  return result


if __name__ == '__main__':
  for i in range(5):
    res=run()
    print res
