#
# dials.algorithms.background.flat_subtractor.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero "luiso" & James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class FlatSubtractor(object):
  ''' The Flat background subtractor '''

  def __init__(self, **kwargs):
    pass

  def __call__(self, reflections):
    layering_and_background_avg(reflections)
    return reflections

def layering_and_background_avg(reflections):
  from dials.algorithms.background import flat_background_flex_2d
  from scitbx.array_family import flex

  shoeboxes = reflections['shoebox']
  for shoebox in shoeboxes:
    #if ref.is_valid():
      data = shoebox.data
      mask = shoebox.mask
      background = shoebox.background
      for i in range(data.all()[0]):
        data2d = data[i:i + 1, :, :]
        mask2d = mask[i:i + 1, :, :]
        data2d.reshape(flex.grid(data.all()[1:]))
        mask2d.reshape(flex.grid(data.all()[1:]))
        background2d = flat_background_flex_2d(data2d, mask2d)
        background2d.reshape(flex.grid(1, background2d.all()[0], background2d.all()[1]))
        background[i:i + 1, :, :] = background2d.as_double()

  return reflections
