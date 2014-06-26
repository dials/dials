#!/usr/bin/env python
#
# flatten_shoebox.py
#
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from dials.array_family import flex

def flatten_shoebox(shoebox, background, mask):
  '''Compress reflection shoebox in rotation direction to make a pseudo-fully
  recorded reflection. N.B. need to handle background masks gracefully. Returns
  flattened shoebox arrays.'''

  # test dimensions of input: if depth only 1 slice then return copy of the data
  # as two-dimensional image

  array_dimensions = shoebox.all()
  assert array_dimensions == background.all()

  if array_dimensions[0] == 1:
    import copy
    return copy.deepcopy(shoebox).reshape(array_dimensions[1:]), \
      copy.deepcopy(background).reshape(array_dimensions[1:]), \
      copy.deepcopy(mask).reshape(array_dimensions[1:])

  # we have a data array which needs to be actually flattened: perform this
  # operation

  from dials.algorithms.shoebox.MaskCode import Foreground, Background, Valid

  # first scan through the slices of the shoebox mask and work out the slices
  # which have valid Foreground pixels in (so we don't add noise)

  # then create a target array and add in the pixel values for those slices which
  # should contribute for the shoebox and background

  shoebox_out = flex.double(flex.grid(array_dimensions[1:]), 0.0)
  background_out = flex.double(flex.grid(array_dimensions[1:]), 0.0)
  mask_out = flex.int(flex.grid(array_dimensions[1:]), 0)

  # we will sum these at the end to reconstruct a proper mask: note carefully
  # that we 'and' background and valid and 'or' foreground

  foreground = flex.bool(flex.grid(array_dimensions[1:]), False)
  background = flex.bool(flex.grid(array_dimensions[1:]), True)
  valid = flex.bool(flex.grid(array_dimensions[1:]), True)

  # FIXME this should be the correct slices not the whole block
  for k in range(array_dimensions[0]):
    shoebox_out += shoebox[k:k+1, :, :].reshape(array_dimensions[1:])
    background_out += background[k:k+1, :, :].reshape(array_dimensions[1:])

    # currently need to work out the mask on a pixel-by-pixel basis
    for j in range(array_dimensions[1]):
      for i in range(array_dimensions[2]):
        foreground[j, i] = foreground[j, i] | (mask[k, j, i] & Foreground)
        background[j, i] = background[j, i] & (mask[k, j, i] & Background)
        valid[j, i] = valid[j, i] & (mask[k, j, i] & Valid)

  # finally work out the appropriate mask for this new flattened shoebox by adding
  # the foreground, background and valid masks...

  mask_out = Foreground * foreground + Background * background + Valid * valid

  return shoebox_out, background_out, mask_out
