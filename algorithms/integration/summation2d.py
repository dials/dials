
#
# integrate2d.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luiso & James
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Summation2d(object):
  '''A class to perform 2D integration'''

  def __init__(self, **kwargs):
    '''Initialise algorithm.'''
    pass
  def __call__(self, sweep, crystal, reflections, reference = None):
    '''Process the reflections.'''
    self.integrate(reflections)
    return reflections

  def integrate(self, reflections):
    flex_2d_layering_n_integrating(reflections)

def flex_2d_layering_n_integrating(ref_table):
  '''
  integrated via summation-integration by layering ,
  frame per frame, each reflection in the table
  '''
  from scitbx.array_family import flex
  from dials.algorithms.integration import raw_2d_cut
  print "Performing summation integration .... "

  from dials.array_family import flex
  # extracting needed info from table
  col_of_shoebox = ref_table['shoebox']
  col_of_its = ref_table['intensity.raw.value']
  col_of_var = ref_table['intensity.raw.variance']
  for row_num in range(ref_table.nrows()):
      local_shoebox = col_of_shoebox[row_num]

      i_r = 0
      i_v = 0
      #if ref.is_valid():
      data = local_shoebox.data
      mask = local_shoebox.mask
      background = local_shoebox.background

      for i in range(data.all()[0]):

        data2d = data[i:i + 1, :, :]
        mask2d = mask[i:i + 1, :, :]
        background2d = background[i:i + 1, :, :]

        data2d.reshape(flex.grid(data.all()[1:]))
        mask2d.reshape(flex.grid(data.all()[1:]))
        background2d.reshape(flex.grid(data.all()[1:]))

        # boosted c++ function
        reslt = raw_2d_cut(data2d, mask2d, background2d)

        i_r += reslt[0]
        i_v += reslt[1]

      col_of_its[row_num] = i_r
      col_of_var[row_num] = i_v

  print "summation integration      ....       Done"

  return ref_table
