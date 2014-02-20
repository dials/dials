
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
  from scitbx.array_family import flex
  from dials.algorithms.integration import raw_2d_cut
  print "performing summation integration .... "


  imported_code = '''
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
  '''

  from dials.array_family import flex


  row_of_shoebox = ref_table['shoebox']
  row_of_its = ref_table['intensity.raw.value']
  row_of_var = ref_table['intensity.raw.variance']
  for row_num in range(ref_table.nrows()):
      local_shoebox = row_of_shoebox[row_num]


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

        reslt = raw_2d_cut(data2d, mask2d, background2d)
        print "reslt[0] =", reslt[0]
        print "reslt[1] =", reslt[1]

        i_r += reslt[0]
        i_v += reslt[1]

      row_of_its[row_num] = i_r
      row_of_var[row_num] = i_v
  old_code = '''
  for ref in reflections:
    if ref.is_valid():
      shoebox = ref.shoebox
      mask = ref.shoebox_mask
      background = ref.shoebox_background
      ref.intensity = 0.0
      ref.intensity_variance = 0.0
      for i in range(shoebox.all()[0]):
        data2d = shoebox[i:i + 1, :, :]
        mask2d = mask[i:i + 1, :, :]
        background2d = background[i:i + 1, :, :]

        data2d.reshape(flex.grid(shoebox.all()[1:]))
        mask2d.reshape(flex.grid(shoebox.all()[1:]))
        background2d.reshape(flex.grid(shoebox.all()[1:]))

        reslt = raw_2d_cut(data2d, mask2d, background2d)

        ref.intensity += reslt[0]
        ref.intensity_variance += reslt[1]
  '''

  print "summation integration .... done"

  return ref_table
