#
# mosflm_like.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from dials.array_family import flex

def from_3D_to_2D_projection(shoebox, background):
  from dials.algorithms.integration import simple_2d_add
  if shoebox.all()[0] == 1:
    #print "No need for adding 3d to convert"
    data2d = shoebox[0:1, :, :]
    background2d = background[0:1, :, :]
    data2d.reshape(flex.grid(shoebox.all()[1:]))
    background2d.reshape(flex.grid(background.all()[1:]))

  else:
    #print "shoebox.all()[0] =", shoebox.all()[0]
    data2d_tot = flex.double(flex.grid(shoebox.all()[1:]),0.0)
    background2d_tot = flex.double(flex.grid(background.all()[1:]),0.0)
    for z_frm in range(shoebox.all()[0]):
      dada2d_to_add = shoebox[z_frm:z_frm + 1, :, :]
      dada2d_to_add.reshape(flex.grid(shoebox.all()[1:]))
      data2d_tot = simple_2d_add(data2d_tot, dada2d_to_add)

      background2d_to_add = background[z_frm:z_frm + 1, :, :]
      background2d_to_add.reshape(flex.grid(background.all()[1:]))
      background2d_tot = simple_2d_add(background2d_tot, background2d_to_add)

    data2d = data2d_tot[:,:]
    background2d = background2d_tot[:, :]

  return data2d, background2d


def from_3D_to_2D_mask_projection(mask):
  from dials.algorithms.integration import mask_add_2d
  if mask.all()[0] == 1:
    mask2d = mask[0:1, :, :]
    mask2d.reshape(flex.grid(mask.all()[1:]))
  else:
    mask2d_tot = flex.int(flex.grid(mask.all()[1:]),1)
    for z_frm in range(mask.all()[0]):
      mask2d_to_add = mask[z_frm:z_frm + 1, :, :]
      mask2d_to_add.reshape(flex.grid(mask.all()[1:]))
      mask2d_tot = mask_add_2d( mask2d_tot, mask2d_to_add)
    mask2d = mask2d_tot[:,:]


  return mask2d


