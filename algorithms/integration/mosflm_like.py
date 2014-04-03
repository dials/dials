#
# mosflm_like.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luiso & James
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class MosflmProfileFitting:
  '''A class to perform 2D mosflm-like profile fitting'''

  def __init__(self, experiment, nblocks):
    self.experiment = experiment
    self.nblocks = nblocks


  #def __call__(self, sweep, crystal, ref_table, reference = None):
  #def __call__(self, sweep, ref_table, reference = None):
  def __call__(self, ref_table, reference = None):

    from dials.array_family import flex
    from dials.algorithms.integration import flex_2d_layering_n_integrating
    from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller

    flex_2d_layering_n_integrating(ref_table)

    xyz = ref_table['xyzcal.px']
    index = sorted(range(len(ref_table)), key=lambda i: xyz[i][2])
    ref_table.reorder(flex.size_t(index))

    new_block_way = '''
    nz_blocks = 5
    dp_lng =len(ref_table)
    zblock_size = dp_lng / nz_blocks
    for block_z_num in range(nz_blocks):
      z_blocks_start = int(block_z_num * zblock_size)
      z_blocks_end = int((block_z_num + 1) * zblock_size)
      local_ref_table = ref_table[z_blocks_start:z_blocks_end]

      xmax, ymax = self.experiment.detector[0].get_image_size()
      local_ref_table = mosflm_caller(local_ref_table, xmax, ymax, self.nblocks)

      ref_table[z_blocks_start:z_blocks_end] = local_ref_table
    #'''

    #old_way = '''
    xmax, ymax = self.experiment.detector[0].get_image_size()
    ref_table = mosflm_caller(ref_table, xmax, ymax, self.nblocks)
    #'''
