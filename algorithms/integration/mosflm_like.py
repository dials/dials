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

  def __call__(self, ref_table, reference = None):

    from dials.array_family import flex
    from dials.algorithms.integration import flex_2d_layering_n_integrating
    from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller

    flex_2d_layering_n_integrating(ref_table)

    xyz = ref_table['xyzcal.px']
    index = sorted(range(len(ref_table)), key=lambda i: xyz[i][2])
    ref_table.reorder(flex.size_t(index))




    #imagin_stuff = '''
    ###############################################################################
    t_intensity = ref_table['intensity.sum.value']
    old_i_table = t_intensity[:]
    ###############################################################################
    #'''




    #new_block_way = '''
    nz_blocks = 12
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

    old_way = '''
    xmax, ymax = self.experiment.detector[0].get_image_size()
    ref_table = mosflm_caller(ref_table, xmax, ymax, self.nblocks)
    #'''




    #imagin_stuff = '''
    ###############################################################################
    t_intensity = ref_table['intensity.prf.value']
    num_ref = len(t_intensity)
    paint_compare = []
    for i in range(num_ref):
      paint_compare.append([ old_i_table[i], t_intensity[i]])
    paint_compare_sort = sorted(paint_compare)
    import numpy
    data1d = numpy.zeros(num_ref, dtype = numpy.float64)
    new_data1d = numpy.zeros(num_ref, dtype = numpy.float64)
    for i in range(num_ref):
      data1d[i] = paint_compare_sort[i][0]
      new_data1d[i] = paint_compare_sort[i][1]

    from matplotlib import pylab
    pylab.plot(data1d)
    pylab.plot(new_data1d)
    pylab.show()
    ###############################################################################
    #'''