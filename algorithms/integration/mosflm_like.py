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



    imagin_stuff = '''
    ###############################################################################
    t_intensity = ref_table['intensity.sum.value']
    old_i_table = t_intensity[:]
    ###############################################################################
    #'''


    #new_multiple_block_way = '''

    data_range_tst = self.experiment.scan.get_oscillation_range()
    #print "self.experiment.scan.get_oscillation_range =", data_range_tst[0]
    #print "self.experiment.scan.get_oscillation_range =", data_range_tst[1]

    nz_blocks = int(abs(data_range_tst[1] - data_range_tst[0]) / 5)
    if(nz_blocks < 1):
      nz_blocks = 1
    print "N(z)blocks = ", nz_blocks


    dp_lng =len(ref_table)
    zblock_size = dp_lng / nz_blocks
    for block_z_num in range(nz_blocks):
      print "Block of images num ", block_z_num
      z_blocks_start = int(block_z_num * zblock_size)
      z_blocks_end = int((block_z_num + 1) * zblock_size)
      local_ref_table = ref_table[z_blocks_start:z_blocks_end]

      xmax, ymax = self.experiment.detector[0].get_image_size()
      local_ref_table = mosflm_caller(local_ref_table, xmax, ymax, self.nblocks)

      ref_table[z_blocks_start:z_blocks_end] = local_ref_table

    #'''


    data_viewing = '''
    print "self.experiment.goniometer =",  dir(self.experiment.goniometer)
    print "self.experiment.scan =",  dir(self.experiment.scan)
    print "self.experiment.imageset =",  dir(self.experiment.imageset)

    [experiment.goniometer] = [ 'from_dict', 'get_fixed_rotation',
     'get_rotation_axis', 'set_fixed_rotation', 'set_rotation_axis', 'to_dict']

    [experiment.scan] = [ 'from_dict', 'get_angle_from_array_index'
    , 'get_angle_from_image_index', 'get_array_index_from_angle'
    , 'get_array_indices_with_angle', 'get_array_range', 'get_epochs'
    , 'get_exposure_times', 'get_image_epoch', 'get_image_index_from_angle'
    , 'get_image_indices_with_angle', 'get_image_oscillation', 'get_image_range'
    , 'get_num_images', 'get_oscillation', 'get_oscillation_range', 'is_angle_valid'
    , 'is_array_index_valid', 'is_image_index_valid', 'set_epochs'
    , 'set_exposure_times', 'set_image_range', 'set_oscillation', 'to_dict']

    (experiment.imageset] = [ '_beam', '_detector', '_get_data_range', '_goniometer'
    , '_image_index', '_indices', '_models', '_reader', '_scan', '_to_array_all'
    , '_to_array_w_range', '_truncate_range', 'complete_set', 'get_array_range'
    , 'get_beam', 'get_detector', 'get_detectorbase', 'get_goniometer'
    , 'get_image_models', 'get_image_size', 'get_path', 'get_scan', 'get_template'
    , 'indices', 'is_valid', 'paths', 'reader', 'set_beam', 'set_detector'
    , 'set_goniometer', 'set_scan', 'to_array']
    '''




    old_single_block_way = '''
    xmax, ymax = self.experiment.detector[0].get_image_size()
    ref_table = mosflm_caller(ref_table, xmax, ymax, self.nblocks)
    #'''


    imagin_stuff = '''
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