#
# __init__.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst & Luis Fuentes Montero
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
# from dials.algorithms.peak_finding.spot_finder_lui import SpotFinderLui
from __future__ import division
from dials.interfaces.peak_finding import SpotFinderInterface
from dials.algorithms.peak_finding.lui_find_peak import *
class SpotFinderLui(SpotFinderInterface):
    '''An interface specification for spot finding classes.'''

    def __init__(self, **kwargs):
        '''Initialise the algorithm with some parameters.

        Params:
            kwargs Key word arguments

        '''
        pass

    def __call__(self, sweep, times, shift):
        '''The main function of the spot finder. Select the pixels from
        the sweep and then group the pixels into spots. Return the data
        in the form of a reflection list.

        Params:
            sweep The sweep object
            some parameters to

        Returns:
            The reflection list

        '''
        reflection_list = do_all_3d(sweep, times, shift)

        return reflection_list

def do_all_3d(sweep, times, shift):
    from scitbx.array_family import flex
    import numpy
    array_3d = sweep.to_array()
    data3d = array_3d.as_numpy_array()
    n_frm = numpy.size(data3d[:, 0:1, 0:1])
    n_row = numpy.size(data3d[0:1, :, 0:1])
    n_col = numpy.size(data3d[0:1, 0:1, :])

    # print "n_frm,n_row,n_col", n_frm, n_row, n_col

    dif3d = numpy.zeros_like(data3d)

    n_blocks_x = 5
    n_blocks_y = 12
    col_block_size = n_col / n_blocks_x
    row_block_size = n_row / n_blocks_y

    for frm_tmp in range(n_frm):
        for tmp_block_x_pos in range(n_blocks_x):
            for tmp_block_y_pos in range(n_blocks_y):

                col_from = int(tmp_block_x_pos * col_block_size)
                col_to = int((tmp_block_x_pos + 1) * col_block_size)
                row_from = int(tmp_block_y_pos * row_block_size)
                row_to = int((tmp_block_y_pos + 1) * row_block_size)

                tmp_dat2d = numpy.copy(data3d[frm_tmp, row_from:row_to, col_from:col_to])
                tmp_dif = find_mask_2d(tmp_dat2d, times, shift)
                dif3d[frm_tmp, row_from:row_to, col_from:col_to] = tmp_dif

    dif_3d_ext = find_ext_mask_3d(dif3d)

    x_from_lst, x_to_lst, y_from_lst, y_to_lst, z_from_lst, z_to_lst = find_bound_3d(dif_3d_ext)


    reflection_list = _create_reflection_list(x_from_lst, x_to_lst, y_from_lst, y_to_lst, z_from_lst, z_to_lst)

    for i, rf_lst in enumerate(reflection_list):
        rf_lst.shoebox = flex.int(numpy.copy(data3d[          \
                                  z_from_lst[i]: z_to_lst[i], \
                                  y_from_lst[i]: y_to_lst[i], \
                                  x_from_lst[i]: x_to_lst[i] ]))
        rf_lst.shoebox_mask = flex.int(numpy.copy(dif_3d_ext[ \
                                  z_from_lst[i]: z_to_lst[i], \
                                  y_from_lst[i]: y_to_lst[i], \
                                  x_from_lst[i]: x_to_lst[i] ]))
    from dials.algorithms.centroid.toy_centroid_Lui import toy_centroid_lui
    centroid = toy_centroid_lui(reflection_list)
    reflection_list = centroid.get_reflections()

    return reflection_list

def _create_reflection_list(x_from_lst, x_to_lst, y_from_lst, y_to_lst, z_from_lst, z_to_lst):

    '''Create a reflection list from the spot data.
    
    Params:
        coords The pixel coordinates
        values The pixel values
        spots The pixel->spot mapping
        bbox The spot bounding boxes
        cpos The centroid position
        cvar The centroid variance
        index The list of valid indices
    
    Returns:
        A list of reflections
    '''
    from dials.model.data import Reflection, ReflectionList
    from dials.algorithms.integration import allocate_reflection_profiles

    # Create the reflection list
    length = len(x_from_lst)

    reflection_list = ReflectionList(length)

    for i, rf_lst in enumerate(reflection_list):
        bbox = [x_from_lst[i], x_to_lst[i], y_from_lst[i], y_to_lst[i], z_from_lst[i], z_to_lst[i]]
        rf_lst.bounding_box = bbox




    return reflection_list
