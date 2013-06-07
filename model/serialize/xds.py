#!/usr/bin/env python
#
# dials.model.serialize.xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def to_crystal(filename):
    ''' Get the crystal model from the xparm file

    Params:
        filename The xparm/or integrate filename

    Return:
        The crystal model

    '''
    from rstbx.cftbx.coordinate_frame_converter import \
        coordinate_frame_converter
    from dials.model.experiment import Crystal

    # Get the real space coordinate frame
    cfc = coordinate_frame_converter(filename)
    real_space_a = cfc.get('real_space_a')
    real_space_b = cfc.get('real_space_b')
    real_space_c = cfc.get('real_space_c')
    space_group = cfc.get('space_group_number')
    mosaicity = cfc.get('mosaicity')

    # Return the crystal model
    return Crystal(
        real_space_a=real_space_a,
        real_space_b=real_space_b,
        real_space_c=real_space_c,
        sg=space_group,
        mosaicity=mosaicity)
