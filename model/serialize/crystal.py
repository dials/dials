from __future__ import division
#!/usr/bin/env python
#
# dials.model.serialize.crystal.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def crystal_to_dict(crystal):
    ''' Convert the crystal model to a dictionary

    Params:
        crystal The crystal model

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict

    # Get the real space vectors
    A = crystal.get_A().inverse()
    real_space_a = (A[0], A[1], A[2])
    real_space_b = (A[3], A[4], A[5])
    real_space_c = (A[6], A[7], A[8])

    # Get the space group Hall symbol
    hall = crystal.get_space_group().info().type().hall_symbol()

    # Get the mosaicity
    mosaicity = crystal.get_mosaicity()

    # Return the information as a python dictionary
    return OrderedDict([
        ('__id__', 'crystal'),
        ('real_space_a', real_space_a),
        ('real_space_b', real_space_b),
        ('real_space_c', real_space_c),
        ('space_group_hall_symbol', hall),
        ('mosaicity', mosaicity)])

def crystal_from_dict(d):
    ''' Convert the dictionary to a crystal model

    Params:
        d The dictionary of parameters

    Returns:
        The crystal model

    '''
    from dials.model.experiment.crystal_model.crystal import Crystal

    # If None, return None
    if d is None:
        return None

    # Check the version and id
    if str(d['__id__']) != "crystal":
        raise ValueError("\"__id__\" does not equal \"crystal\"")

    # Try to get the crystal model from the dictionary
    real_space_a = d['real_space_a']
    real_space_b = d['real_space_b']
    real_space_c = d['real_space_c']
    space_group  = "Hall:" + d['space_group_hall_symbol']
    mosaicity    = d['mosaicity']
    return Crystal(real_space_a, real_space_b, real_space_c,
                   space_group_symbol=space_group,
                   mosaicity=mosaicity)
