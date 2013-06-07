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
    U = crystal.get_U()
    real_space_a = (U[0], U[3], U[6])
    real_space_b = (U[1], U[4], U[7])
    real_space_c = (U[2], U[5], U[8])

    # Get the space group number
    space_group = crystal.get_space_group().info().type().number()

    # Return the information as a python dictionary
    return OrderedDict([
        ('__id__', 'crystal'),
        ('real_space_a', real_space_a),
        ('real_space_b', real_space_b),
        ('real_space_c', real_space_c),
        ('space_group', space_group)])

def crystal_from_dict(d):
    ''' Convert the dictionary to a crystal model

    Params:
        d The dictionary of parameters

    Returns:
        The crystal model

    '''
    from dials.model.experiment.crystal_model.crystal import Crystal

    # If None, return None
    if d == None:
        return None

    # Check the version and id
    if str(d['__id__']) != "crystal":
        raise ValueError("\"__id__\" does not equal \"crystal\"")

    # Try to get the crystal model from the dictionary
    real_space_a = d['real_space_a']
    real_space_b = d['real_space_b']
    real_space_c = d['real_space_c']
    space_group  = d['space_group']
    return Crystal(real_space_a, real_space_b, real_space_c, space_group)
