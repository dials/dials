#!/usr/bin/env python
#
# dials.model.serialize.beam.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def beam_to_dict(beam):
    ''' Convert the beam model to a dictionary

    Params:
        beam The beam model

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    return OrderedDict([
        ('direction', beam.get_direction()),
        ('wavelength', beam.get_wavelength())])

def beam_from_dict(d):
    ''' Convert the dictionary to a beam model

    Params:
        d The dictionary of parameters

    Returns:
        The beam model

    '''
    from dials.model.experiment import Beam
    return Beam(tuple(d['direction']),
                float(d['wavelength']))
