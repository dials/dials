#!/usr/bin/env python
#
# dials.model.serialize.detector.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def panel_to_dict(panel):
    ''' Convert the panel model to a dictionary

    Params:
        panel The panel model

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    return OrderedDict(type=panel.get_type(),
                       fast_axis=panel.get_fast_axis(),
                       slow_axis=panel.get_slow_axis(),
                       origin=panel.get_origin(),
                       pixel_size=panel.get_pixel_size(),
                       image_size=panel.get_image_size(),
                       trusted_range=panel.get_trusted_range())

def panel_from_dict(d):
    ''' Convert the dictionary to a panel model

    Params:
        d The dictionary of parameters

    Returns:
        The panel model

    '''
    from dials.model.experiment import Panel
    return Panel(str(d['type']),
                 tuple(d['fast_axis']),
                 tuple(d['slow_axis']),
                 tuple(d['origin']),
                 tuple(d['image_size']),
                 tuple(d['pixel_size']),
                 tuple(d['trusted_range']))

def detector_to_dict(detector):
    ''' Convert the detector model to a dictionary

    Params:
        detector The detector model

    Returns:
        A dictionary of the parameters

    '''
    return [panel_to_dict(p) for p in detector]

def detector_from_dict(d):
    ''' Convert the dictionary to a detector model

    Params:
        d The dictionary of parameters

    Returns:
        The detector model

    '''
    from dials.model.experiment import Detector, PanelList
    return Detector(PanelList([panel_from_dict(p) for p in d]))
