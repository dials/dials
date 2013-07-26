from __future__ import division
from cctbx.array_family import flex # import dependency
from dials_model_data_ext import *

#reflection_list_attributes = [
#  'miller_index',
#  'entering',
#  'status',
#  'rotation_angle',
#  'beam_vector',
#  'image_coord_px',
#  'image_coord_mm',
#  'frame_number',
#  'panel_number',
#  'bounding_box',
#  'centroid_position',
#  'centroid_variance',
#  'centroid_sq_width',
#  'intensity',
#  'intensity_variance',
#]

#reflection_list_profiles = [
#  'shoebox',
#  'shoebox_mask',
#  'shoebox_background',
#  'transformed_shoebox',
#]


#def reflection_list_pickler(self):
#    state = {}
#    state['size'] = len(self)
#    for attr in reflection_list_attributes + reflection_list_profiles:
#        state[attr] = getattr(self, attr)()
#    return state

#def reflection_list_unpickler(self, state):
#    self.resize(state['size'])
#    for attr in reflection_list_attributes:
#        try:
#            print attr, len(state[attr]), len(self)
#            getattr(self, attr)(state[attr])
#        except KeyError:
#            pass

#    self.shoebox(state['shoebox'])
#    self.shoebox_mask(state['shoebox_mask'])
#    #self.shoebox_background(state['shoebox_background'])
#    self.transformed_shoebox(state['transformed_shoebox'])

#def get_list(rl, attr):
#    try:
#        l = getattr(rl, attr)()
#    except AttributeError:
#        l = [None] * len(rl)
#        for i, r in enumerate(rl):
#            l[i] = getattr(r, attr)
#    return l

#def set_list(rl, state, attr):
#    assert(len(rl) == len(state[attr]))
#    try:
#        getattr(rl, attr)(state[attr])
#    except AttributeError:
#        for r, v in zip(rl, state[attr]):
#            setattr(r, attr, v)
#

#def reflection_list_pickler(self):
#    state = {}
#    state['size'] = len(self)
#    for attr in reflection_list_attributes:
#        state[attr] = get_list(self, attr)
#    return state

#def reflection_list_unpickler(self, state):
#    self.resize(state['size'])
#    for attr in reflection_list_attributes:
#        try:
#            set_list(self, state, attr)
#        except KeyError:
#            pass


# Inject the __getinitargs__ method into ReflectionList
#ReflectionList.__getstate__ = reflection_list_pickler
#ReflectionList.__setstate__ = reflection_list_unpickler
