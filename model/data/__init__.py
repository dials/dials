from __future__ import division
from cctbx.array_family import flex
import boost.python
ext = boost.python.import_ext("dials_model_data_ext")
from dials_model_data_ext import *

def getattrlist(self, name):
  return [getattr(r, name) for r in self]

def getxyzcalmm(self):
  return [r.image_coord_mm + (r.rotation_angle,) for r in self]

def getxyzcalpx(self):
  return [r.image_coord_px + (r.frame_number,) for r in self]

def reflection_list_to_table(self):
  ''' Convert a reflection list to a table. '''
  from dials.array_family import flex

  # Create the table with the length of self
  table = flex.reflection_table(len(self))

  # General properties
  table['flags'] = flex.size_t(getattrlist(self, 'status'))
  table['id']    = flex.size_t(getattrlist(self, 'crystal'))
  table['panel'] = flex.size_t(getattrlist(self, 'panel_number'))

  # Predicted properties
  table['hkl']       = flex.miller_index(getattrlist(self, 'miller_index'))
  table['entering']  = flex.bool(getattrlist(self, 'entering'))
  table['s1']        = flex.vec3_double(getattrlist(self, 'beam_vector'))
  table['xyzcal.mm'] = flex.vec3_double(getxyzcalmm(self))
  table['xyzcal.px'] = flex.vec3_double(getxyzcalpx(self))

  # Observed centroid properties
  table['xyzobs.px.value']    = flex.vec3_double(
    getattrlist(self, 'centroid_position'))
  table['xyzobs.px.variance'] = flex.vec3_double(
    getattrlist(self, 'centroid_variance'))

  # Observed intensity properties
  table['intensity.raw.value']    = flex.double(
    getattrlist(self, 'intensity'))
  table['intensity.raw.variance'] = flex.double(
    getattrlist(self, 'intensity_variance'))
  table['intensity.cor.value']    = flex.double(
    getattrlist(self, 'corrected_intensity'))
  table['intensity.cor.variance'] = flex.double(
    getattrlist(self, 'corrected_intensity_variance'))

  # Shoebox properties
  table['shoebox.bbox'] = flex.int6(getattrlist(self, 'bounding_box'))

  # Return the table
  return table

ReflectionList.to_table = reflection_list_to_table
