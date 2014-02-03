from __future__ import division
from cctbx.array_family import flex
import boost.python
ext = boost.python.import_ext("dials_model_data_ext")
from dials_model_data_ext import *

def getattrlist(self, name):
  return [getattr(r, name) for r in self]

def setattrlist(self, name, data):
  assert(len(self) == len(data))
  for r, d in zip(self, data):
    setattr(r, name, d)

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

@staticmethod
def reflection_list_from_table(table):
  rlist = ReflectionList(table.nrows())
  if 'flags' in table:
    setattrlist(rlist, 'status', table['flags'])
  if 'id' in table:
    setattrlist(rlist, 'crystal', table['id'])
  if 'panel' in table:
    setattrlist(rlist, 'panel_number', table['panel'])
  if 'hkl' in table:
    setattrlist(rlist, 'miller_index', table['hkl'])
  if 'entering' in table:
    setattrlist(rlist, 'entering', table['entering'])
  if 's1' in table:
    setattrlist(rlist, 'beam_vector', table['s1'])
  if 'xyzobs.px.value' in table:
    setattrlist(rlist, 'centroid_position', table['xyzobs.px.value'])
  if 'xyzobs.px.variance' in table:
    setattrlist(rlist, 'centroid_variance', table['xyzobs.px.variance'])
  if 'intensity.raw.value' in table:
    setattrlist(rlist, 'intensity', table['intensity.raw.value'])
  if 'intensity.raw.variance' in table:
    setattrlist(rlist, 'intensity_variance', table['intensity.raw.variance'])
  if 'intensity.cor.value' in table:
    setattrlist(rlist, 'corrected_intensity', table['intensity.cor.value'])
  if 'intensity.cor.variance' in table:
    setattrlist(rlist, 'corrected_intensity_variance',
                table['intensity.cor.variance'])
  if 'shoebox.bbox' in table:
    setattrlist(rlist, 'bounding_box', table['shoebox.bbox'])
  return rlist


ReflectionList.to_table = reflection_list_to_table
ReflectionList.from_table = reflection_list_from_table

