from __future__ import division
from cctbx.array_family import flex
import boost.python
ext = boost.python.import_ext("dials_model_data_ext")
from dials_model_data_ext import *

def reflection_list_to_table(self):
  ''' Convert a reflection list to a table. '''
  from dials.array_family import flex

  # Create the table with the length of self
  table = flex.reflection_table(len(self))

  # General properties
  table['flags'] = flex.int(len(self))
  table['id']    = flex.int(len(self))
  table['panel'] = flex.int(len(self))

  # Predicted properties
  table['hkl']       = flex.int(len(self))
  table['entering']  = flex.int(len(self))
  table['phi']       = flex.int(len(self))
  table['s1']        = flex.int(len(self))
  table['xyzcal.mm'] = flex.int(len(self))
  table['xyzcal.px'] = flex.int(len(self))

  # Observed centroid properties
  table['xyzobs.px.value']    = flex.int(len(self))
  table['xyzobs.px.variance'] = flex.int(len(self))
  table['xyzobs.mm.value']    = flex.int(len(self))
  table['xyzobs.mm.variance'] = flex.int(len(self))

  # Observed intensity properties
  table['intensity.raw.value']    = flex.int(len(self))
  table['intensity.raw.variance'] = flex.int(len(self))
  table['intensity.cor.value']    = flex.int(len(self))
  table['intensity.cor.variance'] = flex.int(len(self))

  # Shoebox properties
  table['shoebox.bbox'] = flex.int(len(self))
  table['shoebox.data'] = flex.int(len(self))
  table['shoebox.bgrd'] = flex.int(len(self))
  table['shoebox.mask'] = flex.int(len(self))

  # Return the table
  return table

ReflectionList.to_table = reflection_list_to_table
