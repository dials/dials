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

def getshoebox(self):
  from dials.model.data import Shoebox
  result = []
  for r in self:
    sbox = Shoebox(r.bounding_box)
    if isinstance(sbox.data, flex.float):
      sbox.data = r.shoebox.as_float()
      sbox.mask = r.shoebox_mask
      sbox.background = r.shoebox_background.as_float()
    else:
      sbox.data = r.shoebox
      sbox.mask = r.shoebox_mask
      sbox.background = r.shoebox_background
    result.append(sbox)
  return result

def setxyzcalmm(self, data):
  assert(len(self) == len(data))
  for r, d in zip(self, data):
    r.image_coord_mm = d[0:2]
    r.rotation_angle = d[2]

def setxyzcalpx(self, data):
  assert(len(self) == len(data))
  for r, d in zip(self, data):
    r.image_coord_px = d[0:2]
    r.frame_number = d[2]

def setshoebox(self, data):
  assert(len(self) == len(data))
  for r, d in zip(self, data):
    if isinstance(d.data, flex.float):
      r.shoebox = d.data.as_double()
      r.shoebox_mask = d.mask
      r.shoebox_background = d.background.as_double()
    else:
      r.shoebox = d.data
      r.shoebox_mask = d.mask
      r.shoebox_background = d.background

_warning_string = '''

=========================================================
The Reflection and ReflectionList classes are deprecated!
  Instances should be removed as soon as possible.
=========================================================

'''
class injector1(object):
    class __metaclass__(ReflectionList.__class__):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)

            old_init = ReflectionList.__init__
            def new_init(self, *args, **kwargs):
              import warnings
              import inspect
              _called = 'Called from: %s (%s:%d)' % \
                        (inspect.stack()[2][3], inspect.stack()[2][1],
                         inspect.stack()[2][2])
              warnings.simplefilter('always', DeprecationWarning)
              warnings.warn(_warning_string + _called, DeprecationWarning)

              old_init(self, *args, **kwargs)
            ReflectionList.__init__ = new_init
            return type.__init__(self, name, bases, dict)
class injector2(object):
    class __metaclass__(Reflection.__class__):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)

            old_init = Reflection.__init__
            def new_init(self, *args, **kwargs):
              import warnings
              import inspect
              _called = 'Called from: "%s" (%s:%d)' % \
                        (inspect.stack()[2][3], inspect.stack()[2][1],
                         inspect.stack()[2][2])
              warnings.simplefilter('always', DeprecationWarning)
              warnings.warn(_warning_string + _called, DeprecationWarning)
              old_init(self, *args, **kwargs)
            Reflection.__init__ = new_init
            return type.__init__(self, name, bases, dict)

class ReflectionAux(injector2, Reflection):
  pass

class ReflectionListAux(injector1, ReflectionList):

  def to_table(self, centroid_is_mm=False):
    ''' Convert a reflection list to a table. '''
    from dials.array_family import flex

    # Create the table with the length of self
    table = flex.reflection_table(len(self))

    # General properties
    table['flags'] = flex.size_t(getattrlist(self, 'status'))
    table['id']    = flex.size_t(getattrlist(self, 'crystal'))
    table['panel'] = flex.size_t(getattrlist(self, 'panel_number'))

    # Predicted properties
    table['miller_index'] = flex.miller_index(getattrlist(self, 'miller_index'))
    table['entering']     = flex.bool(getattrlist(self, 'entering'))
    table['s1']           = flex.vec3_double(getattrlist(self, 'beam_vector'))
    table['xyzcal.mm']    = flex.vec3_double(getxyzcalmm(self))
    table['xyzcal.px']    = flex.vec3_double(getxyzcalpx(self))

    # Observed centroid properties
    if centroid_is_mm:
      table['xyzobs.mm.value']    = flex.vec3_double(
        getattrlist(self, 'centroid_position'))
      table['xyzobs.mm.variance'] = flex.vec3_double(
        getattrlist(self, 'centroid_variance'))
    else:
      table['xyzobs.px.value']    = flex.vec3_double(
        getattrlist(self, 'centroid_position'))
      table['xyzobs.px.variance'] = flex.vec3_double(
        getattrlist(self, 'centroid_variance'))

    # Observed intensity properties
    table['intensity.sum.value']    = flex.double(
      getattrlist(self, 'intensity'))
    table['intensity.sum.variance'] = flex.double(
      getattrlist(self, 'intensity_variance'))

    # Shoebox properties
    table['bbox'] = flex.int6(getattrlist(self, 'bounding_box'))
    table['shoebox'] = flex.shoebox(getshoebox(self))

    # Return the table
    return table

  @staticmethod
  def from_table(table):
    rlist = ReflectionList(table.nrows())
    if 'flags' in table:
      setattrlist(rlist, 'status', table['flags'])
    if 'id' in table:
      setattrlist(rlist, 'crystal', table['id'])
    if 'panel' in table:
      setattrlist(rlist, 'panel_number', table['panel'])
    if 'miller_index' in table:
      setattrlist(rlist, 'miller_index', table['miller_index'])
    if 'entering' in table:
      setattrlist(rlist, 'entering', table['entering'])
    if 's1' in table:
      setattrlist(rlist, 'beam_vector', table['s1'])
    if 'xyzobs.px.value' in table:
      setattrlist(rlist, 'centroid_position', table['xyzobs.px.value'])
    if 'xyzobs.px.variance' in table:
      setattrlist(rlist, 'centroid_variance', table['xyzobs.px.variance'])
    if 'intensity.sum.value' in table:
      setattrlist(rlist, 'intensity', table['intensity.sum.value'])
    if 'intensity.sum.variance' in table:
      setattrlist(rlist, 'intensity_variance', table['intensity.sum.variance'])
    if 'bbox' in table:
      setattrlist(rlist, 'bounding_box', table['bbox'])
    if 'shoebox' in table:
      setshoebox(rlist, table['shoebox'])
    if 'xyzcal.px' in table:
      setxyzcalpx(rlist, table['xyzcal.px'])
    if 'xyzcal.mm' in table:
      setxyzcalmm(rlist, table['xyzcal.mm'])
    return rlist

# _Reflection = Reflection
# _ReflectionList = ReflectionList

# _warning_string = '''

# =========================================================
# The Reflection and ReflectionList classes are deprecated!
#   Instances should be removed as soon as possible.
# =========================================================

# '''
# import functools

# @functools.wraps(_Reflection)
# def Reflection(*args, **kwargs):
#     import warnings
#     warnings.simplefilter('always', DeprecationWarning)
#     warnings.warn(_warning_string, DeprecationWarning)
#     return _Reflection(*args, **kwargs)

# @functools.wraps(ReflectionList)
# def ReflectionList(*args, **kwargs):
#     import warnings
#     warnings.simplefilter('always', DeprecationWarning)
#     warnings.warn(_warning_string, DeprecationWarning)
#     return _ReflectionList(*args, **kwargs)
