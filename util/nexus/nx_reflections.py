
from __future__ import absolute_import, division

schema_url = 'https://github.com/nexusformat/definitions/blob/master/contributed_definitions/NXreflections.nxdl.xml'

def make_dataset(handle, name, dtype, data, description, units=None):
  dset = handle.create_dataset(
    name,
    data.focus(),
    dtype=dtype,
    data=data.as_numpy_array().astype(dtype))
  dset.attrs['description'] = description
  if units is not None:
    dset.attrs['units'] = units
  return dset

def make_uint(handle, name, data, description, units=None):
  return make_dataset(handle, name, "uint64", data, description, units)

def make_int(handle, name, data, description, units=None):
  return make_dataset(handle, name, "int64", data, description, units)

def make_bool(handle, name, data, description, units=None):
  return make_dataset(handle, name, "bool", data, description, units)

def make_float(handle, name, data, description, units=None):
  return make_dataset(handle, name, "float64", data, description, units)

def make_vlen_uint(handle, name, data, description, units=None):
  import h5py
  import numpy as np
  dtype = h5py.special_dtype(vlen=np.dtype("uint64"))
  dset = handle.create_dataset(name, (len(data),), dtype=dtype)
  for i, d in enumerate(data):
    if len(d) > 0:
      dset[i] = d
  dset.attrs['description'] = description
  if units is not None:
    dset.attrs['units'] = units
  return dset

def write(handle, key, data):
  from dials.array_family import flex
  if   key == 'miller_index':
    col1, col2, col3 = zip(*list(data))
    col1 = flex.int(col1)
    col2 = flex.int(col2)
    col3 = flex.int(col3)
    dsc1 = 'The h component of the miller index'
    dsc2 = 'The k component of the miller index'
    dsc3 = 'The l component of the miller index'
    make_int(handle, "h", col1, dsc1)
    make_int(handle, "k", col2, dsc2)
    make_int(handle, "l", col3, dsc3)
  elif key == 'id':
    col = data
    dsc = 'The experiment id'
    make_int(handle, "id", col, dsc)
  elif key == 'partial_id':
    col = data
    desc = 'The reflection id'
    make_uint(handle, "reflection_id", col, desc)
  elif key == 'entering':
    col = data
    dsc = 'Entering or exiting the Ewald sphere'
    make_bool(handle, 'entering', col, dsc)
  elif key == 'flags':
    col = data
    dsc = 'Status of the reflection in processing'
    make_uint(handle, 'flags', col, dsc)
  elif key == 'panel':
    col = data
    dsc = 'The detector module on which the reflection was recorded'
    make_uint(handle, 'det_module', col, dsc)
  elif key == 'd':
    col = data
    dsc = 'The resolution of the reflection'
    make_float(handle, 'd', col, dsc)
  elif key == 'partiality':
    col = data
    dsc = 'The partiality of the reflection'
    make_float(handle, 'partiality', col, dsc)
  elif key == 'xyzcal.px':
    col1, col2, col3 = data.parts()
    dsc1 = 'The predicted bragg peak fast pixel location'
    dsc2 = 'The predicted bragg peak slow pixel location'
    dsc3 = 'The predicted bragg peak frame number'
    make_float(handle, 'predicted_px_x',  col1, dsc1)
    make_float(handle, 'predicted_px_y',  col2, dsc2)
    make_float(handle, 'predicted_frame', col3, dsc3)
  elif key == 'xyzcal.mm':
    col1, col2, col3 = data.parts()
    dsc1 = 'The predicted bragg peak fast millimeter location'
    dsc2 = 'The predicted bragg peak slow millimeter location'
    dsc3 = 'The predicted bragg peak rotation angle number'
    make_float(handle, 'predicted_x', col1, dsc1, units='mm')
    make_float(handle, 'predicted_y', col2, dsc2, units='mm')
    make_float(handle, 'predicted_phi',  col3, dsc3, units='rad')
  elif key == 'bbox':
    d = data.as_int()
    d.reshape(flex.grid((len(data)),6))
    make_int(handle, 'bounding_box', d, 'The reflection bounding box')
  elif key == 'xyzobs.px.value':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel value'
    dsc2 = 'The observed centroid slow pixel value'
    dsc3 = 'The observed centroid frame value'
    make_float(handle, 'observed_px_x',  col1, dsc1)
    make_float(handle, 'observed_px_y',  col2, dsc2)
    make_float(handle, 'observed_frame', col3, dsc3)
  elif key == 'xyzobs.px.variance':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel variance'
    dsc2 = 'The observed centroid slow pixel variance'
    dsc3 = 'The observed centroid frame variance'
    make_float(handle, 'observed_px_x_var',  col1, dsc1)
    make_float(handle, 'observed_px_y_var',  col2, dsc2)
    make_float(handle, 'observed_frame_var', col3, dsc3)
  elif key == 'xyzobs.mm.value':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel value'
    dsc2 = 'The observed centroid slow pixel value'
    dsc3 = 'The observed centroid phi value'
    make_float(handle, 'observed_x', col1, dsc1, units='mm')
    make_float(handle, 'observed_y', col2, dsc2, units='mm')
    make_float(handle, 'observed_phi',  col3, dsc3, units='rad')
  elif key == 'xyzobs.mm.variance':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel variance'
    dsc2 = 'The observed centroid slow pixel variance'
    dsc3 = 'The observed centroid phi variance'
    make_float(handle, 'observed_x_var',  col1, dsc1, units='mm')
    make_float(handle, 'observed_y_var',  col2, dsc2, units='mm')
    make_float(handle, 'observed_phi_var', col3, dsc3, units='rad')
  elif key == 'background.mean':
    col = data
    dsc = 'The mean background value'
    make_float(handle, 'background_mean', col, dsc)
  elif key == 'intensity.sum.value':
    col = data
    dsc = 'The value of the summed intensity'
    make_float(handle, 'int_sum', col, dsc)
  elif key == 'intensity.sum.variance':
    col = data
    dsc = 'The variance of the summed intensity'
    make_float(handle, 'int_sum_var', col, dsc)
  elif key == 'intensity.prf.value':
    col = data
    dsc = 'The value of the profile fitted intensity'
    make_float(handle, 'int_prf', col, dsc)
  elif key == 'intensity.prf.variance':
    col = data
    dsc = 'The variance of the profile fitted intensity'
    make_float(handle, 'int_prf_var', col, dsc)
  elif key == 'profile.correlation':
    col = data
    dsc = 'Profile fitting correlations'
    make_float(handle, 'prf_cc', col, dsc)
  elif key == 'lp':
    col = data
    dsc = 'The lorentz-polarization correction factor'
    make_float(handle, 'lp', col, dsc)
  elif key == 'num_pixels.background':
    col = data
    dsc = 'Number of background pixels'
    make_int(handle, 'num_bg', col, dsc)
  elif key == 'num_pixels.background_used':
    col = data
    dsc = 'Number of background pixels used'
    make_int(handle, 'num_bg_used', col, dsc)
  elif key == 'num_pixels.foreground':
    col = data
    dsc = 'Number of foreground pixels'
    make_int(handle, 'num_fg', col, dsc)
  elif key == 'num_pixels.valid':
    col = data
    dsc = 'Number of valid pixels'
    make_int(handle, 'num_valid', col, dsc)
  elif key == 'profile.rmsd':
    col = data
    dsc = 'Profile rmsd'
    make_float(handle, 'prf_rmsd', col, dsc)
  else:
    raise KeyError('Column %s not written to file' % key)

def read(handle, key):
  from dials.array_family import flex
  import numpy as np
  if   key == 'miller_index':
    h = flex.int(handle['h'][:].astype(np.int32))
    k = flex.int(handle['k'][:].astype(np.int32))
    l = flex.int(handle['l'][:].astype(np.int32))
    return flex.miller_index(h,k,l)
  elif key == 'id':
    return flex.int(handle['id'][:].astype(int))
  elif key == 'partial_id':
    return flex.size_t(handle['reflection_id'][:].astype(int))
  elif key == 'entering':
    return flex.bool(handle['entering'][:])
  elif key == 'flags':
    return flex.size_t(handle['flags'][:].astype(int))
  elif key == 'panel':
    return flex.size_t(handle['det_module'][:].astype(int))
  elif key == 'd':
    return flex.double(handle['d'][:])
  elif key == 'partiality':
    return flex.double(handle['partiality'][:])
  elif key == 'xyzcal.px':
    x = flex.double(handle['predicted_px_x'][:])
    y = flex.double(handle['predicted_px_y'][:])
    z = flex.double(handle['predicted_frame'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzcal.mm':
    x = flex.double(handle['predicted_x'][:])
    y = flex.double(handle['predicted_y'][:])
    z = flex.double(handle['predicted_phi'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'bbox':
    b = flex.int(handle['bounding_box'][:].astype(np.int32))
    return flex.int6(b.as_1d())
  elif key == 'xyzobs.px.value':
    x = flex.double(handle['observed_px_x'][:])
    y = flex.double(handle['observed_px_y'][:])
    z = flex.double(handle['observed_frame'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzobs.px.variance':
    x = flex.double(handle['observed_px_x_var'][:])
    y = flex.double(handle['observed_px_y_var'][:])
    z = flex.double(handle['observed_frame_var'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzobs.mm.value':
    x = flex.double(handle['observed_x'][:])
    y = flex.double(handle['observed_y'][:])
    z = flex.double(handle['observed_phi'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzobs.mm.variance':
    x = flex.double(handle['observed_x_var'][:])
    y = flex.double(handle['observed_y_var'][:])
    z = flex.double(handle['observed_phi_var'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'background.mean':
    return flex.double(handle['background_mean'][:])
  elif key == 'intensity.sum.value':
    return flex.double(handle['int_sum'][:])
  elif key == 'intensity.sum.variance':
    return flex.double(handle['int_sum_var'][:])
  elif key == 'intensity.prf.value':
    return flex.double(handle['int_prf'][:])
  elif key == 'intensity.prf.variance':
    return flex.double(handle['int_prf_var'][:])
  elif key == 'profile.correlation':
    return flex.double(handle['prf_cc'][:])
  elif key == 'lp':
    return flex.double(handle['lp'][:])
  elif key == 'num_pixels.background':
    return flex.int(handle['num_bg'][:].astype(np.int32))
  elif key == 'num_pixels.background_used':
    return flex.int(handle['num_bg_used'][:].astype(np.int32))
  elif key == 'num_pixels.foreground':
    return flex.int(handle['num_fg'][:].astype(np.int32))
  elif key == 'num_pixels.valid':
    return flex.int(handle['num_valid'][:].astype(np.int32))
  elif key == 'profile.rmsd':
    return flex.double(handle['prf_rmsd'][:])
  else:
    raise KeyError('Column %s not read from file' % key)

def dump(entry, reflections, experiments):
  from dials.array_family import flex

  print "Dumping NXreflections"

  # Add the feature
  if "features" in entry:
    features = entry['features']
    assert(features.dtype == 'uint64')
    features.resize((len(features)+1,))
    features[len(features)-1] = 7
  else:
    import numpy as np
    features = entry.create_dataset(
      "features",
      (1,),
      maxshape=(None,),
      dtype=np.uint64)
    features[0] = 7

  # Create the entry
  assert("reflections" not in entry)
  refls = entry.create_group("reflections")
  refls.attrs['NX_class'] = 'NXsubentry'

  # Create the definition
  definition = refls.create_dataset('definition', data='NXreflections')
  definition.attrs['version'] = 1
  definition.attrs['URL'] = schema_url

  refls['experiments'] = experiments

  if reflections is None:
    return

  # For each column in the reflection table dump to file
  for key, data in reflections.cols():
    try:
      write(refls, key, data)
    except KeyError, e:
      print e

  # FIXME Write the overlaps (for testing at the moment)
  overlaps = [[] for i in range(len(reflections))]
  overlaps[0] = [1, 2, 3]
  overlaps[1] = [0, 4]
  overlaps[2] = [0, 3]
  overlaps[3] = [0, 2]
  overlaps[4] = [1]
  make_vlen_uint(refls, "overlaps", overlaps, 'Reflection overlap list')

def load(entry):
  from dials.array_family import flex

  print "Loading NXreflections"

  # Check the feature is present
  assert("features" in entry)
  assert(7 in entry["features"])

  # Get the entry
  refls = entry['reflections']
  assert(refls.attrs['NX_class'] == 'NXsubentry')

  # Get the definition
  definition = refls['definition']
  assert(definition.value == 'NXreflections')
  assert(definition.attrs['version'] == 1)

  # The paths to the experiments
  experiments = list(refls['experiments'])

  # The columns to try
  columns = [
    'miller_index',
    'id',
    'partial_id',
    'entering',
    'flags',
    'panel',
    'd',
    'partiality',
    'xyzcal.px',
    'xyzcal.mm',
    'bbox',
    'xyzobs.px.value',
    'xyzobs.px.variance',
    'xyzobs.mm.value',
    'xyzobs.mm.variance',
    'background.mean',
    'intensity.sum.value',
    'intensity.sum.variance',
    'intensity.prf.value',
    'intensity.prf.variance',
    'profile.correlation',
    'lp',
    'num_pixels.background',
    'num_pixels.foreground',
    'num_pixels.background_used',
    'num_pixels.valid',
    'profile.rmsd',
  ]

  # The reflection table
  table = None

  # For each column in the reflection table dump to file
  for key in columns:
    try:
      col = read(refls, key)
      if table is None:
        table = flex.reflection_table()
      table[key] = col
    except KeyError, e:
      pass

  # Return the table
  return table, experiments
