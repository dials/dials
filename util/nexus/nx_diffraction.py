
from __future__ import division

schema_url = 'https://github.com/nexusformat/definitions/blob/master/contributed_definitions/NXdiffraction.nxdl.xml'

def make_dataset(handle, name, dtype, data, description):
  dset = handle.create_dataset(
    name,
    (len(data),),
    dtype=dtype,
    data=data.as_numpy_array().astype(dtype))
  dset.attrs['description'] = description
  return dset

def make_uint(handle, name, data, description):
  return make_dataset(handle, name, "uint64", data, description)

def make_int(handle, name, data, description):
  return make_dataset(handle, name, "int64", data, description)

def make_bool(handle, name, data, description):
  return make_dataset(handle, name, "bool", data, description)

def make_float(handle, name, data, description):
  return make_dataset(handle, name, "float64", data, description)

def make_vlen_uint(handle, name, data, description):
  import h5py
  import numpy as np
  dtype = h5py.special_dtype(vlen=np.dtype("uint64"))
  dset = handle.create_dataset(name, (len(data),), dtype=dtype)
  for i, d in enumerate(data):
    if len(d) > 0:
      dset[i] = d
  dset.attrs['description'] = description
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
    make_uint(handle, "id", col, dsc)
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
    make_float(handle, 'prd_px_x',  col1, dsc1)
    make_float(handle, 'prd_px_y',  col2, dsc2)
    make_float(handle, 'prd_frame', col3, dsc3)
  elif key == 'xyzcal.mm':
    col1, col2, col3 = data.parts()
    dsc1 = 'The predicted bragg peak fast millimeter location'
    dsc2 = 'The predicted bragg peak slow millimeter location'
    dsc3 = 'The predicted bragg peak rotation angle number'
    make_float(handle, 'prd_mm_x', col1, dsc1)
    make_float(handle, 'prd_mm_y', col2, dsc2)
    make_float(handle, 'prd_phi',  col3, dsc3)
  elif key == 'bbox':
    col1, col2, col3, col4, col5, col6 = data.parts()
    dsc1 = 'The bounding box lower x bound'
    dsc2 = 'The bounding box upper x bound'
    dsc3 = 'The bounding box lower y bound'
    dsc4 = 'The bounding box upper y bound'
    dsc5 = 'The bounding box lower z bound'
    dsc6 = 'The bounding box upper z bound'
    make_int(handle, 'bbx0', col1, dsc1)
    make_int(handle, 'bbx1', col2, dsc2)
    make_int(handle, 'bby0', col3, dsc3)
    make_int(handle, 'bby1', col4, dsc4)
    make_int(handle, 'bbz0', col5, dsc5)
    make_int(handle, 'bbz1', col6, dsc6)
  elif key == 'xyzobs.px.value':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel value'
    dsc2 = 'The observed centroid slow pixel value'
    dsc3 = 'The observed centroid frame value'
    make_float(handle, 'obs_px_x_val',  col1, dsc1)
    make_float(handle, 'obs_px_y_val',  col2, dsc2)
    make_float(handle, 'obs_frame_val', col3, dsc3)
  elif key == 'xyzobs.px.variance':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel variance'
    dsc2 = 'The observed centroid slow pixel variance'
    dsc3 = 'The observed centroid frame variance'
    make_float(handle, 'obs_px_x_var',  col1, dsc1)
    make_float(handle, 'obs_px_y_var',  col2, dsc2)
    make_float(handle, 'obs_frame_var', col3, dsc3)
  elif key == 'xyzobs.mm.value':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel value'
    dsc2 = 'The observed centroid slow pixel value'
    dsc3 = 'The observed centroid phi value'
    make_float(handle, 'obs_mm_x_val', col1, dsc1)
    make_float(handle, 'obs_mm_y_val', col2, dsc2)
    make_float(handle, 'obs_phi_val',  col3, dsc3)
  elif key == 'xyzobs.mm.variance':
    col1, col2, col3 = data.parts()
    dsc1 = 'The observed centroid fast pixel variance'
    dsc2 = 'The observed centroid slow pixel variance'
    dsc3 = 'The observed centroid phi variance'
    make_float(handle, 'obs_mm_x_var',  col1, dsc1)
    make_float(handle, 'obs_mm_y_var',  col2, dsc2)
    make_float(handle, 'obs_phi_var', col3, dsc3)
  elif key == 'background.mean':
    col = data
    dsc = 'The mean background value'
    make_float(handle, 'bkg_mean', col, dsc)
  elif key == 'intensity.sum.value':
    col = data
    dsc = 'The value of the summed intensity'
    make_float(handle, 'int_sum_val', col, dsc)
  elif key == 'intensity.sum.variance':
    col = data
    dsc = 'The variance of the summed intensity'
    make_float(handle, 'int_sum_var', col, dsc)
  elif key == 'intensity.prf.value':
    col = data
    dsc = 'The value of the profile fitted intensity'
    make_float(handle, 'int_prf_val', col, dsc)
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
    return flex.size_t(handle['id'][:].astype(int))
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
    x = flex.double(handle['prd_px_x'][:])
    y = flex.double(handle['prd_px_y'][:])
    z = flex.double(handle['prd_frame'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzcal.mm':
    x = flex.double(handle['prd_mm_x'][:])
    y = flex.double(handle['prd_mm_y'][:])
    z = flex.double(handle['prd_phi'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'bbox':
    x0 = flex.int(handle['bbx0'][:].astype(np.int32))
    x1 = flex.int(handle['bbx1'][:].astype(np.int32))
    y0 = flex.int(handle['bby0'][:].astype(np.int32))
    y1 = flex.int(handle['bby1'][:].astype(np.int32))
    z0 = flex.int(handle['bbz0'][:].astype(np.int32))
    z1 = flex.int(handle['bbz1'][:].astype(np.int32))
    return flex.int6(x0, x1, y0, y1, z0, z1)
  elif key == 'xyzobs.px.value':
    x = flex.double(handle['obs_px_x_val'][:])
    y = flex.double(handle['obs_px_y_val'][:])
    z = flex.double(handle['obs_frame_val'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzobs.px.variance':
    x = flex.double(handle['obs_px_x_var'][:])
    y = flex.double(handle['obs_px_y_var'][:])
    z = flex.double(handle['obs_frame_var'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzobs.mm.value':
    x = flex.double(handle['obs_mm_x_val'][:])
    y = flex.double(handle['obs_mm_y_val'][:])
    z = flex.double(handle['obs_phi_val'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'xyzobs.mm.variance':
    x = flex.double(handle['obs_mm_x_var'][:])
    y = flex.double(handle['obs_mm_y_var'][:])
    z = flex.double(handle['obs_phi_var'][:])
    return flex.vec3_double(x, y, z)
  elif key == 'background.mean':
    return flex.double(handle['bkg_mean'][:])
  elif key == 'intensity.sum.value':
    return flex.double(handle['int_sum_val'][:])
  elif key == 'intensity.sum.variance':
    return flex.double(handle['int_sum_var'][:])
  elif key == 'intensity.prf.value':
    return flex.double(handle['int_prf_val'][:])
  elif key == 'intensity.prf.variance':
    return flex.double(handle['int_prf_var'][:])
  elif key == 'profile.correlation':
    return flex.double(handle['prf_cc'][:])
  elif key == 'lp':
    return flex.double(handle['lp'][:])
  else:
    raise KeyError('Column %s not read from file' % key)

def dump(entry, reflections, experiments):
  from dials.array_family import flex

  print "Dumping NXdiffraction"

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
  assert("diffraction" not in entry)
  diffraction = entry.create_group("diffraction")
  diffraction.attrs['NX_class'] = 'NXsubentry'

  # Create the definition
  definition = diffraction.create_dataset('definition', data='NXdiffraction')
  definition.attrs['version'] = 1
  definition.attrs['URL'] = schema_url

  diffraction['experiments'] = experiments

  if reflections is None:
    return

  # For each column in the reflection table dump to file
  for key, data in reflections.cols():
    try:
      write(diffraction, key, data)
    except KeyError, e:
      print e

  # FIXME Write the overlaps (for testing at the moment)
  overlaps = [[] for i in range(len(reflections))]
  overlaps[0] = [1, 2, 3]
  overlaps[1] = [0, 4]
  overlaps[2] = [0, 3]
  overlaps[3] = [0, 2]
  overlaps[4] = [1]
  make_vlen_uint(diffraction, "overlaps", overlaps, 'Reflection overlap list')

def load(entry):
  from dials.array_family import flex

  print "Loading NXdiffraction"

  # Check the feature is present
  assert("features" in entry)
  assert(7 in entry["features"])

  # Get the entry
  diffraction = entry['diffraction']
  assert(diffraction.attrs['NX_class'] == 'NXsubentry')

  # Get the definition
  definition = diffraction['definition']
  assert(definition.value == 'NXdiffraction')
  assert(definition.attrs['version'] == 1)

  # The paths to the experiments
  experiments = list(diffraction['experiments'])

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
    'lp'
  ]

  # The reflection table
  table = None

  # For each column in the reflection table dump to file
  for key in columns:
    try:
      col = read(diffraction, key)
      if table is None:
        table = flex.reflection_table()
      table[key] = col
    except KeyError, e:
      pass

  # Return the table
  return table, experiments
