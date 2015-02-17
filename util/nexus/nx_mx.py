from __future__ import division

# FIXME - beam direction - Need to fix in dxtbx?

# Extensions to NXMX
#
#  "detector/underload"
#  "detector/timestamp"

schema_url = 'https://github.com/nexusformat/definitions/blob/master/applications/NXmx.nxdl.xml'

def get_nx_class(handle, klass, path):
  if path in handle:
    group = handle[path]
    assert(group.attrs['NX_class'] == klass)
  else:
    group = handle.create_group(path)
    group.attrs["NX_class"] = klass
  return group

def get_nx_sample(handle, path):
  return get_nx_class(handle, "NXsample", path)

def get_nx_beam(handle, path):
  return get_nx_class(handle, "NXbeam", path)

def get_nx_detector(handle, path):
  return get_nx_class(handle, "NXdetector", path)

def get_nx_detector_module(handle, path):
  return get_nx_class(handle, "NXdetector_module", path)

def get_nx_instrument(handle, path):
  return get_nx_class(handle, "NXinstrument", path)

def get_nx_transformations(handle, path):
  return get_nx_class(handle, "NXtransformations", path)

def get_nx_process(handle, path):
  return get_nx_class(handle, "NXprocess", path)

def get_nx_note(handle, path):
  return get_nx_class(handle, "NXnote", path)

def get_nx_data(handle, path, data):
  handle[path] = data
  return handle[path]

def dump_beam(entry, beam):
  ''' Export the beam model. '''
  from scitbx import matrix
  from math import sin, cos, acos, pi

  EPS = 1e-7

  # Make sure that the direction is 0, 0, -1
  assert(matrix.col(beam.get_direction()).dot(matrix.col((0, 0, -1))) < EPS)

  # Get the nx_beam
  nx_sample = get_nx_sample(entry, "sample")
  nx_beam = get_nx_beam(nx_sample, "beam")

  # Generate the stokes polarization parameters
  d = matrix.col(beam.get_direction())
  n = matrix.col(beam.get_polarization_normal())
  p = beam.get_polarization_fraction()

  assert(n.dot(d) < EPS)
  axis = d.cross(matrix.col((0, 1, 0)))
  I = 1.0
  X = 0.0
  W = acos(n.dot(axis) / n.length())
  S0 = I
  S1 = I*p*cos(2*W)*cos(2*X)
  S2 = I*p*sin(2*W)*cos(2*X)
  S3 = I*p*sin(2*X)

  # Set the beam parameters
  nx_beam['direction'] = d.normalize() # FIXME Non-standard
  nx_beam['incident_wavelength'] = beam.get_wavelength()
  nx_beam['incident_polarization_stokes'] = (S0, S1, S2, S3)

def dump_detector(entry, detector, beam, imageset, scan):
  from scitbx import matrix

  EPS = 1e-7

  # Get the detector
  nx_instrument = get_nx_instrument(entry, "instrument")
  nx_detector = get_nx_detector(nx_instrument, "detector")

  # Get some bulk properties
  thickness = detector[0].get_thickness()
  material = detector[0].get_material()
  dtype = detector[0].get_type()
  trusted_range = detector[0].get_trusted_range()

  # Check all panels obey these bulk properties
  assert([abs((p.get_thickness() - thickness) < EPS)
          for p in detector].count(False) == 0)
  assert([p.get_material() == material for p in detector].count(False) == 0)
  assert([p.get_type() == dtype for p in detector].count(False) == 0)
  assert([p.get_trusted_range() == trusted_range for p in detector].count(False) == 0)

  # Take the distance and beam centre from the first panel
  #distance = detector[0].get_distance()
  #beam_centre_x, beam_centre_y = detector[0].get_beam_centre(beam.get_s0())

  # Set the detector properties
  nx_detector['sensor_thickness'] = thickness
  nx_detector['sensor_material'] = material
  nx_detector['type'] = dtype
  #nx_detector['distance'] = distance
  #nx_detector['beam_centre_x'] = beam_centre_x
  #nx_detector['beam_centre_y'] = beam_centre_y
  nx_detector['saturation_value'] = int(trusted_range[1])
  nx_detector['underload'] = int(trusted_range[0]) # FIXME non-standard
  nx_detector['description'] = dtype

  # Make up some fake stuff
  if scan is not None:
    nx_detector['timestamp'] = scan.get_epochs().as_numpy_array() # FIXME non-standard
    nx_detector['dead_time'] = [0] * len(scan)
    nx_detector['count_time'] = [0] * len(scan)
    nx_detector['frame_time'] = scan.get_exposure_times().as_numpy_array()
    nx_detector['detector_readout_time'] = [0] * len(scan)
  else:
    nx_detector['timestamp'] = 0
    nx_detector['dead_time'] = 0
    nx_detector['count_time'] = 0
    nx_detector['frame_time'] = 0
    nx_detector['detector_readout_time'] = 0
  nx_detector['bit_depth_readout'] = 32

  # Create the detector depends on
  nx_detector['depends_on'] = '.'

  # Creat some data for example file
  #data = [im.as_numpy_array() for im in imageset]

  # Create the nx data
  #nx_data = get_nx_data(nx_detector, "data", [[[0]],[[0]],[[0]]])

  # Loop through all the panels
  for i, panel in enumerate(detector):

    # Get some panel attributes
    pixel_size = panel.get_pixel_size()
    image_size = panel.get_image_size()
    origin = matrix.col(panel.get_origin())

    # Get the detector module object
    nx_module = get_nx_detector_module(nx_detector, 'module%d' % i)

    # Set the data size
    nx_module['data_size'] = image_size
    nx_module['data_origin'] = (-1, -1) # FIXME INVALID

    # Set the module offset
    nx_module['module_offset'] = origin.length()
    nx_module['module_offset'].attrs['depends_on'] = '.'
    nx_module['module_offset'].attrs['transformation_type'] = 'translation'
    nx_module['module_offset'].attrs['offset'] = (0.0, 0.0, 0.0)
    nx_module['module_offset'].attrs['vector'] = origin.normalize()

    # The path for items below
    module_offset_path = nx_module['module_offset'].name

    # Write the fast pixel direction
    nx_module['fast_pixel_direction'] = pixel_size[0]
    nx_module['fast_pixel_direction'].attrs['depends_on'] = module_offset_path
    nx_module['fast_pixel_direction'].attrs['transformation_type'] = 'translation'
    nx_module['fast_pixel_direction'].attrs['offset'] = (0.0, 0.0, 0.0)
    nx_module['fast_pixel_direction'].attrs['vector'] = panel.get_fast_axis()

    # Write the slow pixel direction
    nx_module['slow_pixel_direction'] = pixel_size[1]
    nx_module['slow_pixel_direction'].attrs['depends_on'] = module_offset_path
    nx_module['slow_pixel_direction'].attrs['transformation_type'] = 'translation'
    nx_module['slow_pixel_direction'].attrs['offset'] = (0.0, 0.0, 0.0)
    nx_module['slow_pixel_direction'].attrs['vector'] = panel.get_slow_axis()

def dump_goniometer(entry, goniometer, scan):
  ''' Export the goniometer model. '''

  if scan is None or goniometer is None:
    return

  # The angles for each image
  phi0, dphi = scan.get_oscillation(deg=True)
  phi = [phi0+dphi*i for i in range(len(scan))]

  # Write out the rotation axis and oscillation
  nx_sample = get_nx_sample(entry, "sample")
  nx_transformations = get_nx_transformations(nx_sample, "transformations")
  nx_transformations['phi'] = phi
  nx_transformations['phi'].attrs['depends_on'] = '.'
  nx_transformations['phi'].attrs['transformation_type'] = 'rotation'
  nx_transformations['phi'].attrs['offset_units'] = 'mm'
  nx_transformations['phi'].attrs['offset'] = (0.0, 0.0, 0.0)
  nx_transformations['phi'].attrs['vector'] = goniometer.get_rotation_axis()

def dump_crystal(entry, crystal, scan):
  ''' Export the crystal model. '''
  from scitbx.array_family import flex

  # Get the sample
  nx_sample = get_nx_sample(entry, "sample")

  # Set the space group
  nx_sample['unit_cell_group'] = crystal.get_space_group().type().hall_symbol()

  # Get the unit cell and orientation matrix in the case of scan varying and
  # scan static models
  if crystal.num_scan_points:
    num = crystal.num_scan_points
    unit_cell = flex.double(flex.grid(num, 6))
    orientation_matrix = flex.double(flex.grid(num, 9))
    for i in range(num):
      __cell = crystal.get_unit_cell_at_scan_point(i).parameters()
      for j in range(6):
        unit_cell[i,j] = __cell[j]
      __matrix = crystal.get_U_at_scan_point(i)
      for j in range(9):
        orientation_matrix[i,j] = __matrix[j]
    orientation_matrix = [[tuple(orientation_matrix[i:i+1,0:3]),
                           tuple(orientation_matrix[i:i+1,3:6]),
                           tuple(orientation_matrix[i:i+1,6:9])] for i in range(num)]
    unit_cell = [tuple(unit_cell[i:i+1,:]) for i in range(num)]
    average_unit_cell = crystal.get_unit_cell().parameters()
    average_orientation_matrix = crystal.get_U()
    average_orientation_matrix = [
      average_orientation_matrix[0:3],
      average_orientation_matrix[3:6],
      average_orientation_matrix[6:9]]

  else:
    unit_cell = [crystal.get_unit_cell().parameters()]
    orientation_matrix = [crystal.get_U()]
    orientation_matrix = [[tuple(orientation_matrix[0][0:3]),
                           tuple(orientation_matrix[0][3:6]),
                           tuple(orientation_matrix[0][6:9])]]
    average_unit_cell = unit_cell[0]
    average_orientation_matrix = orientation_matrix[0]

  # Save the unit cell data
  nx_sample['name'] = "FROM_DIALS"
  nx_sample['unit_cell'] = unit_cell
  nx_sample['unit_cell'].attrs['angles_units'] = 'deg'
  nx_sample['unit_cell'].attrs['length_units'] = 'angstrom'

  # Save the orientation matrix
  nx_sample['orientation_matrix'] = orientation_matrix

  # Set an average unit cell etc for scan static stuff
  nx_sample['average_unit_cell'] = average_unit_cell
  nx_sample['average_unit_cell'].attrs['angles_units'] = 'deg'
  nx_sample['average_unit_cell'].attrs['length_units'] = 'angstrom'
  nx_sample['average_orientation_matrix'] = average_orientation_matrix

  # Set depends on
  if scan is not None:
    nx_sample['depends_on'] = nx_sample['transformations/phi'].name
  else:
    nx_sample['depends_on'] = '.'

def dump_details(entry):
  from time import strftime

  # Program info
  entry['program_name'] = 'dials.export_nxmx'
  entry['program_name'].attrs['version'] = 1
  entry['program_name'].attrs['configuration'] = ''

  # Set some processing information (each program should add itself)
  nx_process = get_nx_process(entry, "process")
  nx_process['program'] = 'dials'
  nx_process['version'] = 1
  nx_process['date'] = strftime('%Y-%m-%dT%H:%M:%S')

  nx_note = get_nx_note(nx_process, '0_spot_finding')
  nx_note['author'] = 'dials.find_spots'
  nx_note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  nx_note['type'] = 'text/plain'
  nx_note['description'] = 'Spot finding parameters'
  nx_note['data'] = 'dials.find_spots datablock.json'

  nx_note = get_nx_note(nx_process, '1_indexing')
  nx_note['author'] = 'dials.index'
  nx_note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  nx_note['type'] = 'text/plain'
  nx_note['description'] = 'Indexing parameters'
  nx_note['data'] = 'dials.index datablock.json strong.pickle'

  nx_note = get_nx_note(nx_process, '2_refinement')
  nx_note['author'] = 'dials.refine'
  nx_note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  nx_note['type'] = 'text/plain'
  nx_note['description'] = 'Refinement parameters'
  nx_note['data'] = 'dials.refine experiments.json indexed.pickle'

  nx_note = get_nx_note(nx_process, '3_integration')
  nx_note['author'] = 'dials.integrate'
  nx_note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  nx_note['type'] = 'text/plain'
  nx_note['description'] = 'Integration parameters'
  nx_note['data'] = 'dials.integrate refined_experiments.json indexed.pickle'

def load_beam(entry):
  from dxtbx.model import Beam
  from math import sqrt, atan, cos, sin, pi
  from scitbx import matrix

  EPS = 1e-7

  # Get the nx_beam
  nx_sample = get_nx_sample(entry, "sample")
  nx_beam = get_nx_beam(nx_sample, "beam")
  wavelength = nx_beam['incident_wavelength'].value
  direction = matrix.col(nx_beam['direction'])
  S0, S1, S2, S3 = tuple(nx_beam['incident_polarization_stokes'])
  I = S0
  p = sqrt(S1**2 + S2**2 + S3**2)/S0
  W = 0.5*atan(S2/S1)
  X = 0.5*atan(S3/sqrt(S1**2 + S2**2))
  W += pi
  assert(abs(I - 1.0) < EPS)
  assert(abs(X) < EPS)
  assert(p >= 0.0 and p <= 1.0)
  axis1 = direction.cross(matrix.col((0, 1, 0)))
  axis2 = direction.cross(axis1)
  n = (axis2 * cos(W) + axis1 * sin(W)).normalize()
  assert(abs(direction.dot(n)) < EPS)

  # Return the beam model
  return Beam(direction, wavelength, 0, 0, n, p)

def load_detector(entry):
  from dxtbx.model import Detector
  from scitbx import matrix

  # Get the detector module object
  nx_instrument = get_nx_instrument(entry, "instrument")
  nx_detector = get_nx_detector(nx_instrument, "detector")
  assert(nx_detector['depends_on'].value == '.')
  material = nx_detector['sensor_material'].value
  det_type = nx_detector['type'].value
  thickness = nx_detector['sensor_thickness'].value
  trusted_range = (nx_detector['underload'].value, nx_detector['saturation_value'].value)


  # The detector model
  detector = Detector()

  i = 0
  while True:
    try:
      module = get_nx_detector_module(nx_detector, "module%d" % i)
    except Exception:
      break
    # Set the data size
    image_size = module['data_size']

    # Set the module offset
    offset_length = module['module_offset'].value
    assert(module['module_offset'].attrs['depends_on'] == '.')
    assert(module['module_offset'].attrs['transformation_type'] == 'translation')
    assert(tuple(module['module_offset'].attrs['offset']) == (0, 0, 0))
    offset_vector = matrix.col(module['module_offset'].attrs['vector'])
    origin = offset_vector * offset_length

    # Write the fast pixel direction
    module_offset_path = module['module_offset'].name
    pixel_size_x = module['fast_pixel_direction'].value
    assert(module['fast_pixel_direction'].attrs['depends_on'] == module_offset_path)
    assert(module['fast_pixel_direction'].attrs['transformation_type'] == 'translation')
    assert(tuple(module['fast_pixel_direction'].attrs['offset']) == (0, 0, 0))
    fast_axis = tuple(module['fast_pixel_direction'].attrs['vector'])

    # Write the slow pixel direction
    pixel_size_y = module['slow_pixel_direction'].value
    assert(module['slow_pixel_direction'].attrs['depends_on'] == module_offset_path)
    assert(module['slow_pixel_direction'].attrs['transformation_type'] == 'translation')
    assert(tuple(module['slow_pixel_direction'].attrs['offset']) == (0, 0, 0))
    slow_axis = tuple(module['slow_pixel_direction'].attrs['vector'])

    # Get the pixel size and axis vectors
    pixel_size = (pixel_size_x, pixel_size_y)

    # Create the panel
    panel = detector.add_panel()
    panel.set_frame(fast_axis, slow_axis, origin)
    panel.set_pixel_size(pixel_size)
    panel.set_image_size(image_size)
    panel.set_type(det_type)
    panel.set_thickness(thickness)
    panel.set_material(material)
    panel.set_trusted_range(trusted_range)
    i += 1

  # Return the detector and panel
  return detector

def load_goniometer(entry):
  from dxtbx.model import Goniometer

  # Write out the rotation axis and oscillation
  nx_sample = get_nx_sample(entry, "sample")
  try:
    transformations = get_nx_transformations(nx_sample, "transformations")
  except Exception:
    return None
  assert(transformations['phi'].attrs['depends_on'] == '.')
  assert(transformations['phi'].attrs['transformation_type'] == 'rotation')
  assert(transformations['phi'].attrs['offset_units'] == 'mm')
  assert(tuple(transformations['phi'].attrs['offset']) == (0, 0, 0))
  rotation_axis = tuple(transformations['phi'].attrs['vector'])

  # Return the goniometer model
  return Goniometer(rotation_axis)

def load_scan(entry):
  from dxtbx.model import Scan

  # Write out the rotation axis and oscillation
  nx_sample = get_nx_sample(entry, "sample")
  try:
    transformations = get_nx_transformations(nx_sample, "transformations")
  except Exception:
    return None
  phi = transformations['phi']
  assert(transformations['phi'].attrs['depends_on'] == '.')
  assert(transformations['phi'].attrs['transformation_type'] == 'rotation')
  assert(transformations['phi'].attrs['offset_units'] == 'mm')
  assert(tuple(transformations['phi'].attrs['offset']) == (0, 0, 0))
  image_range = (1, len(phi))
  oscillation = (phi[0], phi[1] - phi[0])
  nx_instrument = get_nx_instrument(entry, "instrument")
  nx_detector = get_nx_detector(nx_instrument, "detector")
  exposure_time = nx_detector['frame_time']
  epochs = nx_detector['timestamp']
  return Scan(image_range, oscillation, exposure_time, epochs, deg=True)

def load_crystal(entry):
  from dxtbx.model.crystal import crystal_model as Crystal
  from scitbx.array_family import flex
  from scitbx import matrix
  from cctbx import uctbx
  import numpy

  # Get the sample
  nx_sample = get_nx_sample(entry, "sample")

  # Set the space group
  space_group_symbol = nx_sample['unit_cell_group'].value

  # Get depends on
  if nx_sample['depends_on'].value != '.':
    assert(nx_sample['depends_on'].value == nx_sample['transformations/phi'].name)

  # Read the average unit cell data
  average_unit_cell = flex.double(numpy.array(nx_sample['average_unit_cell']))
  assert(nx_sample['average_unit_cell'].attrs['angles_units'] == 'deg')
  assert(nx_sample['average_unit_cell'].attrs['length_units'] == 'angstrom')
  assert(len(average_unit_cell.all()) == 1)
  assert(len(average_unit_cell) == 6)
  average_orientation_matrix = flex.double(numpy.array(nx_sample['average_orientation_matrix']))
  assert(len(average_orientation_matrix.all()) == 2)
  assert(average_orientation_matrix.all()[0] == 3)
  assert(average_orientation_matrix.all()[1] == 3)

  # Get the real space vectors
  uc = uctbx.unit_cell(tuple(average_unit_cell))
  U = matrix.sqr(average_orientation_matrix)
  B = matrix.sqr(uc.fractionalization_matrix()).transpose()
  A = U * B
  A = A.inverse()
  real_space_a = A[0:3]
  real_space_b = A[3:6]
  real_space_c = A[6:9]

  # Read the unit cell data
  unit_cell = flex.double(numpy.array(nx_sample['unit_cell']))
  assert(nx_sample['unit_cell'].attrs['angles_units'] == 'deg')
  assert(nx_sample['unit_cell'].attrs['length_units'] == 'angstrom')

  # Read the orientation matrix
  orientation_matrix = flex.double(numpy.array(nx_sample['orientation_matrix']))
  assert(len(unit_cell.all()) == 2)
  assert(len(orientation_matrix.all()) == 3)
  assert(unit_cell.all()[0] == orientation_matrix.all()[0])
  assert(unit_cell.all()[1] == 6)
  assert(orientation_matrix.all()[1] == 3)
  assert(orientation_matrix.all()[2] == 3)

  # Construct the crystal model
  crystal = Crystal(
    real_space_a,
    real_space_b,
    real_space_c,
    space_group_symbol)

  # Sort out scan points
  if unit_cell.all()[0] > 1:
    A_list = []
    for i in range(unit_cell.all()[0]):
      uc = uctbx.unit_cell(tuple(unit_cell[i:i+1,:]))
      U = matrix.sqr(tuple(orientation_matrix[i:i+1,:,:]))
      B = matrix.sqr(uc.fractionalization_matrix()).transpose()
      A_list.append(U*B)
    crystal.set_A_at_scan_points(A_list)
  else:
    assert(unit_cell.all_eq(average_unit_cell))
    assert(orientation_matrix.all_eq(average_orientation_matrix))

  # Return the crystal
  return crystal

def dump(entry, experiments):
  from dials.array_family import flex

  # Add the feature
  if "features" in entry:
    features = entry['features']
    assert(features.dtype == 'uint64')
    features.resize((len(features)+1,))
    features[len(features)-1] = 6
  else:
    import numpy as np
    features = entry.create_dataset(
      "features",
      (1,),
      maxshape=(None,),
      dtype=np.uint64)
    features[0] = 6

  # Get the experiment
  for index, experiment in enumerate(experiments):

    # Create the entry
    assert(("experiment_%d" % index) not in entry)
    nxmx = entry.create_group("experiment_%d" % index)
    nxmx.attrs['NX_class'] = 'NXsubentry'
    nxmx.attrs['id'] = index

    # Create the definition
    definition = nxmx.create_dataset('definition', data='NXmx')
    definition.attrs['version'] = 1
    definition.attrs['URL'] = schema_url

    nxmx['title'] = 'FROM_DIALS'

    # Dump the models
    dump_beam(nxmx, experiment.beam)
    dump_detector(nxmx, experiment.detector, experiment.beam, experiment.imageset,
                  experiment.scan)
    dump_goniometer(nxmx, experiment.goniometer, experiment.scan)
    dump_crystal(nxmx, experiment.crystal, experiment.scan)

  # Dump some details
  dump_details(entry)

  # Link the data
  #nxmx['data'] = nxmx['instrument/detector/data']

def find_nx_mx_entries(nx_file, entry):
  '''
  Find NXmx entries

  '''
  hits = []
  def visitor(name, obj):
    if "NX_class" in obj.attrs.keys():
      if obj.attrs["NX_class"] in ["NXentry", "NXsubentry"]:
        if "definition" in obj.keys():
          if obj["definition"].value == "NXmx":
            hits.append(obj)
  nx_file[entry].visititems(visitor)
  return hits

def load(entry):
  from dxtbx.model.experiment.experiment_list import ExperimentList
  from dxtbx.model.experiment.experiment_list import Experiment

  # Check file contains the feature
  assert("features" in entry)
  assert(6 in entry['features'].value)

  experiment_list = ExperimentList()

  # Find all the experiments
  entries = find_nx_mx_entries(entry, ".")
  if len(entries) > 1:
    entries = sorted(entries, key=lambda x: x.attrs['id'])

  for nxmx in entries:

    # Get the definition
    definition = nxmx['definition']
    assert(definition.value == 'NXmx')
    assert(definition.attrs['version'] == 1)

    # Create the experiment
    experiment = Experiment()

    # Read the models
    experiment.beam = load_beam(nxmx)
    experiment.detector = load_detector(nxmx)
    experiment.goniometer = load_goniometer(nxmx)
    experiment.scan = load_scan(nxmx)
    experiment.crystal = load_crystal(nxmx)

    # Return the experiment list
    experiment_list.append(experiment)
  return experiment_list
