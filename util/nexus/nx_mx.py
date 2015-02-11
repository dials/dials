from __future__ import division

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
  from math import sin, cos

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
  I = 1.0
  X = 0.0
  W = n.dot(matrix.col((1, 0, 0))) / n.length()
  S0 = I
  S1 = I*p*cos(2*W)*cos(2*X)
  S2 = I*p*sin(2*W)*cos(2*X)
  S3 = I*p*sin(2*X)

  # Set the beam parameters
  nx_beam['incident_wavelength'] = beam.get_wavelength()
  nx_beam['incident_polarization_stokes'] = (S0, S1, S2, S3)

def dump_detector(entry, detector, beam, imageset):
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
  distance = detector[0].get_distance()
  beam_centre_x, beam_centre_y = detector[0].get_beam_centre(beam.get_s0())

  # Set the detector properties
  nx_detector['sensor_thickness'] = thickness
  nx_detector['sensor_material'] = material
  nx_detector['type'] = dtype
  nx_detector['distance'] = distance
  nx_detector['beam_centre_x'] = beam_centre_x
  nx_detector['beam_centre_y'] = beam_centre_y
  nx_detector['saturation_value'] = int(trusted_range[1]) # FIXME
  nx_detector['description'] = dtype

  # Make up some fake stuff
  nx_detector['count_time'] = 0.0
  nx_detector['dead_time'] = 0.0
  nx_detector['frame_time'] = 0.0
  nx_detector['detector_readout_time'] = 0.0
  nx_detector['bit_depth_readout'] = 32

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

    # Create the detector depends on
    nx_detector['depends_on'] = '.'

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

def dump_crystal(entry, crystal):
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
      __matrix = crystal.get_U_at_scan_point(i).transpose().elems
      for j in range(9):
        orientation_matrix[i,j] = __matrix[j]

  else:
    unit_cell = [crystal.get_unit_cell().parameters()]
    orientation_matrix = [crystal.get_U().transpose()]

  orientation_matrix = [[el[0:3],el[3:6],el[6:9]] for el in orientation_matrix]

  # Save the unit cell data
  nx_sample['name'] = "FROM_DIALS"
  nx_sample['unit_cell'] = unit_cell
  nx_sample['unit_cell'].attrs['angles_units'] = 'deg'
  nx_sample['unit_cell'].attrs['length_units'] = 'angstrom'

  # Save the orientation matrix
  nx_sample['orientation_matrix'] = orientation_matrix

  # Set depends on
  nx_sample['depends_on'] = nx_sample['transformations/phi'].name

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

def dump(entry, experiments):
  from dials.array_family import flex

  # Add the feature
  if "features" in entry:
    assert(entry['features'].dtype == 'uint64')
    entry['features'].append(6)
  else:
    import numpy as np
    features = entry.create_dataset("features", (1,), dtype=np.uint64)
    features[0] = 6

  # Create the entry
  assert(len(experiments) == 1)
  assert("experiment0" not in entry)
  nxmx = entry.create_group("experiment0")
  nxmx.attrs['NX_class'] = 'NXsubentry'

  # Create the definition
  definition = nxmx.create_dataset('definition', data='NXmx')
  definition.attrs['version'] = 1
  definition.attrs['URL'] = schema_url

  nxmx['title'] = 'FROM_DIALS'

  # Get the experiment
  experiment = experiments[0]

  # Dump the models
  dump_beam(nxmx, experiment.beam)
  dump_detector(nxmx, experiment.detector, experiment.beam, experiment.imageset)
  dump_goniometer(nxmx, experiment.goniometer, experiment.scan)
  dump_crystal(nxmx, experiment.crystal)

  # Dump some details
  dump_details(entry)

  # Link the data
  #nxmx['data'] = nxmx['instrument/detector/data']
