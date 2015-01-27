from __future__ import division

def export_beam(outfile, beam):
  ''' Export the beam model. '''
  from scitbx import matrix

  EPS = 1e-7

  # Make sure that the direction is 0, 0, -1
  assert(matrix.col(beam.get_direction()).dot(matrix.col((0, 0, -1))) < EPS)

  # Get the nx_beam
  nx_beam = outfile.entry.sample.beam
  nx_beam['incident_wavelength'] = beam.get_wavelength()


def export_detector(outfile, detector):
  ''' Export the detector model. '''
  from scitbx import matrix

  # Get the panel
  panel = detector[0]

  # Get some panel attributes
  pixel_size = panel.get_pixel_size()
  image_size = panel.get_image_size()
  origin = matrix.col(panel.get_origin())

  # Get the detector module object
  detector = outfile.entry.instrument.detector
  module = outfile.entry.instrument.detector.module['module0']

  # Set the data size
  module['data_size'] = image_size

  # Create the detector translation
  transformations = detector.transformations
  transformations['translation'] = 0
  transformations['translation'].attrs['depends_on'] = '.'
  transformations['translation'].attrs['transformation_type'] = 'translation'
  transformations['translation'].attrs['units'] = 'mm'
  transformations['translation'].attrs['vector'] = (0, 0, 1)

  # Get the path for below
  translation_path = str('%s/%s' % (transformations.path(), 'translation'))

  # Create the detector depends on
  detector['depends_on'] = translation_path

  # Set the module offset
  module['module_offset'] = origin.length()
  module['module_offset'].attrs['depends_on'] = translation_path
  module['module_offset'].attrs['transformation_type'] = 'translation'
  module['module_offset'].attrs['units'] = 'mm'
  module['module_offset'].attrs['vector'] = origin.normalize()

  # The path for items below
  module_offset_path = str('%s/%s' % (module.path(), 'module_offset'))

  # Write the fast pixel direction
  module['fast_pixel_direction'] = pixel_size[0]
  module['fast_pixel_direction'].attrs['depends_on'] = module_offset_path
  module['fast_pixel_direction'].attrs['transformation_type'] = 'translation'
  module['fast_pixel_direction'].attrs['units'] = 'mm'
  module['fast_pixel_direction'].attrs['vector'] = panel.get_fast_axis()

  # Write the slow pixel direction
  module['slow_pixel_direction'] = pixel_size[1]
  module['slow_pixel_direction'].attrs['depends_on'] = module_offset_path
  module['slow_pixel_direction'].attrs['transformation_type'] = 'translation'
  module['slow_pixel_direction'].attrs['units'] = 'mm'
  module['slow_pixel_direction'].attrs['vector'] = panel.get_slow_axis()

  # Write the fast pixel size
  module['fast_pixel_size'] = pixel_size[0]
  module['fast_pixel_size'].attrs['depends_on'] = module_offset_path
  module['fast_pixel_size'].attrs['transformation_type'] = 'translation'
  module['fast_pixel_size'].attrs['units'] = 'mm'
  module['fast_pixel_size'].attrs['vector'] = panel.get_fast_axis()

  # Write the slow pixel size
  module['slow_pixel_size'] = pixel_size[1]
  module['slow_pixel_size'].attrs['depends_on'] = module_offset_path
  module['slow_pixel_size'].attrs['transformation_type'] = 'translation'
  module['slow_pixel_size'].attrs['units'] = 'mm'
  module['slow_pixel_size'].attrs['vector'] = panel.get_slow_axis()

def export_goniometer(outfile, goniometer, scan):
  ''' Export the goniometer model. '''

  # The angles for each image
  phi0, dphi = scan.get_oscillation(deg=True)
  phi = [phi0+dphi*i for i in range(len(scan))]

  # Write out the rotation axis and oscillation
  transformations = outfile.entry.sample.transformations
  transformations['phi'] = phi
  transformations['phi'].attrs['depends_on'] = '.'
  transformations['phi'].attrs['transformation_type'] = 'rotation'
  transformations['phi'].attrs['units'] = 'deg'
  transformations['phi'].attrs['vector'] = goniometer.get_rotation_axis()

def export_crystal(outfile, crystal):
  ''' Export the crystal model. '''

  from scitbx.array_family import flex

  # Get the sample
  sample = outfile.entry.sample

  # Set the space group
  sample['unit_cell_group'] = crystal.get_space_group().type().hall_symbol()

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
    unit_cell = crystal.get_unit_cell().parameters()
    orientation_matrix = crystal.get_U().transpose()

  # Save the unit cell data
  sample['unit_cell'] = unit_cell
  sample['unit_cell'].attrs['angles_units'] = 'deg'
  sample['unit_cell'].attrs['length_units'] = 'angstrom'

  # Save the orientation matrix
  sample['orientation_matrix'] = orientation_matrix

  # Set depends on
  sample['depends_on'] = str('%s/%s' % (sample.transformations.path(), 'phi'))


def export_experiments(outfile, experiments):
  ''' Export the experiments to the NXmx file. '''

  # Ensure only 1 experiment at the moment
  assert(len(experiments) == 1)
  experiment = experiments[0]

  # Ensure only 1 panel at the moment
  assert(len(experiment.detector) == 1)

  # Do crystal change of basis
  space_group_info = experiment.crystal.get_space_group().info()
  cb_op_to_ref = space_group_info.change_of_basis_op_to_reference_setting()
  crystal = experiment.crystal.change_basis(cb_op_to_ref)

  # Export the beam model
  export_beam(outfile, experiment.beam)

  # Export the detector
  export_detector(outfile, experiment.detector)

  # Export the goniometer
  export_goniometer(outfile, experiment.goniometer, experiment.scan)

  # Export the crystal
  export_crystal(outfile, crystal)

def export_reflections(outfile, reflections):
  ''' Export the reflection table to the NXmx file. '''

  # Ensure only a single experiment
  assert(reflections['id'].all_eq(0))

  # Ensure only a single panel
  assert(reflections['panel'].all_eq(0))

  # Select only those reflections which have valid summed intensities and where
  # present, profile fitted intensities. Also only export fully recorded
  # reflections
  selection = reflections['intensity.sum.variance'] <= 0
  if selection.count(True) > 0:
    reflections.del_selected(selection)
    print 'Removing %d reflections with negative variance...' % \
          selection.count(True)

  if 'intensity.prf.variance' in reflections:
    selection = reflections['intensity.prf.variance'] <= 0
    if selection.count(True) > 0:
      reflections.del_selected(selection)
      print 'Removing %d profile reflections with negative variance...' % \
            selection.count(True)

  if 'partiality' in reflections:
    selection = reflections['partiality'] < 0.99
    if selection.count(True) > 0:
      reflections.del_selected(selection)
      print 'Removing %d incomplete reflections...' % \
        selection.count(True)

  # Get the diffraction class
  diffraction = outfile.entry.diffraction

  # Export all the columns
  try:
    col1, col2, col3 = zip(*list(reflections['miller_index']))
    diffraction['h'] = col1
    diffraction['k'] = col2
    diffraction['l'] = col3
    diffraction['h'].attrs['description'] = 'The h component of the miller index'
    diffraction['k'].attrs['description'] = 'The k component of the miller index'
    diffraction['l'].attrs['description'] = 'The l component of the miller index'
  except Exception:
    pass

  try:
    diffraction['id'] = reflections['id']
    diffraction['id'].attrs['description'] = 'The experiment id'
  except Exception:
    pass

  try:
    diffraction['int_sum_val'] = reflections['intensity.sum.value']
    diffraction['int_sum_var'] = reflections['intensity.sum.variance']
    diffraction['int_sum_val'].attrs['description'] = 'The value of the summed intensity'
    diffraction['int_sum_var'].attrs['description'] = 'The variance of the summed intensity'
  except Exception:
    pass

  try:
    diffraction['int_prf_val'] = reflections['intensity.prf.value']
    diffraction['int_prf_var'] = reflections['intensity.prf.variance']
    diffraction['int_prf_val'].attrs['description'] = 'The value of the profile fitted intensity'
    diffraction['int_prf_var'].attrs['description'] = 'The variance of the profile fitted intensity'
  except Exception:
    pass

  try:
    diffraction['lp'] = reflections['lp']
    diffraction['lp'].attrs['description'] = 'The lorentz-polarization correction factor'
  except Exception:
    pass

  try:
    diffraction['det_module'] = reflections['panel']
    diffraction['det_module'].attrs['description'] = 'The detector module on which the reflection was recorded'
  except Exception:
    pass

  try:
    col11, col12, col13, col14, col15, col16 = reflections['bbox'].parts()
    diffraction['bbx0'] = col11
    diffraction['bbx1'] = col12
    diffraction['bby0'] = col13
    diffraction['bby1'] = col14
    diffraction['bbz0'] = col15
    diffraction['bbz1'] = col16
    diffraction['bbx0'].attrs['description'] = 'The bounding box lower x bound'
    diffraction['bbx1'].attrs['description'] = 'The bounding box upper x bound'
    diffraction['bby0'].attrs['description'] = 'The bounding box lower y bound'
    diffraction['bby1'].attrs['description'] = 'The bounding box upper y bound'
    diffraction['bbz0'].attrs['description'] = 'The bounding box lower z bound'
    diffraction['bbz1'].attrs['description'] = 'The bounding box upper z bound'
  except Exception:
    pass

  try:
    col17, col18, col19 = reflections['xyzcal.px'].parts()
    diffraction['prd_px_x'] = col17
    diffraction['prd_px_y'] = col18
    diffraction['prd_frame'] = col19
    diffraction['prd_px_x'].attrs['description'] = 'The predicted bragg peak fast pixel location'
    diffraction['prd_px_y'].attrs['description'] = 'The predicted bragg peak slow pixel location'
    diffraction['prd_frame'].attrs['description'] = 'The predicted bragg peak frame number'
  except Exception:
    pass

  try:
    col20, col21, col22 = reflections['xyzcal.mm'].parts()
    diffraction['prd_mm_x'] = col20
    diffraction['prd_mm_y'] = col21
    diffraction['prd_phi'] = col22
    diffraction['prd_mm_x'].attrs['description'] = 'The predicted bragg peak fast millimeter location'
    diffraction['prd_mm_y'].attrs['description'] = 'The predicted bragg peak slow millimeter location'
    diffraction['prd_phi'].attrs['description'] = 'The predicted bragg peak rotation angle number'
  except Exception:
    pass

  try:
    col23, col24, col25 = reflections['xyzobs.px.value'].parts()
    col26, col27, col28 = reflections['xyzobs.px.variance'].parts()
    diffraction['obs_px_x_val'] = col23
    diffraction['obs_px_x_var'] = col26
    diffraction['obs_px_y_val'] = col24
    diffraction['obs_px_y_var'] = col27
    diffraction['obs_frame_val'] = col25
    diffraction['obs_frame_var'] = col28
    diffraction['obs_px_x_val'].attrs['description'] = 'The observed centroid fast pixel value'
    diffraction['obs_px_x_var'].attrs['description'] = 'The observed centroid fast pixel variance'
    diffraction['obs_px_y_val'].attrs['description'] = 'The observed centroid slow pixel value'
    diffraction['obs_px_y_var'].attrs['description'] = 'The observed centroid slow pixel variance'
    diffraction['obs_frame_val'].attrs['description'] = 'The observed centroid frame value'
    diffraction['obs_frame_var'].attrs['description'] = 'The observed centroid frame variance'
  except Exception:
    pass

  try:
    col29, col30, col31 = reflections['xyzobs.mm.value'].parts()
    col32, col33, col34 = reflections['xyzobs.mm.variance'].parts()
    diffraction['obs_mm_x_val'] = col29
    diffraction['obs_mm_x_var'] = col32
    diffraction['obs_mm_y_val'] = col30
    diffraction['obs_mm_y_var'] = col33
    diffraction['obs_phi_val'] = col31
    diffraction['obs_phi_var'] = col34
    diffraction['obs_mm_x_val'].attrs['description'] = 'The observed centroid fast millimeter value'
    diffraction['obs_mm_x_var'].attrs['description'] = 'The observed centroid fast millimeter variance'
    diffraction['obs_mm_y_val'].attrs['description'] = 'The observed centroid slow millimeter value'
    diffraction['obs_mm_y_var'].attrs['description'] = 'The observed centroid slow millimeter variance'
    diffraction['obs_phi_val'].attrs['description'] = 'The observed centroid phi value'
    diffraction['obs_phi_var'].attrs['description'] = 'The observed centroid phi variance'
  except Exception:
    pass

  try:
    diffraction['partiality'] = reflections['partiality']
    diffraction['partiality'].attrs['description'] = 'The partiality of the reflection'
  except Exception:
    pass

  try:
    diffraction['d'] = reflections['d']
    diffraction['d'].attrs['description'] = 'The resolution of the reflection'
  except Exception:
    pass

  try:
    diffraction['bkg_mean'] = reflections['background.mean']
    diffraction['bkg_mean'].attrs['description'] = 'The mean background value'
  except Exception:
    pass

  try:
    diffraction['entering'] = reflections['entering']
    diffraction['entering'].attrs['description'] = 'Entering or exiting the Ewald sphere'
  except Exception:
    pass

  try:
    diffraction['flags'] = reflections['flags']
    diffraction['flags'].attrs['description'] = 'Status of the reflection in processing'
  except Exception:
    pass

  try:
    diffraction['prf_cc'] = reflections['profile.correlation']
    diffraction['prf_cc'].attrs['description'] = 'Profile fitting correlations'
  except Exception:
    pass

def export_details(outfile):
  from time import strftime

  # Get the entry
  entry = outfile.entry

  # Program info
  entry['program_name'] = 'dials.export_nxmx'
  entry['program_name'].attrs['version'] = 1
  entry['program_name'].attrs['configuration'] = ''

  # Set some processing information (each program should add itself)
  process = entry.process['process']
  process['program'] = 'dials'
  process['version'] = 1
  process['date'] = strftime('%Y-%m-%dT%H:%M:%S')

  note = process.note['0_spot_finding']
  note['author'] = 'dials.find_spots'
  note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  note['type'] = 'text/plain'
  note['description'] = 'Spot finding parameters'
  note['data'] = 'dials.find_spots datablock.json'

  note = process.note['1_indexing']
  note['author'] = 'dials.index'
  note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  note['type'] = 'text/plain'
  note['description'] = 'Indexing parameters'
  note['data'] = 'dials.index datablock.json strong.pickle'

  note = process.note['2_refinement']
  note['author'] = 'dials.refine'
  note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  note['type'] = 'text/plain'
  note['description'] = 'Refinement parameters'
  note['data'] = 'dials.refine experiments.json indexed.pickle'

  note = process.note['3_integration']
  note['author'] = 'dials.integrate'
  note['date'] = strftime('%Y-%m-%dT%H:%M:%S')
  note['type'] = 'text/plain'
  note['description'] = 'Integration parameters'
  note['data'] = 'dials.integrate refined_experiments.json indexed.pickle'

def import_beam(infile):
  from dxtbx.model import Beam

  # Get the nx_beam
  nx_beam = infile.entry.sample.beam
  wavelength = nx_beam['incident_wavelength'].value

  # Return the beam model
  return Beam((0, 0, -1), wavelength)

def import_detector(infile):
  from dxtbx.model import Detector
  from scitbx import matrix

  # Get the detector module object
  detector = infile.entry.instrument.detector
  module = detector.module['module0']

  # Set the data size
  image_size = module['data_size']

  # Create the detector translation
  transformations = detector.transformations
  assert(transformations['translation'].value == 0)
  assert(transformations['translation'].attrs['depends_on'] == '.')
  assert(transformations['translation'].attrs['transformation_type'] == 'translation')
  assert(transformations['translation'].attrs['units'] == 'mm')
  assert(tuple(transformations['translation'].attrs['vector']) == (0, 0, 1))

  # Get the path for below
  translation_path = str('%s/%s' % (transformations.path(), 'translation'))

  # Create the detector depends on
  assert(detector['depends_on'].value == translation_path)

  # Set the module offset
  offset_length = module['module_offset'].value
  assert(module['module_offset'].attrs['depends_on'] == translation_path)
  assert(module['module_offset'].attrs['transformation_type'] == 'translation')
  assert(module['module_offset'].attrs['units'] == 'mm')
  offset_vector = matrix.col(module['module_offset'].attrs['vector'])
  origin = offset_vector * offset_length

  # The path for items below
  module_offset_path = str('%s/%s' % (module.path(), 'module_offset'))

  # Write the fast pixel direction
  pixel_size_x1 = module['fast_pixel_direction'].value
  assert(module['fast_pixel_direction'].attrs['depends_on'] == module_offset_path)
  assert(module['fast_pixel_direction'].attrs['transformation_type'] == 'translation')
  assert(module['fast_pixel_direction'].attrs['units'] == 'mm')
  fast_axis1 = tuple(module['fast_pixel_direction'].attrs['vector'])

  # Write the slow pixel direction
  pixel_size_y1 = module['slow_pixel_direction'].value
  assert(module['slow_pixel_direction'].attrs['depends_on'] == module_offset_path)
  assert(module['slow_pixel_direction'].attrs['transformation_type'] == 'translation')
  assert(module['slow_pixel_direction'].attrs['units'] == 'mm')
  slow_axis1 = tuple(module['slow_pixel_direction'].attrs['vector'])

  # Write the fast pixel size
  pixel_size_x2 = module['fast_pixel_size'].value
  assert(module['fast_pixel_size'].attrs['depends_on'] == module_offset_path)
  assert(module['fast_pixel_size'].attrs['transformation_type'] == 'translation')
  assert(module['fast_pixel_size'].attrs['units'] == 'mm')
  fast_axis2 = tuple(module['fast_pixel_size'].attrs['vector'])

  # Write the slow pixel size
  pixel_size_y2 = module['slow_pixel_size'].value
  assert(module['slow_pixel_size'].attrs['depends_on'] == module_offset_path)
  assert(module['slow_pixel_size'].attrs['transformation_type'] == 'translation')
  assert(module['slow_pixel_size'].attrs['units'] == 'mm')
  slow_axis2 = tuple(module['slow_pixel_size'].attrs['vector'])

  # Get the pixel size and axis vectors
  assert(pixel_size_x1 == pixel_size_x2)
  assert(pixel_size_y1 == pixel_size_y2)
  assert(fast_axis1 == fast_axis2)
  assert(slow_axis1 == slow_axis2)
  pixel_size = (pixel_size_x1, pixel_size_y1)
  fast_axis = fast_axis1
  slow_axis = slow_axis1

  # FIXME stuff not set
  det_type = ''
  name = ''
  thickness = 0
  material = ''
  trusted_range = (0, 0)

  # Return the detector and panel
  detector = Detector()
  panel = detector.add_panel()
  panel.set_frame(fast_axis, slow_axis, origin)
  panel.set_pixel_size(pixel_size)
  panel.set_image_size(image_size)
  return detector

def import_goniometer(infile):
  from dxtbx.model import Goniometer

  # Write out the rotation axis and oscillation
  transformations = infile.entry.sample.transformations
  assert(transformations['phi'].attrs['depends_on'] == '.')
  assert(transformations['phi'].attrs['transformation_type'] == 'rotation')
  assert(transformations['phi'].attrs['units'] == 'deg')
  rotation_axis = tuple(transformations['phi'].attrs['vector'])

  # Return the goniometer model
  return Goniometer(rotation_axis)

def import_scan(infile):
  from dxtbx.model import Scan

  # Write out the rotation axis and oscillation
  transformations = infile.entry.sample.transformations
  phi = transformations['phi']
  assert(transformations['phi'].attrs['depends_on'] == '.')
  assert(transformations['phi'].attrs['transformation_type'] == 'rotation')
  assert(transformations['phi'].attrs['units'] == 'deg')
  image_range = (0, len(phi))
  oscillation = (phi[0], phi[1] - phi[0])
  return Scan(image_range, oscillation, deg=True)

def import_crystal(infile):
  from dxtbx.model.crystal import crystal_model as Crystal
  from scitbx.array_family import flex
  from scitbx import matrix

  # Get the sample
  sample = infile.entry.sample

  # Set the space group
  space_group_symbol = sample['unit_cell_group'].value

  # Save the unit cell data
  unit_cell = sample['unit_cell']
  assert(sample['unit_cell'].attrs['angles_units'] == 'deg')
  assert(sample['unit_cell'].attrs['length_units'] == 'angstrom')

  # Save the orientation matrix
  orientation_matrix = sample['orientation_matrix']

  if len(unit_cell.shape) == 1:
    assert(len(orientation_matrix.shape) == 1)
    unit_cell = unit_cell.value
    orientation_matrix = matrix.sqr(orientation_matrix.value)
  else:
    assert(len(unit_cell.shape) == 2)
    assert(len(orientation_matrix.shape) == 2)
    assert(unit_cell.shape[1] == 6)
    assert(orientation_matrix.shape[1] == 9)
    num = unit_cell.shape[0]
    assert(orientation_matrix.shape[0] == num)
    unit_cell = flex.double(unit_cell)
    orientation_matrix = flex.double(orientation_matrix)
    assert(False)

  U = orientation_matrix.inverse()
  real_space_a = U[0:3]
  real_space_b = U[3:6]
  real_space_c = U[6:9]

  # Set depends on
  assert(sample['depends_on'].value == str('%s/%s' % (sample.transformations.path(), 'phi')))

  # Construct the crystal model
  crystal = Crystal(
    real_space_a,
    real_space_b,
    real_space_c,
    space_group_symbol)

  # Return the crystal
  return crystal

def import_experiments(infile):
  from dxtbx.model.experiment.experiment_list import ExperimentList
  from dxtbx.model.experiment.experiment_list import Experiment

  # Create the experiment
  experiment = Experiment()

  # Import the beam model
  experiment.beam = import_beam(infile)

  # Import the detector
  experiment.detector = import_detector(infile)

  # Import the goniometer
  experiment.goniometer = import_goniometer(infile)

  # Import the scan
  experiment.scan = import_scan(infile)

  # Import the crystal
  experiment.crystal = import_crystal(infile)

  # Return the experiment list
  experiment_list = ExperimentList()
  experiment_list.append(experiment)
  return experiment_list

def import_reflections(outfile):
  ''' Import the reflection table to the NXmx file. '''
  from dials.array_family import flex

  # Create the reflection table
  reflections = flex.reflection_table()

  # Get the diffraction class
  diffraction = outfile.entry.diffraction

  # Export all the columns
  try:
    col1 = diffraction['h'].value
    col2 = diffraction['k'].value
    col3 = diffraction['l'].value
    reflections['miller_index'] = flex.miller_index(
      flex.int(col1),
      flex.int(col2),
      flex.int(col3))
  except Exception:
    pass

  try:
    reflections['id'] = flex.size_t(diffraction['id'].value)
  except Exception:
    pass

  try:
    reflections['intensity.sum.value'] = flex.double(diffraction['int_sum_val'].value)
    reflections['intensity.sum.variance'] = flex.double(diffraction['int_sum_var'].value)
  except Exception:
    pass

  try:
    reflections['intensity.prf.value'] = flex.double(diffraction['int_prf_val'].value)
    reflections['intensity.prf.variance'] = flex.double(diffraction['int_prf_var'].value)
  except Exception:
    pass

  try:
    reflections['lp'] = flex.double(diffraction['lp'].value)
  except Exception:
    pass

  try:
    reflections['panel'] = flex.size_t(diffraction['det_module'].value)
  except Exception:
    pass

  try:
    col11 = diffraction['bbx0'].value
    col12 = diffraction['bbx1'].value
    col13 = diffraction['bby0'].value
    col14 = diffraction['bby1'].value
    col15 = diffraction['bbz0'].value
    col16 = diffraction['bbz1'].value
    reflections['bbox'] = flex.int6(
      flex.int(col11),
      flex.int(col12),
      flex.int(col13),
      flex.int(col14),
      flex.int(col15),
      flex.int(col16))
  except Exception:
    pass

  try:
    col17 = diffraction['prd_px_x'].value
    col18 = diffraction['prd_px_y'].value
    col19 = diffraction['prd_frame'].value
    reflections['xyzcal.px'] = flex.vec3_double(
      flex.double(col17),
      flex.double(col18),
      flex.double(col19))
  except Exception:
    pass

  try:
    col20 = diffraction['prd_mm_x'].value
    col21 = diffraction['prd_mm_y'].value
    col22 = diffraction['prd_phi'].value
    reflections['xyzcal.mm'] = flex.vec3_double(
      flex.double(col20),
      flex.double(col21),
      flex.double(col22))
  except Exception:
    pass

  try:
    col23 = diffraction['obs_px_x_val'].value
    col26 = diffraction['obs_px_x_var'].value
    col24 = diffraction['obs_px_y_val'].value
    col27 = diffraction['obs_px_y_var'].value
    col25 = diffraction['obs_frame_val'].value
    col28 = diffraction['obs_frame_var'].value
    reflections['xyzobs.px.value'] = flex.vec3_double(
      flex.double(col23),
      flex.double(col24),
      flex.double(col25))
    reflections['xyzobs.px.variance'] = flex.vec3_double(
      flex.double(col26),
      flex.double(col27),
      flex.double(col28))
  except Exception:
    pass

  try:
    col29 = diffraction['obs_mm_x_val'].value
    col32 = diffraction['obs_mm_x_var'].value
    col30 = diffraction['obs_mm_y_val'].value
    col33 = diffraction['obs_mm_y_var'].value
    col31 = diffraction['obs_phi_val'].value
    col34 = diffraction['obs_phi_var'].value
    reflections['xyzobs.mm.value'] = flex.vec3_double(
      flex.double(col29),
      flex.double(col30),
      flex.double(col31))
    reflections['xyzobs.mm.variance'] = flex.vec3_double(
      flex.double(col32),
      flex.double(col33),
      flex.double(col34))
  except Exception:
    pass

  try:
    reflections['partiality'] = flex.double(diffraction['partiality'].value)
  except Exception:
    pass

  try:
    reflections['d'] = flex.double(diffraction['d'].value)
  except Exception:
    pass

  try:
    reflections['background.mean'] = flex.double(diffraction['bkg_mean'].value)
  except Exception:
    pass

  try:
    reflections['entering'] = flex.bool(map(bool,list(diffraction['entering'].value)))
  except Exception:
    pass

  try:
    reflections['flags'] = flex.int(diffraction['flags'].value)
  except Exception:
    pass

  try:
    reflections['profile.correlation'] = flex.double(diffraction['prf_cc'].value)
  except Exception:
    pass

  return reflections

def export(experiments, reflections, filename):
  ''' Export the experiments and reflections as an NXmx file. '''
  from dials.scratch.jmp.mtz2 import mtz2

  # Open the NXmx file for writing
  outfile = mtz2.File(filename, 'w')

  # Export the experiments
  print 'Exporting experimental models...'
  export_experiments(outfile, experiments)

  # Export the reflection table data
  print 'Exporting reflection data...'
  export_reflections(outfile, reflections)

  # Export some extra details
  print 'Exporting some additional processing details...'
  export_details(outfile)

  # Flush the file
  outfile.flush()
  print 'Wrote NXmx file %s' % filename

  # FIXME The following items have not been set from export_mtz
  # o.set_divhd(0.0).set_divvd(0.0)
  # o.set_bbfac(0.0).set_bscale(1.0)
  # o.set_sdbfac(0.0).set_sdbscale(0.0).set_nbscal(0)
  # o.set_lcrflg(0)
  # o.set_datum(flex.float((0.0, 0.0, 0.0)))
  # o.set_misflg(0)
  # o.set_jumpax(0)
  # o.set_ldtype(2)

def load(filename):
  ''' Import the experiments and reflections as an NXmx file. '''
  from dials.scratch.jmp.mtz2 import mtz2

  # Open the NXmx file for reading
  infile = mtz2.File(filename, 'r')

  # Import the experiments
  print 'Importing experimental models...'
  experiments = import_experiments(infile)

  # Import the reflection table data
  print 'Importing reflection data...'
  reflections = import_reflections(infile)

  print "Loaded %d reflections" % len(reflections)

  # Return experiment and reflections
  return experiments, reflections
