from __future__ import absolute_import, division, print_function

import math
import os

def test(dials_regression):
  from dxtbx.serialize import load
  from dials.algorithms.profile_model.gaussian_rs import Model
  from dials.algorithms.profile_model.gaussian_rs import MaskCalculator3D
  from dxtbx.model.experiment_list import Experiment, ExperimentList

  # Set the sweep filename and load the sweep
  sweep_filename = os.path.join(dials_regression, 'centroid_test_data', 'sweep.json')
  crystal_filename = os.path.join(dials_regression, 'centroid_test_data', 'crystal.json')

  # Load the sweep
  sweep = load.imageset(sweep_filename)
  crystal = load.crystal(crystal_filename)
  beam = sweep.get_beam()
  detector = sweep.get_detector()
  goniometer = sweep.get_goniometer()
  scan = sweep.get_scan()
  delta_d = 3 * beam.get_sigma_divergence(deg=False)
  try:
    mosaicity = crystal.get_mosaicity(deg=False)
  except AttributeError:
    mosaicity = 0
  delta_m = 3 * mosaicity
  nsigma = 3
  profile_model = Model(
    None,
    nsigma,
    beam.get_sigma_divergence(deg=False),
    mosaicity)
  experiment = ExperimentList()
  experiment.append(Experiment(
    imageset = sweep,
    beam = beam,
    detector = detector,
    goniometer = goniometer,
    scan = scan,
    crystal = crystal,
    profile = profile_model))

  assert len(detector) == 1

  # Get the function object to mask the foreground
  mask_foreground = MaskCalculator3D(beam, detector,
      goniometer, scan, delta_d, delta_m)

  from scitbx import matrix
  from scitbx.array_family import flex
  from dials.algorithms.shoebox import MaskCode
  from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
  s0 = beam.get_s0()
  m2 = goniometer.get_rotation_axis()
  s0_length = matrix.col(beam.get_s0()).length()
  width, height = detector[0].get_image_size()
  zrange = scan.get_array_range()
  phi0, dphi = scan.get_oscillation(deg=False)

  # Generate some reflections
  reflections = generate_reflections(detector, beam, scan, experiment, 10)

  # Mask the foreground in each
  mask_foreground(
    reflections['shoebox'],
    reflections['s1'],
    reflections['xyzcal.px'].parts()[2],
    reflections['panel'])

  # Loop through all the reflections and check the mask values
  shoebox = reflections['shoebox']
  beam_vector = reflections['s1']
  rotation_angle = reflections['xyzcal.mm'].parts()[2]
  for l in range(len(reflections)):
    mask = shoebox[l].mask
    x0, x1, y0, y1, z0, z1 = shoebox[l].bbox
    s1 = beam_vector[l]
    phi = rotation_angle[l]
    cs = CoordinateSystem(m2, s0, s1, phi)

    def rs_coord(i, j, k):
      s1d = detector[0].get_pixel_lab_coord((i,j))
      s1d = matrix.col(s1d).normalize() * s0_length
      e1, e2 = cs.from_beam_vector(s1d)
      e3 = cs.from_rotation_angle_fast(phi0 + (k - zrange[0]) * dphi)
      return e1, e2, e3

    new_mask = flex.int(mask.accessor(), 0)
    for k in range(z1 - z0):
      for j in range(y1 - y0):
        for i in range(x1 - x0):
          #value1 = mask[k, j, i]
          e11, e12, e13 = rs_coord(x0 + i,     y0 + j,     z0 + k)
          e21, e22, e23 = rs_coord(x0 + i + 1, y0 + j,     z0 + k)
          e31, e32, e33 = rs_coord(x0 + i,     y0 + j + 1, z0 + k)
          e41, e42, e43 = rs_coord(x0 + i,     y0 + j,     z0 + k + 1)
          e51, e52, e53 = rs_coord(x0 + i + 1, y0 + j + 1, z0 + k)
          e61, e62, e63 = rs_coord(x0 + i + 1, y0 + j,     z0 + k + 1)
          e71, e72, e73 = rs_coord(x0 + i,     y0 + j + 1, z0 + k + 1)
          e81, e82, e83 = rs_coord(x0 + i + 1, y0 + j + 1, z0 + k + 1)
          de1 = (e11/delta_d)**2+(e12/delta_d)**2#+(e13/delta_m)**2
          de2 = (e21/delta_d)**2+(e22/delta_d)**2#+(e23/delta_m)**2
          de3 = (e31/delta_d)**2+(e32/delta_d)**2#+(e33/delta_m)**2
          de4 = (e41/delta_d)**2+(e42/delta_d)**2#+(e43/delta_m)**2
          de5 = (e51/delta_d)**2+(e52/delta_d)**2#+(e53/delta_m)**2
          de6 = (e61/delta_d)**2+(e62/delta_d)**2#+(e63/delta_m)**2
          de7 = (e71/delta_d)**2+(e72/delta_d)**2#+(e73/delta_m)**2
          de8 = (e81/delta_d)**2+(e82/delta_d)**2#+(e83/delta_m)**2
          de = math.sqrt(min([de1, de2, de3, de4, de5, de6, de7, de8]))
          gx = min([e11, e21, e31, e41, e51, e61, e71, e81])
          gy = min([e12, e22, e32, e42, e52, e62, e72, e82])
          gz = min([e13, e23, e33, e43, e53, e63, e73, e83])
          dee = (gx/delta_d)**2 + (gy/delta_d)**2# + (gz/delta_m)**2
          #print math.sqrt(dee), de
          if (x0 + i < 0 or y0 + j < 0 or
              x0 + i >= width or y0 + j >= height or
              z0 + k < zrange[0] or z0 + k >= zrange[1]):
            value2 = MaskCode.Valid
          else:
            if de <= 1.0:
              value2 = MaskCode.Valid | MaskCode.Foreground
            else:
              value2 = MaskCode.Valid | MaskCode.Background
          new_mask[k,j,i] = value2

    try:
      assert(all(m1 == m2 for m1, m2 in zip(mask, new_mask)))
    except Exception:
      import numpy
      numpy.set_printoptions(threshold=10000)
      diff = (mask == new_mask).as_numpy_array()
      print(diff.astype(numpy.int))
      #print mask.as_numpy_array()
      #print new_mask.as_numpy_array()
      #print (new_mask.as_numpy_array()[:,:,:] %2) * (new_mask.as_numpy_array() == 5)
      raise

def generate_reflections(detector, beam, scan, experiment, num):
  from random import randint, seed
  from scitbx import matrix
  from dials.array_family import flex
  from dials.algorithms.shoebox import MaskCode
  seed(0)
  assert(len(detector) == 1)
  beam_vector = flex.vec3_double(num)
  xyzcal_px = flex.vec3_double(num)
  xyzcal_mm = flex.vec3_double(num)
  panel = flex.size_t(num)
  s0_length = matrix.col(beam.get_s0()).length()
  for i in range(num):
    x = randint(0, 2000)
    y = randint(0, 2000)
    z = randint(0, 8)
    s1 = detector[0].get_pixel_lab_coord((x, y))
    s1 = matrix.col(s1).normalize() * s0_length
    phi = scan.get_angle_from_array_index(z, deg=False)
    beam_vector[i] = s1
    xyzcal_px[i] = (x, y, z)
    (x, y)  = detector[0].pixel_to_millimeter((x, y))
    xyzcal_mm[i] = (x, y, phi)
    panel[i] = 0

  sigma_b = experiment[0].beam.get_sigma_divergence(deg=False)
  try:
    sigma_m = experiment[0].crystal.get_mosaicity(deg=False)
  except AttributeError:
    sigma_m = 0

  rlist = flex.reflection_table()
  rlist['id'] = flex.int(len(beam_vector), 0)
  rlist['s1'] = beam_vector
  rlist['panel'] = panel
  rlist['xyzcal.px'] = xyzcal_px
  rlist['xyzcal.mm'] = xyzcal_mm
  rlist['bbox'] = rlist.compute_bbox(experiment)
  index = []
  image_size = experiment[0].detector[0].get_image_size()
  array_range = experiment[0].scan.get_array_range()
  bbox = rlist['bbox']
  for i in range(len(rlist)):
    x0, x1, y0, y1, z0, z1 = bbox[i]
    if (x0 < 0 or x1 > image_size[0] or
        y0 < 0 or y1 > image_size[1] or
        z0 < array_range[0] or z1 > array_range[1]):
      index.append(i)
  rlist.del_selected(flex.size_t(index))
  rlist['shoebox'] = flex.shoebox(
    rlist['panel'], rlist['bbox'])
  rlist['shoebox'].allocate_with_value(MaskCode.Valid)
  return rlist
