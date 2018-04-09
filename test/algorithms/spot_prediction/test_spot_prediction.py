from __future__ import absolute_import, division, print_function

import os
import math
import pytest

class SpotPredictor:
  def __init__(self, dials_regression):
    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import ScanStaticRayPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from iotbx.xds import xparm, integrate_hkl
    from dials.util import ioutil
    import dxtbx
    from rstbx.cftbx.coordinate_frame_converter import \
        coordinate_frame_converter
    from scitbx import matrix

    # The XDS files to read from
    integrate_filename = os.path.join(dials_regression, 'data', 'sim_mx', 'INTEGRATE.HKL')
    gxparm_filename = os.path.join(dials_regression, 'data', 'sim_mx', 'GXPARM.XDS')

    # Read the XDS files
    self.integrate_handle = integrate_hkl.reader()
    self.integrate_handle.read_file(integrate_filename)
    self.gxparm_handle = xparm.reader()
    self.gxparm_handle.read_file(gxparm_filename)

    # Get the parameters we need from the GXPARM file
    models = dxtbx.load(gxparm_filename)
    self.beam = models.get_beam()
    self.gonio = models.get_goniometer()
    self.detector = models.get_detector()
    self.scan = models.get_scan()

    assert len(self.detector) == 1

    #print self.detector

    # Get crystal parameters
    self.space_group_type = ioutil.get_space_group_type_from_xparm(
        self.gxparm_handle)
    cfc = coordinate_frame_converter(gxparm_filename)
    a_vec = cfc.get('real_space_a')
    b_vec = cfc.get('real_space_b')
    c_vec = cfc.get('real_space_c')
    self.unit_cell = cfc.get_unit_cell()
    self.ub_matrix = matrix.sqr(a_vec + b_vec + c_vec).inverse()

    # Get the minimum resolution in the integrate file
    self.d_min = self.detector[0].get_max_resolution_at_corners(
        self.beam.get_s0())

    # Get the number of frames from the max z value
    xcal, ycal, zcal = zip(*self.integrate_handle.xyzcal)
    self.scan.set_image_range((self.scan.get_image_range()[0],
                               self.scan.get_image_range()[0] +
                                int(math.ceil(max(zcal)))))

    # Create the index generator
    generate_indices = IndexGenerator(self.unit_cell, self.space_group_type,
                                      self.d_min)

    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    fixed_rotation = self.gonio.get_fixed_rotation()
    setting_rotation = self.gonio.get_setting_rotation()
    UB = self.ub_matrix
    dphi = self.scan.get_oscillation_range(deg=False)

    # Create the ray predictor
    self.predict_rays = ScanStaticRayPredictor(s0, m2, fixed_rotation,
                                               setting_rotation, dphi)

    # Predict the spot locations
    self.reflections = self.predict_rays(
                                    generate_indices.to_array(), UB)

    # Calculate the intersection of the detector and reflection frames
    success = ray_intersection(self.detector, self.reflections)
    self.reflections.select(success)

@pytest.fixture(scope="session")
def spotpredictor(dials_regression):
  return SpotPredictor(dials_regression)

def test_dmin(spotpredictor):
  """Ensure calculated d_min < d_min in integrate file"""
  d = [spotpredictor.unit_cell.d(h) for h in spotpredictor.integrate_handle.hkl]
  d_min = min(d)
  assert spotpredictor.d_min <= d_min

def test_miller_index_set(spotpredictor):
  """Ensure we have the whole set of miller indices"""
  gen_hkl = {}
  for r in spotpredictor.reflections:
    gen_hkl[r['miller_index']] = True
  for hkl in spotpredictor.integrate_handle.hkl:
    assert(gen_hkl[hkl])

def test_rotation_angles(spotpredictor):
  """Ensure the rotation angles agree with XDS"""

  # Create a dict of lists of xy for each hkl
  gen_phi = {}
  for r in spotpredictor.reflections:
    hkl = r['miller_index']
    phi = r['phi']
    try:
      a = gen_phi[hkl]
      a.append(phi)
      gen_phi[hkl] = a
    except KeyError:
      gen_phi[hkl] = [phi]

  for hkl, xyz in zip(spotpredictor.integrate_handle.hkl,
                      spotpredictor.integrate_handle.xyzcal):

    xds_phi = spotpredictor.scan.get_oscillation(deg=False)[0] + \
              xyz[2]*spotpredictor.scan.get_oscillation(deg=False)[1]

    # Select the nearest xy to use if there are 2
    my_phi = gen_phi[hkl]
    if len(my_phi) == 2:
      my_phi0 = my_phi[0]
      my_phi1 = my_phi[1]
      diff0 = abs(xds_phi - my_phi0)
      diff1 = abs(xds_phi - my_phi1)
      if diff0 < diff1:
        my_phi = my_phi0
      else:
        my_phi = my_phi1
    else:
      my_phi = my_phi[0]

    assert xds_phi == pytest.approx(my_phi, abs=0.1)

def test_beam_vectors(spotpredictor):
  """Ensure |s1| == |s0|"""
  from scitbx import matrix
  s0_length = matrix.col(spotpredictor.beam.get_s0()).length()
  for r in spotpredictor.reflections:
    s1 = r['s1']
    s1_length = matrix.col(s1).length()
    assert s0_length == pytest.approx(s1_length, abs=1e-7)

def test_image_coordinates(spotpredictor):
  """Ensure the image coordinates agree with XDS"""
  from scitbx import matrix

  # Create a dict of lists of xy for each hkl
  gen_xy = {}
  for r in spotpredictor.reflections:
    hkl = r['miller_index']
    xy  = r['xyzcal.mm'][0:2]
    xy = spotpredictor.detector[0].millimeter_to_pixel(xy)
    try:
      a = gen_xy[hkl]
      a.append(xy)
      gen_xy[hkl] = a
    except KeyError:
      gen_xy[hkl] = [xy]

  for hkl, xyz in zip(spotpredictor.integrate_handle.hkl,
                      spotpredictor.integrate_handle.xyzcal):

    xds_xy = (xyz[0] - 0.5, xyz[1] - 0.5)

    # Select the nearest xy to use if there are 2
    my_xy = gen_xy[hkl]
    if len(my_xy) == 2:
      my_xy0 = my_xy[0]
      my_xy1 = my_xy[1]
      diff0 = (matrix.col(xds_xy) - matrix.col(my_xy0)).length()
      diff1 = (matrix.col(xds_xy) - matrix.col(my_xy1)).length()
      if diff0 < diff1:
        my_xy = my_xy0
      else:
        my_xy = my_xy1
    else:
      my_xy = my_xy[0]

    assert xds_xy[0] == pytest.approx(my_xy[0], abs=0.1), (xds_xy, gen_xy[hkl])
    assert xds_xy[1] == pytest.approx(my_xy[1], abs=0.1), (xds_xy, gen_xy[hkl])
