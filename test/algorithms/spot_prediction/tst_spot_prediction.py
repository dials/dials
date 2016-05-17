from __future__ import division

class TestSpotPredictor:

  def __init__(self):
    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import ScanStaticRayPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from iotbx.xds import xparm, integrate_hkl
    from dials.util import ioutil
    from math import ceil
    from os.path import realpath, dirname, join
    import dxtbx
    from rstbx.cftbx.coordinate_frame_converter import \
        coordinate_frame_converter
    from scitbx import matrix

    # The XDS files to read from
    test_path = dirname(dirname(dirname(realpath(__file__))))
    integrate_filename = join(test_path, 'data/sim_mx/INTEGRATE.HKL')
    gxparm_filename = join(test_path, 'data/sim_mx/GXPARM.XDS')

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

    assert(len(self.detector) == 1)

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
                                int(ceil(max(zcal)))))

    # Create the index generator
    generate_indices = IndexGenerator(self.unit_cell, self.space_group_type,
                                      self.d_min)

    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    fixed_rotation = self.gonio.get_fixed_rotation()
    UB = self.ub_matrix
    dphi = self.scan.get_oscillation_range(deg=False)

    # Create the ray predictor
    self.predict_rays = ScanStaticRayPredictor(s0, m2, fixed_rotation, dphi)

    # Predict the spot locations
    self.reflections = self.predict_rays(
                                    generate_indices.to_array(), UB)

    # Calculate the intersection of the detector and reflection frames
    success = ray_intersection(self.detector, self.reflections)
    self.reflections.select(success)

  def test_dmin(self):
    """Ensure calculated d_min < d_min in integrate file"""
    d = [self.unit_cell.d(h) for h in self.integrate_handle.hkl]
    d_min = min(d)
    assert(self.d_min <= d_min)
    print "OK"

  def test_miller_index_set(self):
    """Ensure we have the whole set of miller indices"""
    gen_hkl = {}
    for r in self.reflections:
      gen_hkl[r['miller_index']] = True
    for hkl in self.integrate_handle.hkl:
      assert(gen_hkl[hkl])

    print "OK"

  def test_rotation_angles(self):
    """Ensure the rotation angles agree with XDS"""

    # Create a dict of lists of xy for each hkl
    gen_phi = {}
    for r in self.reflections:
      hkl = r['miller_index']
      phi = r['phi']
      try:
        a = gen_phi[hkl]
        a.append(phi)
        gen_phi[hkl] = a
      except KeyError:
        gen_phi[hkl] = [phi]

    for hkl, xyz in zip(self.integrate_handle.hkl,
                        self.integrate_handle.xyzcal):

      xds_phi = self.scan.get_oscillation(deg=False)[0] + \
                xyz[2]*self.scan.get_oscillation(deg=False)[1]

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

      assert(abs(xds_phi - my_phi) < 0.1)

    print "OK"

  def test_beam_vectors(self):
    """Ensure |s1| == |s0|"""
    from scitbx import matrix
    s0_length = matrix.col(self.beam.get_s0()).length()
    for r in self.reflections:
      s1 = r['s1']
      s1_length = matrix.col(s1).length()
      assert(abs(s0_length - s1_length) < 1e-7)

    print "OK"

  def test_image_coordinates(self):
    """Ensure the image coordinates agree with XDS"""
    from scitbx import matrix

    # Create a dict of lists of xy for each hkl
    gen_xy = {}
    for r in self.reflections:
      hkl = r['miller_index']
      xy  = r['xyzcal.mm'][0:2]
      xy = self.detector[0].millimeter_to_pixel(xy)
      try:
        a = gen_xy[hkl]
        a.append(xy)
        gen_xy[hkl] = a
      except KeyError:
        gen_xy[hkl] = [xy]

    for hkl, xyz in zip(self.integrate_handle.hkl,
                        self.integrate_handle.xyzcal):

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

      if (abs(xds_xy[0] - my_xy[0]) > 0.1 or
          abs(xds_xy[1] - my_xy[1]) > 0.1):
        print xds_xy, gen_xy[hkl]
      assert(abs(xds_xy[0] - my_xy[0]) < 0.1)
      assert(abs(xds_xy[1] - my_xy[1]) < 0.1)

    print "OK"

  def run(self):
    self.test_dmin()
    self.test_miller_index_set()
    self.test_rotation_angles()
    self.test_beam_vectors()
    self.test_image_coordinates()
    print "OK"

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = TestSpotPredictor()
    test.run()
