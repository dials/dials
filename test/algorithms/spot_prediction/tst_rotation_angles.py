from __future__ import absolute_import, division

import libtbx.load_env
from os.path import join, isdir
have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression", test=isdir)

def run():
  if not have_dials_regression:
    print "Skipping test: dials_regression not available."
    return

  from scitbx import matrix
  from iotbx.xds import xparm, integrate_hkl
  from dials.util import ioutil
  from math import ceil
  from dials.algorithms.spot_prediction import RotationAngles
  import dxtbx
  from rstbx.cftbx.coordinate_frame_converter import \
      coordinate_frame_converter

  # The XDS files to read from
  integrate_filename = join(dials_regression, 'data/sim_mx/INTEGRATE.HKL')
  gxparm_filename = join(dials_regression, 'data/sim_mx/GXPARM.XDS')

  # Read the XDS files
  integrate_handle = integrate_hkl.reader()
  integrate_handle.read_file(integrate_filename)
  gxparm_handle = xparm.reader()
  gxparm_handle.read_file(gxparm_filename)

  # Get the parameters we need from the GXPARM file
  models = dxtbx.load(gxparm_filename)
  beam = models.get_beam()
  gonio = models.get_goniometer()
  detector = models.get_detector()
  scan = models.get_scan()

  # Get the crystal parameters
  space_group_type = ioutil.get_space_group_type_from_xparm(gxparm_handle)
  cfc = coordinate_frame_converter(gxparm_filename)
  a_vec = cfc.get('real_space_a')
  b_vec = cfc.get('real_space_b')
  c_vec = cfc.get('real_space_c')
  unit_cell = cfc.get_unit_cell()
  UB = matrix.sqr(a_vec + b_vec + c_vec).inverse()
  ub_matrix = UB

  # Get the minimum resolution in the integrate file
  d = [unit_cell.d(h) for h in integrate_handle.hkl]
  d_min = min(d)

  # Get the number of frames from the max z value
  xcal, ycal, zcal = zip(*integrate_handle.xyzcal)
  num_frames = int(ceil(max(zcal)))
  scan.set_image_range((scan.get_image_range()[0],
                      scan.get_image_range()[0] + num_frames - 1))

  # Create the rotation angle object
  ra = RotationAngles(beam.get_s0(), gonio.get_rotation_axis())

  # Setup the matrices
  ub = matrix.sqr(ub_matrix)
  s0 = matrix.col(beam.get_s0())
  m2 = matrix.col(gonio.get_rotation_axis())

  # For all the miller indices
  for h in integrate_handle.hkl:
    h = matrix.col(h)

    # Calculate the angles
    angles = ra(h, ub)

    # For all the angles
    for phi in angles:
      r = m2.axis_and_angle_as_r3_rotation_matrix(angle=phi)
      pstar = r * ub * h
      s1 = s0 + pstar
      assert(abs(s1.length() - s0.length()) < 1e-7)

  print "OK"

  # Create a dict of lists of xy for each hkl
  gen_phi = {}
  for h in integrate_handle.hkl:

    # Calculate the angles
    angles = ra(h, ub)
    gen_phi[h] = angles
#        for phi in angles:
#            try:
#                a = gen_phi[h]
#                a.append(phi)
#                gen_phi[h] = a
#            except KeyError:
#                gen_phi[h] = [phi]

  # For each hkl in the xds file
  for hkl, xyz in zip(integrate_handle.hkl,
                      integrate_handle.xyzcal):

    # Calculate the XDS phi value
    xds_phi = scan.get_oscillation(deg=False)[0] + \
              xyz[2]*scan.get_oscillation(deg=False)[1]

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

    # Check the Phi values are the same
    assert(abs(xds_phi - my_phi) < 0.1)

  # Test Passed
  print "OK"

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
