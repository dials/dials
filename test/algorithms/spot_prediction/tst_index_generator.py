from __future__ import division

import libtbx.load_env
from os.path import realpath, dirname, join, isdir
have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression", test=isdir)

def run():
  if not have_dials_regression:
    print "Skipping test: dials_regression not available."
    return

  from iotbx.xds import xparm, integrate_hkl
  from dials.util import ioutil
  from dials.algorithms.spot_prediction import IndexGenerator
  import numpy
  from rstbx.cftbx.coordinate_frame_converter import \
      coordinate_frame_converter
  from scitbx import matrix

  # The XDS files to read from
  integrate_filename = join(dials_regression, 'data/sim_mx/INTEGRATE.HKL')
  gxparm_filename = join(dials_regression, 'data/sim_mx/GXPARM.XDS')

  # Read the XDS files
  integrate_handle = integrate_hkl.reader()
  integrate_handle.read_file(integrate_filename)
  gxparm_handle = xparm.reader()
  gxparm_handle.read_file(gxparm_filename)

  # Get the parameters we need from the GXPARM file
  d_min = 1.6
  space_group_type = ioutil.get_space_group_type_from_xparm(gxparm_handle)
  cfc = coordinate_frame_converter(gxparm_filename)
  a_vec = cfc.get('real_space_a')
  b_vec = cfc.get('real_space_b')
  c_vec = cfc.get('real_space_c')
  unit_cell = cfc.get_unit_cell()
  UB = matrix.sqr(a_vec + b_vec + c_vec).inverse()
  ub_matrix = UB

  # Generate the indices
  index_generator = IndexGenerator(unit_cell, space_group_type, d_min)
  miller_indices = index_generator.to_array()

  # Get individual generated hkl
  gen_h = [hkl[0] for hkl in miller_indices]
  gen_k = [hkl[1] for hkl in miller_indices]
  gen_l = [hkl[2] for hkl in miller_indices]

  # Get individual xds generated hkl
  xds_h = [hkl[0] for hkl in integrate_handle.hkl]
  xds_k = [hkl[1] for hkl in integrate_handle.hkl]
  xds_l = [hkl[2] for hkl in integrate_handle.hkl]

  # Get min/max generated hkl
  min_gen_h, max_gen_h = numpy.min(gen_h), numpy.max(gen_h)
  min_gen_k, max_gen_k = numpy.min(gen_k), numpy.max(gen_k)
  min_gen_l, max_gen_l = numpy.min(gen_l), numpy.max(gen_l)

  # Get min/max xds generated hkl
  min_xds_h, max_xds_h = numpy.min(xds_h), numpy.max(xds_h)
  min_xds_k, max_xds_k = numpy.min(xds_k), numpy.max(xds_k)
  min_xds_l, max_xds_l = numpy.min(xds_l), numpy.max(xds_l)

  # Ensure we have the whole xds range  in the generated set
  assert(min_gen_h <= min_xds_h and max_gen_h >= max_xds_h)
  assert(min_gen_k <= min_xds_k and max_gen_k >= max_xds_k)
  assert(min_gen_l <= min_xds_l and max_gen_l >= max_xds_l)

  # Test Passed
  print "OK"

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
