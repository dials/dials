import pycbf
from scitbx import matrix
import libtbx.load_env
import os
import math

have_dials_regression = libtbx.env.has_module("dials_regression")
assert(have_dials_regression)

dials_regression = libtbx.env.find_in_repositories(
  relative_path="dials_regression",
  test=os.path.isdir)

test_file = os.path.join(dials_regression, 'image_examples', 'DLS_I19',
                         'I19_P300k_00001.cbf')

handle = pycbf.cbf_handle_struct()
handle.read_widefile(test_file, pycbf.MSG_DIGEST)
crystal = handle.get_crystal_id()
gonio = handle.construct_goniometer()

# go find gonio axis names

handle.rewind_datablock()
handle.find_category('diffrn_measurement_axis')
assert(handle.count_rows() == gonio.axes)

gonio_axes = []

for j in range(handle.count_rows()):
  handle.rewind_column()
  handle.select_row(j)
  handle.find_column('axis_id')
  gonio_axes.append(handle.get_value())

# now get the axis information for these

axis_table = { }
scan_axis = None

for axis in gonio_axes:
  start, delta = handle.get_axis_setting(axis)
  axis_no = handle.count_axis_ancestors(axis)
  axis_table[axis_no] = {'vector':matrix.col(handle.get_axis_vector(axis)),
                         'start':start, 'delta':delta, 'name':axis}
  if delta:
    assert(not scan_axis), (scan_axis, axis, delta)
    scan_axis = axis

# now compute the pre-rotation, scan-rotation, post-rotation - FIXME need to 
# make sure that orientation of axis is being correctly applied...

total_r_before = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
total_r_scan = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
total_r_after = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

d2r = math.pi / 180.0

scan = False

for j in sorted(axis_table):
  v = axis_table[j]['vector']
  s = d2r * axis_table[j]['start']
  if axis_table[j]['delta']:
    total_r_scan *= v.axis_and_angle_as_r3_rotation_matrix(s)
    scan = True
  elif scan:
    total_r_after *= v.axis_and_angle_as_r3_rotation_matrix(s)
  else:
    total_r_before *= v.axis_and_angle_as_r3_rotation_matrix(s)

print 'Before'
print total_r_before
print 'Scan'
print total_r_scan
print 'After'
print total_r_after

