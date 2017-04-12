from __future__ import absolute_import, division
import os
import libtbx.load_env
from libtbx import easy_run

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)


def run(args):
  if not have_dials_regression:
    print "Skipping tst_goniometer_calibration.py: dials_regression not available."
    return

  data_dir = os.path.join(dials_regression, 'rotation_calibration')
  o0_k0_p0 = os.path.join(data_dir, 'experiments_o0_k0_p0.json')
  o0_k0_p48 = os.path.join(data_dir, 'experiments_o0_k0_p48.json')
  o0_k48_p48 = os.path.join(data_dir, 'experiments_o0_k48_p48.json')
  o48_k48_p48 = os.path.join(data_dir, 'experiments_o48_k48_p48.json')

  args = ["dials.goniometer_calibration",
          o0_k0_p0, o0_k0_p48, o0_k48_p48, o48_k48_p48,
          "xoalign=xoalign_config.py"]

  command = " ".join(args)
  print command
  result = easy_run.fully_buffered(command=command).raise_if_errors()

  expected_output = '''
Goniometer axes and angles (ImgCIF coordinate system):
GON_PHI:  rotation of 157.741 degrees about axis (-0.44164,-0.59119,0.67487)
GON_KAPPA:  rotation of 143.755 degrees about axis (-0.51623,-0.85046,0.10111)
GON_OMEGA:  rotation of 48.006 degrees about axis (1.00000,0.00002,-0.00049)

Goniometer axes and angles (MOSFLM coordinate system):
GON_PHI:  rotation of 157.741 degrees about axis (-0.67446,-0.59119,-0.44226)
GON_KAPPA:  rotation of 143.755 degrees about axis (-0.10064,-0.85046,-0.51633)
GON_OMEGA:  rotation of 48.006 degrees about axis (-0.00043,0.00002,1.00000)

ImgCIF _axis loop template:
loop_
  _axis.id
  _axis.type
  _axis.equipment
  _axis.depends_on
  _axis.vector[1]
  _axis.vector[2]
  _axis.vector[3]
  _axis.offset[1]
  _axis.offset[2]
  _axis.offset[3]
  GON_PHI    rotation  goniometer  GON_KAPPA  -0.4416  -0.5912   0.6749  .  .  .
  GON_KAPPA  rotation  goniometer  GON_OMEGA  -0.5162  -0.8505   0.1011  .  .  .
  GON_OMEGA  rotation  goniometer  .           1.0000   0.0000  -0.0005  .  .  .
'''
  for line in expected_output.splitlines():
    if not line: continue
    assert line in result.stdout_lines, line

  assert os.path.exists('xoalign_config.py')
  expected_xoalign_config = '''\
GONIOMETER_AXES_NAMES = ('GON_OMEGA', 'GON_KAPPA', 'GON_PHI')
GONIOMETER_AXES = [(-0.00043, 0.00002, 1.00000), (-0.10064, -0.85046, -0.51633), (-0.67446, -0.59119, -0.44226)]
GONIOMETER_DATUM = (0,0,0) # in degrees
'''
  with open('xoalign_config.py', 'rb') as f:
    text = f.read()
    from libtbx.test_utils import show_diff
    assert not show_diff(text, expected_xoalign_config)

  print "OK"

if __name__ == '__main__':
  import sys
  from dials.test import cd_auto
  with cd_auto(__file__):
    run(sys.argv[1:])
