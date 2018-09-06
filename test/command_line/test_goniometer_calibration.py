from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run

def test_goniometer_calibration(dials_regression):
  data_dir = os.path.join(dials_regression, 'rotation_calibration')
  o0_k0_p0 = os.path.join(data_dir, 'experiments_o0_k0_p0.json')
  o0_k0_p48 = os.path.join(data_dir, 'experiments_o0_k0_p48.json')
  o0_k48_p48 = os.path.join(data_dir, 'experiments_o0_k48_p48.json')
  o48_k48_p48 = os.path.join(data_dir, 'experiments_o48_k48_p48.json')

  args = ["dials.goniometer_calibration",
          o0_k0_p0, o0_k0_p48, o0_k48_p48, o48_k48_p48,
          "space_group=P422",
          "xoalign=xoalign_config.py"]

  command = " ".join(args)
  print(command)
  result = easy_run.fully_buffered(command=command).raise_if_errors()

  expected_output = '''
Goniometer axes and angles (ImgCIF coordinate system):
GON_PHI:  rotation of 47.976 degrees about axis (0.99997,-0.00588,-0.00586)
GON_KAPPA:  rotation of 47.773 degrees about axis (0.91314,0.27943,-0.29681)
GON_OMEGA:  rotation of 48.000 degrees about axis (1.00000,0.00000,0.00000)

Goniometer axes and angles (MOSFLM coordinate system):
GON_PHI:  rotation of 47.976 degrees about axis (0.00494,-0.00588,0.99997)
GON_KAPPA:  rotation of 47.773 degrees about axis (0.29597,0.27943,0.91341)
GON_OMEGA:  rotation of 48.000 degrees about axis (-0.00092,0.00000,1.00000)

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
  GON_OMEGA  rotation  goniometer  .          1.0000   0.0000   0.0000  .  .  .
  GON_KAPPA  rotation  goniometer  GON_OMEGA  0.9131   0.2794  -0.2968  .  .  .
  GON_PHI    rotation  goniometer  GON_KAPPA  1.0000  -0.0059  -0.0059  .  .  .
'''
  assert '\n'.join(result.stdout_lines[12:]) == ''.join(expected_output)

  assert os.path.exists('xoalign_config.py')
  expected_xoalign_config = '''\
GONIOMETER_AXES_NAMES = ('GON_OMEGA', 'GON_KAPPA', 'GON_PHI')
GONIOMETER_AXES = [(-0.00092, 0.00000, 1.00000), (0.29597, 0.27943, 0.91341), (0.00494, -0.00588, 0.99997)]
GONIOMETER_DATUM = (0,0,0) # in degrees
'''
  with open('xoalign_config.py', 'rb') as f:
    text = f.read()
    assert text == expected_xoalign_config
