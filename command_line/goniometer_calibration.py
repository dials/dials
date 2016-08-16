from __future__ import division
from cctbx.array_family import flex

import iotbx.phil

help_message = '''

Examples::

  dials.goniometer_calibration experiments_1.json experiments_2.json axis.name=GONIO_PHI axis.name=GONIO_OMEGA

'''


phil_scope = iotbx.phil.parse(
'''
space_group = None
  .type = space_group
axis
  .multiple = True
{
  name = None
    .type = str
}
''')


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from libtbx.utils import Sorry
  import libtbx.load_env

  usage = "%s [options] experiments.json" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  if len(experiments) <= 1:
    parser.print_help()
    return

  assert len(params.axis) == len(experiments) - 1

  from dials.algorithms.indexing.compare_orientation_matrices import \
       difference_rotation_matrix_axis_angle
  from scitbx import matrix
  crystals = []
  for experiment in experiments:
    crystal = experiment.crystal
    gonio = experiment.goniometer
    scan = experiment.scan
    fixed_rotation = matrix.sqr(gonio.get_fixed_rotation())
    setting_rotation = matrix.sqr(gonio.get_setting_rotation())
    rotation_axis = matrix.col(gonio.get_rotation_axis())
    rotation_matrix = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
      scan.get_oscillation()[0], deg=True)
    U = matrix.sqr(crystal.get_U())
    U = setting_rotation * rotation_matrix * fixed_rotation * U
    crystal.set_U(U)
    if params.space_group is not None:
      crystal.set_space_group(params.space_group.group())
    crystals.append(crystal)

  rows = []

  for i in range(len(crystals) - 1):
    R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
      crystals[i], crystals[i+1])
    print "Rotation of %.3f degrees" %angle, "about axis (%.3f, %.3f, %.3f)" %axis
    depends_on = '.'
    if i+1 < len(params.axis):
      depends_on = params.axis[i+1].name
    rows.insert(0, (
      params.axis[i].name, 'rotation', 'goniometer', depends_on,
      '%.4f' %axis[0], '%.4f' %axis[1], '%.4f' %axis[2], 0, 0, 0))

  from iotbx import cif
  loop = cif.model.loop(
    header=['_axis.id', '_axis.type', '_axis.equipment', '_axis.depends_on',
            '_axis.vector[1]', '_axis.vector[2]', '_axis.vector[3]',
            '_axis.offset[1]', '_axis.offset[2]', '_axis.offset[3]'])
  for row in reversed(rows):
    loop.add_row(row)

  print loop

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
