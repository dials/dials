from __future__ import division

import iotbx.phil
master_phil_scope = iotbx.phil.parse("""
unit_cell = None
  .type = unit_cell
space_group = None
  .type = space_group
d_min = None
  .type = float(value_min=0)
n_samples = 1000
  .type = int(value_min=0)
random_seed = 42
  .type = int
plot = False
  .type = bool
""")


def random_rotation():
  import random
  from scitbx.math import euler_angles_as_matrix
  return euler_angles_as_matrix([random.uniform(0,360) for i in xrange(3)])


def run(args):

  from libtbx.phil import command_line

  from dials.util.command_line import Importer
  from scitbx.array_family import flex
  importer = Importer(args)
  sweeps = importer.imagesets
  assert len(sweeps) == 1
  sweep = sweeps[0]

  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=importer.unhandled_arguments)
  working_phil.show()

  params = working_phil.extract()
  assert params.unit_cell is not None
  assert params.space_group is not None
  unit_cell = params.unit_cell
  space_group = params.space_group.group()

  import random
  from cctbx.crystal.crystal_model import crystal_model
  from cctbx import crystal, miller
  from scitbx.array_family import flex
  from scitbx import matrix

  flex.set_random_seed(params.random_seed)
  random.seed(params.random_seed)

  crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                      space_group=space_group)

  # the reciprocal matrix
  B = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()

  n_predicted = flex.double()

  for i in range(params.n_samples):
    U = random_rotation()
    A = U * B
    direct_matrix = A.inverse()

    cryst_model = crystal_model(direct_matrix[0:3],
                                direct_matrix[3:6],
                                direct_matrix[6:9],
                                space_group=space_group)

    from dials.algorithms.integration import ReflectionPredictor
    predictor = ReflectionPredictor()
    reflections = predictor(sweep, cryst_model)
    predicted_reflections = reflections
    miller_indices = predicted_reflections.miller_index()
    miller_set = miller.set(
      crystal_symmetry, miller_indices, anomalous_flag=True)
    if params.d_min is not None:
      resolution_sel = miller_set.d_spacings().data() > params.d_min
      predicted_reflections = predicted_reflections.select(resolution_sel)
    n_predicted.append(predicted_reflections.size())
    print n_predicted[-1]

  print "Basic statistics:"
  from scitbx.math import basic_statistics
  stats = basic_statistics(n_predicted)
  stats.show()

  print "Histogram:"
  hist = flex.histogram(n_predicted, n_slots=20)
  hist.show()

  print "Raw spot counts:"
  print list(n_predicted)

  if params.plot:
    from matplotlib import pyplot
    pyplot.bar(hist.slot_centers(), hist.slots(), width=hist.slot_width())
    pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
