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
nproc = 1
  .type = int(value_min=1)
""")


def random_rotation():
  import random
  from scitbx.math import euler_angles_as_matrix
  return euler_angles_as_matrix([random.uniform(0,360) for i in xrange(3)])


def run(args):

  from libtbx.phil import command_line

  from dials.util.command_line import Importer
  from dials.array_family import flex
  print args
  importer = Importer(args, check_format=False)
  assert len(importer.datablocks) == 1
  sweeps = importer.datablocks[0].extract_imagesets()
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
  from scitbx import matrix

  flex.set_random_seed(params.random_seed)
  random.seed(params.random_seed)

  crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                      space_group=space_group)

  # the reciprocal matrix
  B = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()

  n_predicted = flex.double()

  def predict_once(args):
    from dials.model.experiment.experiment_list import Experiment
    U = args[0]
    A = U * B
    direct_matrix = A.inverse()
    cryst_model = crystal_model(direct_matrix[0:3],
                                direct_matrix[3:6],
                                direct_matrix[6:9],
                                space_group=space_group)
    experiment = Experiment(imageset=sweep,
                            beam=sweep.get_beam(),
                            detector=sweep.get_detector(),
                            goniometer=sweep.get_goniometer(),
                            scan=sweep.get_scan(),
                            crystal=cryst_model)
    predicted_reflections = flex.reflection_table.from_predictions(
      experiment)
    miller_indices = predicted_reflections['miller_index']
    miller_set = miller.set(
      crystal_symmetry, miller_indices, anomalous_flag=True)
    if params.d_min is not None:
      resolution_sel = miller_set.d_spacings().data() > params.d_min
      predicted_reflections = predicted_reflections.select(resolution_sel)
    return len(predicted_reflections)

  from libtbx import easy_mp
  args = [(random_rotation(),) for i in range(params.n_samples)]
  results = easy_mp.parallel_map(
    func=predict_once,
    iterable=args,
    processes=params.nproc,
    preserve_order=True,
    preserve_exception_message=True)
  n_predicted = flex.double(results)

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
    from matplotlib.backends.backend_pdf import PdfPages

    red, blue = '#B2182B', '#2166AC'
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.bar(hist.slot_centers(), hist.slots(), width=0.75*hist.slot_width(),
           color=blue, edgecolor=blue)
    ax.set_xlabel('Spot count')
    ax.set_ylabel('Frequency')
    pdf = PdfPages("predicted_count_histogram.pdf")
    pdf.savefig(fig)
    pdf.close()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
