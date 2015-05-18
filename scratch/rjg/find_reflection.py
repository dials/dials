from __future__ import division

from libtbx.phil import command_line
from libtbx.utils import Sorry
import iotbx.phil
# from dials.util.command_line import Importer
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_experiments
from dials.array_family import flex

help_message = '''

'''

phil_scope = iotbx.phil.parse("""
miller_index = None
  .type = ints(size=3)
  .multiple = True
expand_to_p1 = False
  .type = bool
""", process_includes=True)


def show_predictions(predictions):
  for p in predictions:
    print "Miller index: %8i %8i %8i" %p['miller_index']
    print "  xyzcal.mm:  %8.2f %8.2f %8.2f" %p['xyzcal.mm']
    print "  xyzcal.px:  %8.2f %8.2f %8.2f" %p['xyzcal.px']
    print


def run(args):
  import libtbx.load_env
  usage = "%s experiments.json [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  if len(experiments) == 0:
    parser.print_help()
    return
  elif len(experiments) > 1:
    raise Sorry("More than one experiment present")

  assert len(params.miller_index), "Must specify at least one miller_index to predict."

  experiment = experiments[0]

  reflections = flex.reflection_table()
  miller_indices = flex.miller_index()
  entering_flags = flex.bool()
  for mi in params.miller_index:
    miller_indices.append(mi)
    miller_indices.append(mi)
    entering_flags.append(True)
    entering_flags.append(False)
  reflections['miller_index'] = miller_indices
  reflections['entering'] = entering_flags
  reflections['id'] = flex.size_t(len(reflections), 0)

  if params.expand_to_p1:
    from cctbx.miller import expand_to_p1_iselection
    proxy = expand_to_p1_iselection(
      experiment.crystal.get_space_group(),
      anomalous_flag=True,
      indices=miller_indices,
      build_iselection=True)
    reflections = reflections.select(proxy.iselection)
    reflections['miller_index'] = proxy.indices

  from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictor
  predictor = ExperimentsPredictor([experiment])
  predicted = predictor.predict(reflections)

  zmin, zmax = experiment.scan.get_array_range()
  z = predicted['xyzcal.px'].parts()[2]
  predicted = predicted.select((z >= zmin) & (z <= zmax))

  show_predictions(predicted)



if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])

