from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME dev.dials.duelling_profiles

import iotbx.phil

phil_scope = iotbx.phil.parse("""\
  method = *example nonsense flat
    .type = choice
""", process_includes=True)

help_message = '''

Examples::

  dev.dials.duelling_profiles experiments.json integrated.pickle

'''

def model_background(shoebox, mean_bg):
  from scitbx.random import variate, poisson_distribution
  dz, dy, dx = shoebox.focus()
  g = variate(poisson_distribution(mean = mean_bg))
  for k in range(dz):
    for j in range(dy):
      for i in range(dx):
        shoebox[k, j, i] += g.next()
  return

def model_reflection_example(reflection, experiment):
  hkl = reflection['miller_index']
  i0 = reflection['intensity.sum.value'] / reflection['dqe']
  s1 = reflection['s1']
  xyz = reflection['xyzcal.px']
  pixels = reflection['shoebox']
  mean_bg = reflection['background.mean']
  crystal = experiment.crystal
  profile = experiment.profile
  Amat = crystal.get_A_at_scan_point(int(xyz[2]))
  return

def model_reflection_nonsense(reflection, experiment):
  return

def model_reflection_flat(reflection, experiment):
  pixels = reflection['shoebox']
  pixels.flatten()
  data = pixels.data
  dz, dy, dx = data.focus()
  return

def main(reflections, experiment, method):
  nref0 = len(reflections)

  if 'intensity.prf.variance' in reflections:
    selection = reflections.get_flags(
      reflections.flags.integrated,
      all=True)
  else:
    selection = reflections.get_flags(
      reflections.flags.integrated_sum)
  reflections = reflections.select(selection)

  nref1 = len(reflections)

  print 'Removed %d invalid reflections, %d remain' % (nref0 - nref1, nref1)

  for j, reflection in enumerate(reflections):
    globals()['model_reflection_%s' % method](reflection, experiment)

  return

def run(args):
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] integrated.pickle experiments.json" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) != 1 or len(reflections) != 1:
    parser.print_help()
    exit()

  if not 'shoebox' in reflections[0]:
    print 'Please add shoeboxes to reflection pickle'
    exit()

  main(reflections[0], experiments[0], params.method)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
