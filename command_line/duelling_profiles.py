from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME dev.dials.duelling_profiles

import iotbx.phil

phil_scope = iotbx.phil.parse("""\
""", process_includes=True)

help_message = '''

Examples::

  dev.dials.duelling_profiles experiments.json integrated.pickle

'''

def model_reflection(reflection, crystal, profile):
  hkl = reflection['miller_index']
  i0 = reflection['intensity.sum.value'] / reflection['dqe']
  s1 = reflection['s1']
  xyz = reflection['xyzcal.px']
  pixels = reflection['shoebox']
  Amat = crystal.get_A_at_scan_point(int(xyz[2]))
  return

def main(reflections, experiments):
  nref0 = len(reflections)
  crystal = experiments.crystal
  profile = experiments.profile

  if 'intensity.prf.variance' in reflections:
    selection = reflections.get_flags(
      reflections.flags.integrated,
      all=True)
  else:
    selection = reflections.get_flags(
      reflections.flags.integrated_sum)
  reflections = reflections.select(selection)

  nref1 = len(reflections)

  print 'Removed %d invalid reflections, %d remain' % (nref0-nref1, nref1)

  for j, reflection in enumerate(reflections):
    model_reflection(reflection, crystal, profile)

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

  main(reflections[0], experiments[0])

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
