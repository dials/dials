# LIBTBX_SET_DISPATCHER_NAME dev.dials.griddex
from __future__ import absolute_import, division, print_function

import libtbx.phil
import libtbx.load_env

help_message = '''

Cross reference indexing solutions.

Examples::

  %s expts0.json refl0.json

''' % libtbx.env.dispatcher_name

phil_scope = libtbx.phil.parse("""
  d_min = None
    .type = float(value_min=0.0)
""")

def test_index(experiment, reflections):

  # map reflections to reciprocal space from image space

  reflections.centroid_px_to_mm(experiment.detector, experiment.scan)

  reflections.map_centroids_to_reciprocal_space(
    experiment.detector, experiment.beam, experiment.goniometer)

  # now compute fractional indices - in Python rather than trying to push
  # everything to C++ for the moment

  from scitbx import matrix
  ub = matrix.sqr(experiment.crystal.get_A())
  rub = ub.inverse()

  from dials.array_family import flex
  hkl_real = flex.vec3_double(len(reflections))

  for j, rlp in enumerate(reflections['rlp']):
    hkl_real[j] = rub * rlp

  hkl = hkl_real.iround()

  ms = 0.0
  for (_h, _k, _l), (_hr, _kr, _lr) in zip(hkl, hkl_real):
    ms += (_hr - _h) ** 2 + (_kr - _k) ** 2 + (_lr - _l) ** 2

  import math
  return math.sqrt(ms / len(reflections))

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] experiments.json reflections.pickle" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  assert len(experiments) == len(reflections)

  nn = len(experiments)

  # FIXME check that all the crystals are in the primitive setting...

  # now compute grid of reciprocal RMSD's
  result = { }

  for j, expt in enumerate(experiments):
    for k, refl in enumerate(reflections):
      result[j, k] = test_index(expt, refl)

  # print matrix of results
  print('        ' + ''.join(['%7d' % j for j in range(nn)]))
  for k in range(nn):
    record = ''.join([' %6.3f' % result[j, k] for j in range(nn)])
    print('%8d' % k + record)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
