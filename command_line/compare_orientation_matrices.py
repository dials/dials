from __future__ import division
from cctbx.array_family import flex

import iotbx.phil

master_phil_scope = iotbx.phil.parse(
"""
hkl = None
  .type = ints(size=3)
  .multiple=True
""")


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from libtbx.utils import Sorry

  parser = OptionParser(
    phil=master_phil_scope,
    read_experiments=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  if len(experiments) <= 1:
    raise Sorry('more than 1 experiment is required')

  hkl = flex.miller_index(params.hkl)

  from dials.algorithms.indexing.compare_orientation_matrices import \
       show_rotation_matrix_differences
  show_rotation_matrix_differences(experiments.crystals(),
                                   miller_indices=hkl)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
