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
  from libtbx.phil import command_line
  from dials.util.command_line import Importer
  importer = Importer(args, check_format=False)
  assert importer.reflections is None
  assert len(importer.experiments) > 1

  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=importer.unhandled_arguments)
  working_phil.show()
  params = working_phil.extract()
  hkl = flex.miller_index(params.hkl)

  from dials.algorithms.indexing.compare_orientation_matrices import \
       show_rotation_matrix_differences
  show_rotation_matrix_differences(importer.experiments.crystals(),
                                   miller_indices=hkl)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
