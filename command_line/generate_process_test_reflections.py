# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_process_test_reflections

from __future__ import division

def main(args):
  from dials.algorithms.simulation.generate_test_reflections import main \
    as _main
  from dials.algorithms.simulation.generate_test_reflections import \
    master_phil
  from libtbx.phil import command_line
  cmd = command_line.argument_interpreter(master_params = master_phil)
  working_phil = cmd.process_and_fetch(args = sys.argv[1:])
  _main(working_phil.extract())

if __name__ == '__main__':
  import sys
  main(sys.argv[1:])
