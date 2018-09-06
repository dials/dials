# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_test_reflections

from __future__ import absolute_import, division, print_function

def main(args):
  from dials.algorithms.simulation.generate_test_reflections import \
     simple_gaussian_spots
  from dials.algorithms.simulation.generate_test_reflections import \
    master_phil
  from libtbx.phil import command_line
  cmd = command_line.argument_interpreter(master_params = master_phil)
  working_phil = cmd.process_and_fetch(args = sys.argv[1:])
  params = working_phil.extract()
  rlist = simple_gaussian_spots(params)
  import six.moves.cPickle as pickle
  if params.output.all:
    with open(params.output.all, 'wb') as fh:
      pickle.dump(rlist, fh, pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
  import sys
  main(sys.argv[1:])
