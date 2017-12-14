#!/usr/bin/env python
#
# show_integration_memory_info.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dev.dials.show_integration_memory_info

from __future__ import absolute_import, division
from libtbx.phil import parse

# Create the help message
help_message = '''

  Show some information about possible integration memory usage

'''

# Create the phil parameters
phil_scope = parse('''

  use_dynamic_mask = True
    .type = bool
    .help = "Use the dynamic mask"

''')


class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] /path/to/experiments.json" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      epilog=help_message,
      usage=usage,
      phil=phil_scope,
      read_experiments=True)

  def run(self):
    ''' Run the script. '''
    from dials.util.command_line import Command
    from dials.util.options import flatten_experiments
    from libtbx.utils import Sorry
    from libtbx.introspection import machine_memory_info
    from math import floor
    memory_info = machine_memory_info()
    total_memory = memory_info.memory_total()

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) == 0:
      self.parser.print_help()
      return
    elif len(experiments) > 1:
      raise RuntimeError("Only 1 experiment can be specified")

    imageset = experiments[0].imageset

    from dials.algorithms.integration.parallel_integrator import MultiThreadedIntegrator

    max_block_size = MultiThreadedIntegrator.compute_max_block_size(
      imageset,
      params.use_dynamic_mask,
      total_memory)

    print("Total memory:   %.2f GB" % (total_memory / 1e9))
    print("Max Block Size: %d" % max_block_size)

    for i in range(1, min(max_block_size, len(imageset)+1)):
      memory = MultiThreadedIntegrator.compute_required_memory(
        imageset[0:i],
        params.use_dynamic_mask)
      print("Block Size = %d, Memory = %.2f GB" % (i, memory / 1e9))

if __name__  == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
