#!/usr/bin/env dials.python
from __future__ import division

import sys, os

from libtbx.phil import command_line, parse

if len(sys.argv) != 2: exit("please pass the path to a phil file")
#with(open(sys.argv[1])) as f:
#  phil = f.read()

#print phil
phil = sys.argv[1]

master_phil = parse("""
  input
    .multiple = True
  {
    experiments = None
      .type = path
    reflections = None
      .type = path
  }
  """)

cmd_line = command_line.argument_interpreter(master_params=master_phil)
working_phil = cmd_line.process_and_fetch(args=(phil,))
working_params = working_phil.extract()
for input in working_params.input:
  print input.experiments, input.reflections
