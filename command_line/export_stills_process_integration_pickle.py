#!/usr/bin/env python
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
# LIBTBX_SET_DISPATCHER_NAME dev.dials.export_stills_process_integration_pickle

from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
import os
from xfel.command_line.frame_extractor import ConstructFrame
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
import sys

if __name__ == '__main__':

  experiments = ExperimentListFactory.from_json_file(sys.argv[1])
  reflections = flex.reflection_table.from_pickle(sys.argv[2])

  if len(sys.argv) == 4:
    outfile_template = sys.argv[3]
  else:
    outfile_template = "int-%03d.pickle"

  # Split everything into separate experiments for pickling
  for e_number in range(len(experiments)):
    experiment = experiments[e_number]
    e_selection = reflections['id'] == e_number
    subset = reflections.select(e_selection)

    frame = ConstructFrame(subset, experiment).make_frame()
    frame["pixel_size"] = experiment.detector[0].get_pixel_size()[0]

    outfile = outfile_template % (e_number+1)

    print("Writing %s" % outfile)
    easy_pickle.dump(outfile, frame)

