#!/usr/bin/env python
#
# detector_max_resolution.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

if __name__ == '__main__':
  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks, flatten_experiments

  # Initialise the option parser
  usage = "usage: %prog [options] {datablock.json|experiments.json}"
  parser = OptionParser(
    usage=usage,
    read_datablocks=True,
    read_experiments=True)

  # Parse the command line arguments
  params, options = parser.parse_args(show_diff_phil=True)

  # Get the experiments and datablocks
  if (len(params.input.datablock) == 0 and len(params.input.experiments) == 0):
    parser.print_help()
    exit(0)

  # Extract into a list
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)

  # Loop through data block checking detectors
  for i, db in enumerate(datablocks):
    for j, imageset in enumerate(db.extract_imagesets()):
      beam = imageset.get_beam()
      detector = imageset.get_detector()
      d_min = detector.get_max_resolution(beam.get_s0())
      print 'Datablock %d, imageset %d: max resolution = %f' % (i, j, d_min)

  # Loop through the experments checking detectors
  for i, ex in enumerate(experiments):
    beam = ex.beam
    detector = ex.detector
    d_min = detector.get_max_resolution(beam.get_s0())
    print 'Experiment %d: max resolution = %f' % (i, d_min)
