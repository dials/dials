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
  from optparse import OptionParser
  from dials.util.command_line import Importer

  # Initialise the option parser
  usage = "usage: %prog [options] {datablock.json|experiments.json}"
  parser = OptionParser(usage)

  # Parse the command line arguments
  options, args = parser.parse_args()

  # Load the datablock or experiment list
  importer = Importer(args, include=['datablocks', 'experiments'])

  # Print unhandled arguments
  if len(importer.unhandled_arguments) > 0:
    print '-' * 80
    print 'The following arguments were not handled'
    for a in importer.unhandled_arguments:
      print '  %s' % a
    print '-' * 80

  # Get the experiments and datablocks
  datablocks = importer.datablocks
  experiments = importer.experiments
  if ((datablocks == None or len(datablocks) == 0) and
      (experiments == None or len(experiments) == 0)):
    parser.print_help()
    exit(0)

  # Loop through data block checking detectors
  if datablocks:
    for i, db in enumerate(datablocks):
      for j, imageset in enumerate(db.extract_imagesets()):
        beam = imageset.get_beam()
        detector = imageset.get_detector()
        d_min = detector.get_max_resolution(beam.get_s0())
        print 'Datablock %d, imageset %d: max resolution = %f' % (i, j, d_min)

  # Loop through the experments checking detectors
  if experiments:
    for i, ex in enumerate(experiments):
      beam = ex.beam
      detector = ex.detector
      d_min = detector.get_max_resolution(beam.get_s0())
      print 'Experiment %d: max resolution = %f' % (i, d_min)
