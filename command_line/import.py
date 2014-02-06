#!/usr/bin/env python
#
# import.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dxtbx.datablock import DataBlockFactory, DataBlockDumper
from dials.model.experiment.experiment_list import ExperimentListFactory
from dials.model.experiment.experiment_list import ExperimentListDumper


class ImageFileImporter(object):
  ''' Import a data block from files. '''

  def __init__(self, options):
    ''' Initialise with the options'''
    self.options = options

  def __call__(self, args):

    # Check we have some filenames
    if len(args) == 0:
      self.parser.print_help()
      exit(0)

    # Sort arguments
    if self.options.sort:
      args = sorted(args)

    # Get the data blocks from the input files
    # We've set verbose to print out files as they're tested.
    unhandled = []
    datablocks = DataBlockFactory.from_args(args,
      verbose=self.options.verbose, unhandled=unhandled)

    # Print out any unhandled files
    if len(unhandled) > 0:
      print '-' * 80
      print 'The following command line arguments were not handled:'
      for filename in unhandled:
        print '  %s' % filename

    # Loop through the data blocks
    for i, datablock in enumerate(datablocks):

      # Extract any sweeps
      sweeps = datablock.extract_sweeps()

      # Extract any stills
      stills = datablock.extract_stills()
      if not stills:
        num_stills = 0
      else:
        num_stills = len(stills)

      # Print some data block info
      print "-" * 80
      print "DataBlock %d" % i
      print "  format: %s" % str(datablock.format_class())
      print "  num images: %d" % datablock.num_images()
      print "  num sweeps: %d" % len(sweeps)
      print "  num stills: %d" % num_stills

      # Loop through all the sweeps
      if self.options.verbose > 1:
        for j, sweep in enumerate(sweeps):
          print ""
          print "Sweep %d" % j
          print "  length %d" % len(sweep)
          print sweep.get_beam()
          print sweep.get_goniometer()
          print sweep.get_detector()
          print sweep.get_scan()

    # Write the datablock to a JSON or pickle file
    if self.options.output:
      print "-" * 80
      print 'Writing datablocks to %s' % self.options.output
      dump = DataBlockDumper(datablocks)
      dump.as_file(self.options.output, compact=options.compact)


class XDSFileImporter(object):
  ''' Import a data block from xds. '''

  def __init__(self, options):
    ''' Initialise with the options'''
    self.options = options

  def __call__(self, args):
    import os
    # Get the XDS.INP file
    xds_inp = os.path.join(self.options.xds_dir, 'XDS.INP')
    if options.xds_file is None:
      xds_file = XDSFileImporter.find_best_xds_file(self.options.xds_dir)
    else:
      xds_file = os.path.join(self.options.xds_dir, self.options.xds_file)

    # Check a file is given
    if xds_file is None:
      raise RuntimeError('No XDS file found')

    # Load the experiment list
    unhandled = []
    experiments = ExperimentListFactory.from_xds(xds_inp, xds_file)

    # Print out any unhandled files
    if len(unhandled) > 0:
      print '-' * 80
      print 'The following command line arguments were not handled:'
      for filename in unhandled:
        print '  %s' % filename

    # Print some general info
    print '-' * 80
    print 'Read %d experiments' % len(experiments)

    # Loop through the data blocks
    for i, exp in enumerate(experiments):

      # Print some experiment info
      print "-" * 80
      print "Experiment %d" % i
      print "  format: %s" % str(exp.imageset.reader().get_format_class())
      print "  type: %s" % type(exp.imageset)
      print "  num images: %d" % len(exp.imageset)

      # Print some model info
      if options.verbose > 1:
        print ""
        if exp.beam:       print exp.beam
        else:              print "no beam!"
        if exp.detector:   print exp.detector
        else:              print "no detector!"
        if exp.goniometer: print exp.goniometer
        else:              print "no goniometer!"
        if exp.scan:       print exp.scan
        else:              print "no scan!"
        if exp.crystal:    print exp.crystal
        else:              print "no crystal!"

    # Write the experiment list to a JSON or pickle file
    if options.output:
      print "-" * 80
      print 'Writing experiments to %s' % options.output
      dump = ExperimentListDumper(experiments)
      dump.as_file(options.output, split=options.split, compact=options.compact)

    # Optionally save as a data block
    if options.xds_datablock:
      print "-" * 80
      print "Writing data block to %s" % options.xds_datablock
      dump = DataBlockDumper(experiments.to_datablocks())
      dump.as_file(options.xds_datablock, compact=options.compact)

  @staticmethod
  def find_best_xds_file(xds_dir):
    ''' Find the best available file.'''
    from os.path import exists, join

    # The possible files to check
    paths = [join(xds_dir, 'XDS_ASCII.HKL'),
             join(xds_dir, 'INTEGRATE.HKL'),
             join(xds_dir, 'GXPARM.XDS'),
             join(xds_dir, 'XPARM.XDS')]

    # Return the first path that exists
    for p in paths:
      if exists(p):
        return p

    # If no path exists, return None
    return None


class ParseOptions(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''

    from optparse import OptionParser
    usage = "usage: %prog [options] /path/to/image/files"
    self.parser = OptionParser(usage)

    # Print verbose output
    self.parser.add_option(
      "-v", "--verbose",
      dest = "verbose",
      action = "count", default = 0,
      help = "Set the verbosity level (-vv gives a verbosity level of 2)")

    # Write the datablock to JSON or Pickle
    self.parser.add_option(
      "-o", "--output",
      dest = "output",
      type = "string", default = None,
      help = "The output JSON or pickle file (filename.json | filename.pickle)")

    # Write the datablock to JSON or Pickle
    self.parser.add_option(
      "-c", "--compact",
      dest = "compact",
      action = "store_true", default = False,
      help = "For JSON output, use compact representation")

    # Write the models to different JSON files
    self.parser.add_option(
      "-s", "--split",
      dest = "split",
      action = "store_true", default = False,
      help = "For JSON output, split models into separate files")

    # Don't sort input filenames
    self.parser.add_option(
      "-n", "--no-sort",
      dest = "sort",
      action = "store_false", default = True,
      help = "Don't sort input files (default is True)")

    # Import from XDS directory
    self.parser.add_option(
      '--xds',
      dest = 'xds_dir',
      type = 'string', default = None,
      help = 'Directory containing XDS files.')

    # Specify the file to use
    self.parser.add_option(
      '--xds-file',
      dest = 'xds_file',
      type = 'string', default = None,
      help = 'Explicitly specify file to use (fname=xds_dir/xds_file)')

    # Add an option to output a datablock with xds as well.
    self.parser.add_option(
      '--xds-datablock',
      dest = 'xds_datablock',
      type = 'string', default = None,
      help = 'Output filename of data block with xds')

  def __call__(self):
    ''' Parse the options. '''
    # Parse the command line arguments
    (options, args) = self.parser.parse_args()

    # Check the options are consistent
    self.check_consistent(options)

    # Return the options and args
    return options, args

  def check_consistent(selfi, options):
    ''' Check consistency of options. '''
    ignored = []
    if options.xds_dir:
      pass
    else:
      if options.xds_file:
        ignored.append('--xds-file')
      if options.xds_datablock:
        ignored.append('--xds-datablock')
      if options.split:
        ignored.append('--split')
    if len(ignored) > 0:
      print ''
      print '-' * 80
      print 'The following command line options were ignored:'
      for line in ignored:
        print '  ' + line
      print '-' * 80
      print ''


if __name__ == '__main__':

  # Parse the command line options
  parse = ParseOptions()
  (options, args) = parse()

  # Choose the importer to use
  if options.xds_dir:
    importer = XDSFileImporter(options)
  else:
    importer = ImageFileImporter(options)

  # Import the data
  importer(args)
