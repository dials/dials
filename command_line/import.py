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
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dxtbx.model.experiment.experiment_list import ExperimentListDumper

help_message = '''

This program is used to import image data files into a format that can be used
within dials. The program looks at the metadata for each image along with the
filenames to determine the relationship between sets of images. Once all the
images have been analysed, a datablock object is written to file which specifies
the relationship between files. For example if two sets of images which belong
to two rotation scans have been given, two image sweeps will be saved. Images to
be processed are specified as command line arguments. Sometimes, there is a
maximum number of arguments that can be given on the command line and the number
of files may exceed this. In this case image filenames can be input on stdin
delmited by a new line using the -i option (see below for examples).

Examples:

  dials.import image_*.cbf

  dials.import image_1_*.cbf image_2_*.cbf

  find . -name "image_*.cbf" | dials.import -i

'''


class ImageFileImporter(object):
  ''' Import a data block from files. '''

  def __init__(self, options):
    ''' Initialise with the options'''
    self.options = options

  def __call__(self, args):
    import sys

    # Check we have some filenames
    if len(args) == 0 and not self.options.stdin:
      self.parser.print_help()
      exit(0)

    # Try reading from stdin as well
    if self.options.stdin:
      args.extend([l.strip() for l in sys.stdin.readlines()])

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


class ParseOptions(object):
  ''' Class to parse the command line options. '''
  from optparse import OptionParser

  class Parser(OptionParser):
    def format_epilog(self, formatter):
      return self.epilog

  def __init__(self):
    ''' Set the expected options. '''

    usage = "usage: %prog [options] /path/to/image/files"
    self.parser = ParseOptions.Parser(usage, epilog=help_message)

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
      type = "string", default = "datablock.json",
      help = "The output JSON or pickle file (filename.json | filename.pickle)")

    # Standard input read files, rather than from the command-line
    self.parser.add_option(
      "-i", "--stdin",
      dest = "stdin",
      action = "store_true",
      default = False,
      help = "Read filenames from standard input rather than command-line")

    # Write the datablock to JSON or Pickle
    self.parser.add_option(
      "-c", "--compact",
      dest = "compact",
      action = "store_true", default = False,
      help = "For JSON output, use compact representation")

    # Don't sort input filenames
    self.parser.add_option(
      "-n", "--no-sort",
      dest = "sort",
      action = "store_false", default = True,
      help = "Don't sort input files (default is True)")

  def __call__(self):
    ''' Parse the options. '''
    # Parse the command line arguments
    (options, args) = self.parser.parse_args()

    # Return the options and args
    return options, args


if __name__ == '__main__':

  # Parse the command line options
  parse = ParseOptions()
  (options, args) = parse()

  # Choose the importer to use
  importer = ImageFileImporter(options)
  importer.parser = parse.parser

  # Import the data
  importer(args)
