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

  def __init__(self, params, options):
    ''' Initialise with the options'''
    self.params = params
    self.options = options

  def __call__(self, args):
    import sys

    # Sort arguments
    if self.params.sort:
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
    if self.params.output:
      print "-" * 80
      print 'Writing datablocks to %s' % self.params.output
      dump = DataBlockDumper(datablocks)
      dump.as_file(self.params.output, compact=self.params.compact)


class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # Create the phil parameters
    phil_scope = parse('''

      output = datablock.json
        .type = str
        .help = "The output JSON or pickle file"

      compact = False
        .type = bool
        .help = "For JSON output use compact representation"

      sort = False
        .type = bool
        .help = "Sort input files"

    ''')

    # Create the option parser
    usage = "usage: %prog [options] /path/to/image/files"
    self.parser = OptionParser(
      usage=usage,
      stdin_options=True,
      phil=phil_scope,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''

    # Parse the command line arguments
    params, options, args = self.parser.parse_args(
      show_diff_phil=True,
      return_unhandled=True)

    # Check we have some filenames
    if len(args) == 0:
      self.parser.print_help()
      exit(0)

    # Choose the importer to use
    importer = ImageFileImporter(params, options)

    # Import the data
    importer(args)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
