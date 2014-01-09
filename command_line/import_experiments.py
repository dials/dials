from __future__ import division
from dials.model.experiment.experiment_list import ExperimentListFactory
from dials.model.experiment.experiment_list import ExperimentListDumper

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

if __name__ == '__main__':
  import os.path
  from optparse import OptionParser

  usage = "usage: %prog [options] /path/to/image/files"
  parser = OptionParser(usage)

  # Print verbose output
  parser.add_option(
    "-v", "--verbose",
    dest = "verbose",
    action = "count", default = 0,
    help = "Set the verbosity level (-vv gives a verbosity level of 2)")

  # Write the experiment list to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = None,
    help = "The output JSON or pickle file (filename.json | filename.pickle)")

  # Write the experiment list to JSON or Pickle
  parser.add_option(
    "-c", "--compact",
    dest = "compact",
    action = "store_true", default = False,
    help = "For JSON output, use compact representation")

  # Write the models to different JSON files
  parser.add_option(
    "-s", "--split",
    dest = "split",
    action = "store_true", default = False,
    help = "For JSON output, split models into separate files")

  # Import from XDS directory
  parser.add_option(
    '--xds',
    dest = 'xds_dir',
    type = 'string', default = None,
    help = 'Directory containing XDS files.')

  # Specify the file to use
  parser.add_option(
    '--xds-file',
    dest = 'xds_file',
    type = 'string', default = None,
    help = 'Explicitly specify file to use (fname=xds_dir/xds_file)')

  # Parse the command line arguments
  (options, args) = parser.parse_args()

  # Either import from XDS or other file input
  if options.xds_dir is not None:
    xds_inp = os.path.join(options.xds_dir, 'XDS.INP')
    if options.xds_file is None:
      xds_file = find_best_xds_file(options.xds_dir)
    else:
      xds_file = os.path.join(options.xds_dir, options.xds_file)

    # Check a file is given
    if xds_file is None:
      raise RuntimeError('No XDS file found')

    # Load the experiment list
    unhandled = []
    experiments = ExperimentListFactory.from_xds(xds_inp, xds_file)

  else:

    # Check the XDS file option is not set
    if options.xds_file is not None:
      print 'Error! \'xds\' option must be set to use \'xds-file\' option'
      print ''
      parser.print_help()
      exit(0)

    # Check some files are given
    if len(args) == 0:
      parser.print_help()
      exit(0)

    # Get the experiments from the input files
    unhandled = []
    experiments = ExperimentListFactory.from_args(args,
      verbose=options.verbose, unhandled=unhandled)

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
