from __future__ import division
from dials.model.experiment.experiment_list import ExperimentListFactory
from dials.model.experiment.experiment_list import ExperimentListDumper

if __name__ == '__main__':

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

  # Parse the command line arguments
  (options, args) = parser.parse_args()
  if len(args) == 0:
    parser.print_help()

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
