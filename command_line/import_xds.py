#!/usr/bin/env python
#
# import_xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class SpotXDSImporter(object):
  ''' Class to import a spot.xds file to a reflection table. '''

  def __init__(self, spot_xds):
    self._spot_xds = spot_xds

  def __call__(self, options):
    ''' Import the spot.xds file. '''
    from iotbx.xds import spot_xds
    from dials.util.command_line import Command
    from dials.array_family import flex
    import dxtbx

    # Read the SPOT.XDS file
    Command.start('Reading SPOT.XDS')
    handle = spot_xds.reader()
    handle.read_file(self._spot_xds)
    centroid = handle.centroid
    intensity = handle.intensity
    try:
      miller_index = handle.miller_index
    except AttributeError:
      miller_index = None
    Command.end('Read {0} spots from SPOT.XDS file.'.format(len(centroid)))

    # Create the reflection list
    Command.start('Creating reflection list')
    table = flex.reflection_table()
    table['id'] = flex.size_t(len(centroid), 0)
    table['panel'] = flex.size_t(len(centroid), 0)
    if miller_index:
      table['hkl'] = flex.miller_index(miller_index)
    table['xyzobs.px.value'] = flex.vec3_double(centroid)
    table['intensity.raw.value'] = flex.double(intensity)
    Command.end('Created reflection list')

    # Remove invalid reflections
    Command.start('Removing invalid reflections')
    if miller_index and options.remove_invalid:
      flags = flex.bool([h != (0, 0, 0) for h in table['hkl']])
      table = table.select(flags)
    Command.end('Removed invalid reflections, %d remaining' % len(table))

    # Output the table to pickle file
    if options.output is None: options.output = 'spot_xds.pickle'
    Command.start('Saving reflection table to %s' % options.output)
    table.as_pickle(options.output)
    Command.end('Saved reflection table to %s' % options.output)


class IntegrateHKLImporter(object):
  ''' Class to import an integrate.hkl file to a reflection table. '''

  def __init__(self, integrate_hkl):
    self._integrate_hkl = integrate_hkl

  def __call__(self, options):
    ''' Import the integrate.hkl file. '''

    from iotbx.xds import integrate_hkl
    from dials.array_family import flex
    from dials.util.command_line import Command
    import dxtbx
    from math import pi
    from scitbx import matrix

    # Read the SPOT.XDS file
    Command.start('Reading INTEGRATE.HKL')
    handle = integrate_hkl.reader()
    handle.read_file(self._integrate_hkl)
    hkl    = handle.hkl
    xyzcal = handle.xyzcal
    xyzobs = handle.xyzobs
    iobs   = handle.iobs
    sigma  = handle.sigma
    Command.end('Read %d reflections from INTEGRATE.HKL file.' % len(hkl))

    # Create the reflection list
    Command.start('Creating reflection table')
    table = flex.reflection_table()
    table['id'] = flex.size_t(len(hkl), 0)
    table['panel'] = flex.size_t(len(hkl), 0)
    table['hkl'] = flex.miller_index(hkl)
    table['xyzcal.px'] = flex.vec3_double(xyzcal)
    table['xyzobs.px.value'] = flex.vec3_double(xyzobs)
    table['intensity.cor.value'] = flex.double(iobs)
    table['intensity.cor.variance'] = flex.double(sigma)**2
    Command.end('Created table with {0} reflections'.format(len(table)))

    # Output the table to pickle file
    if options.output is None: options.output = 'integrate_hkl.pickle'
    Command.start('Saving reflection table to %s' % options.output)
    table.as_pickle(options.output)
    Command.end('Saved reflection table to %s' % options.output)


class XDSFileImporter(object):
  ''' Import a data block from xds. '''

  def __init__(self, args):
    ''' Initialise with the options'''
    self.args = args

  def __call__(self, options):
    from dials.model.experiment.experiment_list import ExperimentListFactory
    from dials.model.experiment.experiment_list import ExperimentListDumper
    import os
    # Get the XDS.INP file
    xds_inp = os.path.join(self.args[0], 'XDS.INP')
    if options.xds_file is None:
      xds_file = XDSFileImporter.find_best_xds_file(self.args[0])
    else:
      xds_file = os.path.join(self.args[0], options.xds_file)

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
    if options.output is None:
      options.output = 'experiments.json'
    print "-" * 80
    print 'Writing experiments to %s' % options.output
    dump = ExperimentListDumper(experiments)
    dump.as_file(options.output)

    # Optionally save as a data block
    if options.xds_datablock:
      print "-" * 80
      print "Writing data block to %s" % options.xds_datablock
      dump = DataBlockDumper(experiments.to_datablocks())
      dump.as_file(options.xds_datablock)

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


def select_importer(args):
  from os.path import split
  import libtbx.load_env
  path, filename = split(args[0])
  if filename == 'SPOT.XDS':
    return SpotXDSImporter(args[0])
  elif filename == 'INTEGRATE.HKL':
    return IntegrateHKLImporter(args[0])
  else:
    raise RuntimeError('expected (SPOT.XDS|INTEGRATE.HKL), got %s' % filename)


if __name__ == '__main__':

  from optparse import OptionParser

  # The option parser
  usage = "usage: %prog [options] (SPOT.XDS|INTEGRATE.HKL)"
  parser = OptionParser(usage)

  # The thing to import
  parser.add_option(
    "-i", "--input",
    dest = "input",
    type = "choice", choices=["experiment", "reflections"],
    default = "experiment",
    help = "The input method")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = None,
    help = "The output file")

  # Specify the file to use
  parser.add_option(
    '--xds-file',
    dest = 'xds_file',
    type = 'string', default = None,
    help = 'Explicitly specify file to use (fname=xds_dir/xds_file)')

  # Add an option to output a datablock with xds as well.
  parser.add_option(
    '--xds-datablock',
    dest = 'xds_datablock',
    type = 'string', default = None,
    help = 'Output filename of data block with xds')

  # Remove invalid reflections
  parser.add_option(
    "-r", "--remove-invalid",
    dest = "remove_invalid",
    action = "store_true", default = False,
    help = "Remove non-index reflections (if miller indices are present)")

  # Print verbose output
  parser.add_option(
    "-v", "--verbose",
    dest = "verbose",
    action = "count", default = 0,
    help = "Set the verbosity level (-vv gives a verbosity level of 2)")

  # Parse the command line arguments
  (options, args) = parser.parse_args()

  # Check number of arguments
  if len(args) != 1:
    parser.print_help()
    exit(0)

  # Select the importer class
  if options.input == 'experiment':
    importer = XDSFileImporter(args)
    pass
  else:
    importer = select_importer(args)

  # Import the XDS data
  importer(options)
