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

  def __call__(self, filename, remove_invalid=False, **kwargs):
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
    if miller_index and remove_invalid:
      flags = flex.bool([h != (0, 0, 0) for h in table['hkl']])
      table = table.select(flags)
    Command.end('Removed invalid reflections, %d remaining' % len(table))

    # Output the table to pickle file
    if filename is None: filename = 'spot_xds.pickle'
    Command.start('Saving reflection table to %s' % filename)
    table.as_pickle(filename)
    Command.end('Saved reflection table to %s' % filename)


class IntegrateHKLImporter(object):
  ''' Class to import an integrate.hkl file to a reflection table. '''

  def __init__(self, integrate_hkl):
    self._integrate_hkl = integrate_hkl

  def __call__(self, filename, **kwargs):
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
    if filename is None: filename = 'integrate_hkl.pickle'
    Command.start('Saving reflection table to %s' % filename)
    table.as_pickle(filename)
    Command.end('Saved reflection table to %s' % filename)


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

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = None,
    help = "The output file")

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
  importer = select_importer(args)

  # Import the XDS data
  importer(options.output, remove_invalid=options.remove_invalid)
