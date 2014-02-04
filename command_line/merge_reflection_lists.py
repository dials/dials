#!/usr/bin/env python
#
# merge_reflection_lists.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

if __name__  == '__main__':

  from dials.array_family import flex
  from optparse import OptionParser
  from dials.util.command_line import Command

  # The script usage
  usage = "usage: %prog [options] /path/to/image/reflection/files"
  parser = OptionParser(usage)

  # Print verbose output
  parser.add_option(
    "-v", "--verbose",
    dest = "verbose",
    action = "count", default = 0,
    help = "Set the verbosity level (-vv gives a verbosity level of 2)")

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = 'merged.pickle',
    help = "The output file")

  # Add merging method option
  parser.add_option(
    "-m", "--method",
    dest = "method",
    type = "choice", choices=["update", "extend"], default="update",
    help = "The method of merging")

  # Parse the command line arguments
  (options, args) = parser.parse_args()

  # Load the reflection lists
  tables = []
  for filename in args:
    Command.start('Loading reflection list from %s' % filename)
    t = flex.reflection_table.from_pickle(filename)
    Command.end('Loaded %d reflections from %s' % (len(t), filename))
    tables.append(t)

  # Check the number of tables
  if len(tables) < 2:
    raise RuntimeError('Need atleast 2 tables to merge, got %d' % len(tables))

  # Get the number of rows and columns
  nrows = [t.nrows() for t in tables]
  ncols = [t.ncols() for t in tables]

  # Merge the reflection lists
  if options.method == "update":
    assert(all(n == nrows[0] for n in nrows[1:]))
    table = tables[0]
    for t in tables[1:]:
      table.update(t)
  elif options.method == "extend":
    assert(all(n == ncols[0] for n in ncols[1:]))
    table = tables[0]
    for t in tables[1:]:
      table.extend(t)
  else:
    raise RuntimeError('unknown method, %s' % options.method)

  # Write the reflections to the file
  Command.start('Writing %d reflections to %s' % (len(table), options.output))
  table.as_pickle(options.output)
  Command.end('Wrote %d reflections to %s' % (len(table), options.output))
