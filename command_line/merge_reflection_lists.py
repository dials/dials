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

class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    # Create the phil parameters
    phil_scope = parse('''

      output = merged.pickle
        .type = str
        .help = "The output file"

      method = *update extend
        .type = choice
        .help = "The method of merging"

    ''')

    # The script usage
    usage = "usage: %prog [options] /path/to/image/reflection/files"
    self.parser = OptionParser(usage=usage, phil=phil_scope)

  def run(self):
    ''' Run the script. '''
    from dials.array_family import flex
    from dials.util.command_line import Command

    # Parse the command line arguments
    params, options, args = self.parser.parse_args(show_diff_phil=True)
    if len(args) == 0:
      self.parser.print_help()
      exit(0)

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
    if params.method == "update":
      assert(all(n == nrows[0] for n in nrows[1:]))
      table = tables[0]
      for t in tables[1:]:
        table.update(t)
    elif params.method == "extend":
      assert(all(n == ncols[0] for n in ncols[1:]))
      table = tables[0]
      for t in tables[1:]:
        table.extend(t)
    else:
      raise RuntimeError('unknown method, %s' % params.method)

    # Write the reflections to the file
    Command.start('Writing %d reflections to %s' % (len(table), params.output))
    table.as_pickle(params.output)
    Command.end('Wrote %d reflections to %s' % (len(table), params.output))


if __name__  == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
