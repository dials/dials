#!/usr/bin/env python
#
# dials.sort_reflections.py
#
#  Copyright (C) 2013 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Sort(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    phil_scope = parse('''

      key = 'miller_index'
        .type = str
        .help = "The chosen sort key. This should be an attribute of "
                "the Reflection object."

      reverse = False
        .type = bool
        .help = "Reverse the sort direction"

      output = None
        .type = str
        .help = "The output reflection filename"

    ''')

    # The script usage
    usage  = """
      usage: %prog [options] reflections.pickle

      Example: %prog key=miller_index output=sorted.pickle
    """

    # Initialise the base class
    self.parser = OptionParser(usage=usage, phil=phil_scope)

  def run(self):
    '''Execute the script.'''
    from dials.array_family import flex # import dependency

    # Parse the command line
    params, options, args = self.parser.parse_args(show_diff_phil=True)

    # Check the number of arguments is correct
    if len(args) != 1:
      self.parser.print_help()
      return

    # Load the reflections
    reflections = flex.reflection_table.from_pickle(args[0])

    # Check the key is valid
    assert(params.key in reflections)

    # Sort the reflections
    print "Sorting by %s with reverse=%r" % (params.key, params.reverse)
    reflections.sort(params.key, params.reverse)

    if options.verbose > 0:
      print "Head of sorted list " + attr + ":"
      n = min(len(reflections), 10)
      for i in range(10):
        print (reflections[i][attr])

    # Save sorted reflections to file
    if params.output:
      print "Saving reflections to {0}".format(params.output)
      reflections.as_pickle(params.output)

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Sort()
    script.run()
  except Exception as e:
    halraiser(e)
