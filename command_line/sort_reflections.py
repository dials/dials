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
from operator import attrgetter
from dials.util.script import ScriptRunner

class Sort(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = """usage: %prog [options] reflections.pickle

Example: %prog -k "miller_index:2,1,0" -o sorted.pickle"""

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output filename option
    self.config().add_option(
        '-o', '--output-reflections-filename',
        dest = 'output_reflections_filename',
        type = 'string', default = None,
        help = 'Set the filename for reflections predicted by the refined '
               'model.')

    # Reverse sort option
    self.config().add_option(
        '-r', '--reverse',
        dest = 'reverse_sort',
        action = "store_true",
        default = False,
        help = 'Whether to reverse the sort direction.')

    # Sort key option
    self.config().add_option(
        '-k', '--key',
        dest = 'key',
        type = 'string', default = 'hkl',
        help = ('The chosen sort key. This should be an attribute of '
                'the Reflection object. Array-like attributes may have a sort '
                'order specified, using this format "miller_index:0,1,2"'))

    # Verbosity option
    self.config().add_option(
        '-v', '--verbose',
        dest = 'verbose',
        action = "store_true",
        default = False,
        help = 'Whether to print any output.')

  @staticmethod
  def sort_refs_by_attr(reflections, attr, reverse=False):
    return sorted(reflections, key=attrgetter(attr), reverse=reverse)

  @staticmethod
  def sort_refs_by_attr_elems(reflections, attr, order, reverse=False):
    refs_sorted = reflections
    for i in order:
      refs_sorted = sorted(refs_sorted,
                           key=lambda x: attrgetter(attr)(x)[i],
                           reverse=reverse)
    return refs_sorted

  def main(self, params, options, args):
    '''Execute the script.'''
    import string
    from dials.model.serialize import load
    from cctbx.array_family import flex # import dependency
    import cPickle as pickle

    # Check the number of arguments is correct
    if len(args) != 1:
      self.config().print_help()
      return

    reflections = pickle.load(open(args[0], 'rb'))

    # Process input string to get attribute and optional order
    (attr, _, order) = options.key.partition(":")
    order = order.split(",")
    try:
      order = map(int, order)
    except ValueError:
      order = None

    # Is the requested attribute valid?
    try:
      val = attrgetter(attr)(reflections[0])
    except AttributeError as e:
      print "Please check your sort key."
      raise

    # Try to get a default order, e.g. (0,1,2) for array-like attributes
    if order is None:
      try:
        order = range(len(val))
      except TypeError:
        pass

    if order:
      if options.verbose:
        if options.reverse_sort:
          msg = "Sorting by " + attr + " elements " + str(order) + \
                " with sorts descending"
        else:
          msg = "Sorting by " + attr + " by elements " + str(order) + \
                " with sorts ascending"
        print msg
      rlist = Sort.sort_refs_by_attr_elems(reflections, attr, order,
                                           options.reverse_sort)
    else:
      if options.verbose:
        if options.reverse_sort:
          msg = "Sorting by " + attr + " with sorts descending"
        else:
          msg = "Sorting by " + attr + " with sorts ascending"
        print msg
      rlist = Sort.sort_refs_by_attr(reflections, attr, options.reverse_sort)

    if options.verbose:
      print "Head of sorted list " + attr + ":"
      n = min(len(rlist), 10)
      for i in range(10):
        print attrgetter(attr)(rlist[i])

    # Save sorted reflections to file
    output_reflections_filename = options.output_reflections_filename
    if output_reflections_filename:
      print "Saving reflections to {0}".format(output_reflections_filename)
      pickle.dump(rlist, open(output_reflections_filename, 'wb'),
          pickle.HIGHEST_PROTOCOL)

    return

if __name__ == '__main__':
  script = Sort()
  script.run()
