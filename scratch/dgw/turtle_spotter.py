#!/usr/bin/env python
#
# turtle_spotter.py
#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

"""A turtle diffraction spot viewer"""

class ScriptRunner(object):
  """Class to run script."""

  def __init__(self, reflections_filename):
    """Setup the script."""

    # Filename data
    self.reflections_filename = reflections_filename

  def __call__(self):
    """Run the script."""
    import cPickle as pickle
    from dials.model.data import ReflectionList # import dependency
    from dials.util.command_line import Command

    # Read the pickle file
    Command.start('Reading reflection file.')
    with open(self.reflections_filename, 'rb') as f:
      self.reflections = pickle.load(f)

    Command.end('Read {0} spots from reflection file.'.format(
        len(self.reflections)))

    self.view()

  def view(self):
    import turtle
    coords = [ref.image_coord_mm for ref in self.reflections]
    x, y = zip(*coords)
    min_x, max_x = min(x), max(x)
    min_y, max_y = min(y), max(y)
    low = min(min_x, min_y)
    high = max(max_x, max_y)
    turtle.title("Reflections from " + self.reflections_filename)
    turtle.setworldcoordinates(low, low, high, high)
    turtle.pen(speed=0,pensize=2)
    turtle.hideturtle()
    for ref in self.reflections:
      (x, y) = ref.image_coord_mm
      turtle.penup()
      turtle.setx(x)
      turtle.sety(y)
      turtle.pendown()
      turtle.circle(1.0, steps=8)
    turtle.done()



if __name__ == '__main__':

  from optparse import OptionParser

  # Specify the command line options
  usage  = "usage: %prog [options] " \
           "/path/to/reflections.pickle "

  # Create an option parser
  parser = OptionParser(usage)

  # Parse the arguments
  options, args = parser.parse_args()

  # Print help if no arguments specified, otherwise call function
  if len(args) < 1:
    parser.print_help()
  else:
    # Initialise the script runner
    runner = ScriptRunner(
        reflections_filename=args[0])

    # Run the script
    runner()
