#!/usr/bin/env python
#
# dials.two_theta_offset.py
#
#  Copyright (C) 2016 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from libtbx.phil import parse

help_message = '''

dials.two_theta_offset experiment_one.json experiment_two.json

'''

phil_scope = parse('''
''', process_includes=True)

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] experiment_one.json experiment_two.json" \
            % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      check_format=False,
      read_experiments=True)

  def run(self):
    '''Execute the script.'''
    from dials.util.command_line import Command
    from dials.array_family import flex
    from dials.util.options import flatten_experiments

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    # Check the number of experiments
    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) < 2:
      self.parser.print_help()
      return

    detectors = [experiment.detector[0] for experiment in experiments]

    from scitbx import matrix

    # FIXME make this a phil parameter in x, y
    offset = 100

    # FIXME iterate over pairs of experiments, if close in two-theta discard

    # pick two positions, at nominal origin offset in fast, slow

    x1 = matrix.col(detectors[0].get_origin())
    y1 = matrix.col(detectors[0].get_origin()) + \
      offset * matrix.col(detectors[0].get_fast_axis()) + \
      offset * matrix.col(detectors[0].get_slow_axis())

    x2 = matrix.col(detectors[1].get_origin())
    y2 = matrix.col(detectors[1].get_origin()) + \
      offset * matrix.col(detectors[1].get_fast_axis()) + \
      offset * matrix.col(detectors[1].get_slow_axis())

    centre, axis = find_centre_of_rotation(x1, x2, y1, y2)
    print 'Centre of axis: %7.4f %7.4f %7.4f' % centre.elems, \
      '  axis:  %7.4f %7.4f %7.4f' % axis.elems

def component(a, n):
  return a - a.dot(n) * n

def find_centre_of_rotation(x1, x2, y1, y2):
  '''Find centre of rotation which takes postion x1 -> x2 and y1 -> y2'''

  # chords of rotation of x, y

  cx = x2 - x1
  cy = y2 - y1

  # know axis is perpendicular to both of these -> is cross product

  ncx = cx.normalize()
  ncy = cy.normalize()

  axis = cx.cross(cy).normalize()

  # normal vectors

  nx = component(cx, axis).normalize().cross(axis)
  ny = component(cy, axis).normalize().cross(axis)

  # origin of normal vectors

  ox = component(x1 + 0.5 * cx, axis)
  oy = component(y1 + 0.5 * cy, axis)

  import math

  h = (oy - ox).dot(ncx)
  a = (-ny).angle(ncx)
  d = h / math.cos(a)
  return oy + d * ny, axis

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
