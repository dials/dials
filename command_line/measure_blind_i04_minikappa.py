#!/usr/bin/env python
# LIBTBX_SET_DISPATCHER_NAME dev.dials.measure_blind_i04_minikappa

from __future__ import absolute_import, division
from libtbx.phil import parse

help_message = '''

dev.dials.measure_blind_i04_minikappa wavelength=0.979 resolution=1.6

'''

phil_scope = parse('''
wavelength = 0.979
  .type = float
  .help = 'wavelength (for two-theta calculation)'
resolution = 1.6
  .type = float
  .help = 'resolution of diffraction'
''')

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] " \
            % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      check_format=False,
      read_experiments=False)

  def run(self):
    '''Execute the script.'''
    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    import math
    from scitbx import matrix
    from dials.algorithms.refinement import rotation_decomposition

    # FIXME get this from experiments.json
    e1 = matrix.col((1, 0, 0))
    e2 = matrix.col((0.914, 0.279, -0.297))
    e3 = matrix.col((1, 0, 0))

    # FIXME get this from experiments.json
    s0n = matrix.col((0, 0, -1))

    # rotate blind region about beam by +/- two theta
    two_theta = 2.0 * math.asin(0.5 * params.wavelength / params.resolution)
    R_ptt = s0n.axis_and_angle_as_r3_rotation_matrix(two_theta)
    R_ntt = s0n.axis_and_angle_as_r3_rotation_matrix(-two_theta)

    # now decompose to kappa, phi

    sol_plus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
      R_ptt, e1, e2, e3, return_both_solutions=True, deg=True)

    sol_minus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
      R_ntt, e1, e2, e3, return_both_solutions=True, deg=True)

    solutions = []
    if sol_plus:
      solutions.extend(sol_plus)
    if sol_minus:
      solutions.extend(sol_minus)

    print '%d solutions found, any one should do' % len(solutions)

    # FIXME work out extent to which this one will be shadowed... add **
    # if so...

    print '  Kappa     Phi'
    for s in solutions:
      print '%8.3f %8.3f' % (s[1], s[2])

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
