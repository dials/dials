#!/usr/bin/env python
# LIBTBX_SET_DISPATCHER_NAME dev.dials.measure_blind_i04_minikappa

from __future__ import absolute_import, division
from libtbx.phil import parse

help_message = '''

dev.dials.measure_blind_i04_minikappa resolution=1.6 experiments.json

'''

phil_scope = parse('''
resolution = 0.0
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
      check_format=True,
      read_experiments=True)

  def run(self):
    '''Execute the script.'''
    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    import math
    from scitbx import matrix
    from dials.algorithms.refinement import rotation_decomposition
    from dials.util.options import flatten_experiments
    from dials.array_family import flex

    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) != 1:
      self.parser.print_help()
      return

    expt = experiments[0]

    axes = expt.goniometer.get_axes()
    beam = expt.beam
    det = expt.detector

    if params.resolution:
      resolution = params.resolution
    else:
      resolution = det.get_max_inscribed_resolution(expt.beam.get_s0())

    # at this point, predict all of the reflections in the scan possible (i.e.
    # extend scan to 360 degrees) - this points back to expt

    scan = self.make_scan_360(expt.scan)

    # now get a full set of all unique miller indices
    from cctbx import miller
    from cctbx import crystal

    all = miller.build_set(crystal_symmetry=crystal.symmetry(
      space_group=expt.crystal.get_space_group(),
      unit_cell=expt.crystal.get_unit_cell()), anomalous_flag=True,
      d_min=resolution)

    obs = self.predict_to_miller_set_with_shadow(expt, resolution)

    print 'Fraction of unique observations at datum: %.3f' % \
      (len(obs.indices()) / len(all.indices()))

    missing = all.lone_set(other=obs)

    print '%d unique reflections in blind region' % len(missing.indices())

    e1 = matrix.col(axes[0])
    e2 = matrix.col(axes[1])
    e3 = matrix.col(axes[2])

    s0n = matrix.col(beam.get_s0()).normalize()

    # rotate blind region about beam by +/- two theta
    two_theta = 2.0 * math.asin(0.5 * beam.get_wavelength() / resolution)
    R_ptt = s0n.axis_and_angle_as_r3_rotation_matrix(two_theta)
    R_ntt = s0n.axis_and_angle_as_r3_rotation_matrix(-two_theta)

    # now decompose to kappa, phi
    sol_plus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
      R_ptt, e1, e2, e3, return_both_solutions=True, deg=False)

    sol_minus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
      R_ntt, e1, e2, e3, return_both_solutions=True, deg=False)

    solutions = []
    if sol_plus:
      solutions.extend(sol_plus)
    if sol_minus:
      solutions.extend(sol_minus)

    if not solutions:
      print 'Maximum two theta: %.3f,' % (two_theta * 180.0 / math.pi),
      print 'sorry, this is impossible with this goniometer'
      return

    print 'Maximum two theta: %.3f,' % (two_theta * 180.0 / math.pi),
    print '%d solutions found' % len(solutions)

    print '  Kappa     Phi    #new'
    for s in solutions:
      expt.goniometer.set_angles((0, s[1], s[2]))
      obs = self.predict_to_miller_set_with_shadow(expt, resolution)
      new = missing.common_set(obs)

      print '%8.3f %8.3f %d' % (s[1] * 180. / math.pi, s[2] * 180. / math.pi, \
                                len(new.indices()))

  def make_scan_360(self, scan):
    epochs = scan.get_epochs()
    exposure_times = scan.get_exposure_times()
    image_range = scan.get_image_range()
    oscillation = scan.get_oscillation()

    current = 1 + image_range[1] - image_range[0]
    turn = int(round(360./oscillation[1]))
    extra = turn - current

    for j in range(extra):
      epochs.append(0.0)
      exposure_times.append(0.0)
    image_range = image_range[0], image_range[1]+extra

    scan.set_image_range(image_range)
    scan.set_epochs(epochs)
    scan.set_exposure_times(exposure_times)
    return scan

  def predict_to_miller_set(self, expt, resolution):
    from dials.array_family import flex

    predicted = flex.reflection_table.from_predictions(expt, dmin=resolution)
    hkl = predicted['miller_index']

    # now get a full set of all unique miller indices
    from cctbx import miller
    from cctbx import crystal

    obs = miller.set(crystal_symmetry=crystal.symmetry(
      space_group=expt.crystal.get_space_group(),
      unit_cell=expt.crystal.get_unit_cell()), anomalous_flag=True,
      indices=hkl).unique_under_symmetry()

    return obs

  def predict_to_miller_set_with_shadow(self, expt, resolution):
    from dials.array_family import flex
    from dials.util import is_inside_polygon
    from dials.command_line.check_strategy import filter_shadowed_reflections
    masker = expt.imageset.reader().get_format().get_goniometer_shadow_masker()
    predicted = flex.reflection_table.from_predictions(expt, dmin=resolution)

    # transmogrify this to an ExperimentList from an Experiment
    from dxtbx.model import ExperimentList
    experiments = ExperimentList()
    experiments.append(expt)
    predicted['id'] = flex.int(predicted.size(), 0)

    shadowed = filter_shadowed_reflections(experiments, predicted)
    predicted = predicted.select(~shadowed)

    hkl = predicted['miller_index']

    # now get a full set of all unique miller indices
    from cctbx import miller
    from cctbx import crystal

    obs = miller.set(crystal_symmetry=crystal.symmetry(
      space_group=expt.crystal.get_space_group(),
      unit_cell=expt.crystal.get_unit_cell()), anomalous_flag=True,
      indices=hkl).unique_under_symmetry()

    return obs


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
