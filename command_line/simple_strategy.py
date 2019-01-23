# LIBTBX_SET_DISPATCHER_NAME dev.dials.simple_strategy

from __future__ import absolute_import, division, print_function

import copy
import math

import iotbx.phil
from dials.array_family import flex
from scitbx import matrix

help_message = '''

'''

phil_scope= iotbx.phil.parse("""

include scope dials.algorithms.profile_model.factory.phil_scope

d_min = None
  .type = float
  .help = "Minimum d-spacing of predicted reflections"

unit_cell_scale = 2
  .type = float

space_group = None
  .type = space_group

degrees_per_bin = 5
  .type = float(value_min=0)

minimum_fraction_new = 0.001
  .type = float(value_min=0, value_max=1)

""")


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  import libtbx.load_env

  usage = "%s [options] datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)

  if len(experiments) == 0:
    parser.print_help()
    exit(0)

  imagesets = experiments.imagesets()

  predicted_all = None
  dose = flex.size_t()

  for i_expt, expt in enumerate(experiments):
    if params.space_group is not None:
      expt.crystal.set_space_group(params.space_group.group())
    strategy = Strategy(expt, d_min=params.d_min,
                        unit_cell_scale=params.unit_cell_scale,
                        degrees_per_bin=params.degrees_per_bin,
                        min_frac_new=params.minimum_fraction_new)
    strategy.plot(prefix='strategy1_')

    expt2 = copy.deepcopy(expt)
    scan = expt2.scan
    gonio = expt2.goniometer
    angles = gonio.get_angles()
    theta_max = strategy.theta_max
    fixed_rotation = matrix.sqr(gonio.get_fixed_rotation())
    setting_rotation = matrix.sqr(gonio.get_setting_rotation())
    rotation_axis = matrix.col(gonio.get_rotation_axis_datum())
    rotation_matrix = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
      scan.get_oscillation()[0], deg=True)
    D_p = (setting_rotation * rotation_matrix * fixed_rotation)

    beam = expt2.beam
    s0 = matrix.col(beam.get_unit_s0())

    # rotate crystal by at least 2 * theta_max around axis perpendicular to
    # goniometer rotation axis
    n = rotation_axis.cross(s0)

    solutions = flex.vec3_double()
    rotation_angles = flex.double()
    for sign in (1, -1):
      i = 0
      while True:
        rot_angle = 2 * sign * (theta_max + i)
        i += 1

        R = n.axis_and_angle_as_r3_rotation_matrix(rot_angle, deg=True)

        axes = gonio.get_axes()
        assert len(axes) == 3
        e1, e2, e3 = (matrix.col(e) for e in reversed(axes))

        from dials.algorithms.refinement import rotation_decomposition
        solns = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
          R * D_p, e1, e2, e3, return_both_solutions=True, deg=True)
        if solns is not None:
          solutions.extend(flex.vec3_double(solns))
          for i in range(len(solns)):
            rotation_angles.append(rot_angle)
          break

    for rot_angle, solution in zip(rotation_angles, solutions):
      angles = reversed(solution)
      gonio.set_angles(angles)

      print()
      print("Goniometer settings to rotate crystal by %.2f degrees:" %rot_angle)
      for name, angle in zip(gonio.get_names(), gonio.get_angles()):
        print("%s: %.2f degrees" %(name, angle))
      print()

    strategy2 = Strategy(expt2, d_min=params.d_min,
                         unit_cell_scale=params.unit_cell_scale,
                         degrees_per_bin=params.degrees_per_bin,
                         min_frac_new=params.minimum_fraction_new)
    strategy2.plot(prefix='strategy2_')

    stats = ComputeStats([strategy, strategy2])
    stats.show()
    plot_statistics(stats, prefix='multi_strategy_')


class Strategy(object):

  def __init__(self, experiment, other=None, d_min=None, unit_cell_scale=1, degrees_per_bin=5,
               min_frac_new=0.001):
    print(experiment.goniometer)
    print(experiment.scan)
    self.experiment = copy.deepcopy(experiment)
    self.other = other
    self.unit_cell_scale = unit_cell_scale
    self.degrees_per_bin = degrees_per_bin
    self.min_frac_new = min_frac_new

    scan = self.experiment.scan
    detector = self.experiment.detector
    beam = self.experiment.beam
    crystal = self.experiment.crystal

    from cctbx import uctbx
    s = self.unit_cell_scale
    assert s > 0
    uc = crystal.get_unit_cell()
    a, b, c, alpha, beta, gamma = uc.parameters()
    uc_mod = uctbx.unit_cell((a/s, b/s, c/s, alpha, beta, gamma))
    crystal.set_B(matrix.sqr(uc_mod.fractionalization_matrix()).transpose())

    self.d_min = d_min
    if self.d_min is None:
      self.d_min = detector.get_max_inscribed_resolution(beam.get_s0())

    # bragg's law
    sin_theta = 0.5 * beam.get_wavelength()/self.d_min
    theta_max_rad = math.asin(sin_theta)
    self.theta_max = theta_max_rad * 180/math.pi

    print("theta_max (degrees): %.2f" %self.theta_max)

    Btot = 1 - 3 * (4 * theta_max_rad - math.sin(4 * theta_max_rad))/(32 * (math.sin(theta_max_rad)**3))
    print(Btot)

    # Section 2.9, Dauter Acta Cryst. (1999). D55, 1703-1717
    #max_rotation = 360 + 2 * self.theta_max
    max_rotation = 360

    image_range = scan.get_image_range()
    oscillation = scan.get_oscillation()
    scan.set_image_range((image_range[0],
                          image_range[0]+int(max_rotation/oscillation[1])))
    self.predict()
    x, y, z = self.predicted['xyzcal.px'].parts()
    self.predicted['dose'] = flex.size_t(
      list(flex.ceil(z * oscillation[1]/self.degrees_per_bin).iround()))
    self.stats = ComputeStats([self], degrees_per_bin=self.degrees_per_bin)
    self.ieither_completeness = self.stats.ieither_completeness
    self.iboth_completeness = self.stats.iboth_completeness
    self.frac_new_ref = self.stats.frac_new_ref
    self.frac_new_pairs = self.stats.frac_new_pairs
    self.determine_cutoffs(self.min_frac_new)
    self.show()

  def predict(self):
    # Populate the reflection table with predictions
    self.predicted = flex.reflection_table.from_predictions(
      self.experiment,
      force_static=True,
      dmin=self.d_min)
    self.predicted['id'] = flex.int(len(self.predicted), 0)

  def determine_fraction_new_cutoff(self, fraction_new, cutoff):
    imax = flex.max_index(fraction_new)
    isel = (fraction_new < cutoff).iselection()
    return isel[(isel > imax).iselection()[0]]

  def determine_cutoffs(self, min_frac_new=None):
    if min_frac_new is None:
      self.min_frac_new = min_frac_new
    self.cutoff_non_anom = self.degrees_per_bin * self.determine_fraction_new_cutoff(
      self.frac_new_ref, min_frac_new)
    self.cutoff_anom = self.degrees_per_bin * self.determine_fraction_new_cutoff(
      self.frac_new_pairs, min_frac_new)

  def show(self):
    print("Warning: Use of the .show() method is deprecated. Use print(object) instead.")
    print(str(self))

  def __str__(self):
    output = [ str(self.stats) ]
    output.append("Suggested cutoff (non-anom): %.2f degrees" % self.cutoff_non_anom)
    output.append("  (completeness: %.2f %%)" % (
      100 * self.ieither_completeness[int(self.cutoff_non_anom/self.degrees_per_bin)]))
    output.append("Suggested cutoff (anom): %.2f degrees" % self.cutoff_anom)
    output.append("  (completeness: %.2f %%)" % (
      100 * self.iboth_completeness[int(self.cutoff_anom/self.degrees_per_bin)]))
    return "\n".join(output)

  def plot(self, prefix=''):
    plot_statistics(self.stats, prefix=prefix, degrees_per_bin=self.degrees_per_bin,
                    cutoff_anom=self.cutoff_anom,
                    cutoff_non_anom=self.cutoff_non_anom)

class ComputeStats(object):

  def __init__(self, strategies, n_bins=8, degrees_per_bin=5):
    from cctbx import crystal, miller

    sg = strategies[0].experiment.crystal.get_space_group() \
      .build_derived_reflection_intensity_group(anomalous_flag=True)
    cs = crystal.symmetry(
      unit_cell=strategies[0].experiment.crystal.get_unit_cell(), space_group=sg)

    for i, strategy in enumerate(strategies):
      if i == 0:
        predicted = copy.deepcopy(strategy.predicted)
      else:
        predicted_ = copy.deepcopy(strategy.predicted)
        predicted_['dose'] += (flex.max(predicted['dose']) + 1)
        predicted.extend(predicted_)
    ms = miller.set(cs, indices=predicted['miller_index'], anomalous_flag=True)
    ma = miller.array(ms, data=flex.double(ms.size(),1),
                      sigmas=flex.double(ms.size(), 1))
    if 1:
      merging = ma.merge_equivalents()
      o = merging.array().customized_copy(
        data=merging.redundancies().data().as_double()).as_mtz_dataset('I').mtz_object()
      o.write('predicted.mtz')

    d_star_sq = ma.d_star_sq().data()

    binner = ma.setup_binner_d_star_sq_step(
      d_star_sq_step=(flex.max(d_star_sq)-flex.min(d_star_sq)+1e-8)/n_bins)

    dose = predicted['dose']
    range_width = 1
    range_min = flex.min(dose) - range_width
    range_max = flex.max(dose)
    n_steps = 2 + int((range_max - range_min) - range_width)

    binner_non_anom = ma.as_non_anomalous_array().use_binning(
      binner)
    self.n_complete = flex.size_t(binner_non_anom.counts_complete()[1:-1])

    from xia2.Modules.PyChef2 import ChefStatistics
    chef_stats = ChefStatistics(
      ma.indices(), ma.data(), ma.sigmas(),
      ma.d_star_sq().data(), dose, self.n_complete, binner,
      ma.space_group(), ma.anomalous_flag(), n_steps)

    def fraction_new(completeness):
      # Completeness so far at end of image
      completeness_end = completeness[1:]
      # Completeness so far at start of image
      completeness_start = completeness[:-1]
      # Fraction of unique reflections observed for the first time on each image
      return completeness_end - completeness_start

    self.dose = dose
    self.ieither_completeness = chef_stats.ieither_completeness()
    self.iboth_completeness = chef_stats.iboth_completeness()
    self.frac_new_ref = fraction_new(self.ieither_completeness) / degrees_per_bin
    self.frac_new_pairs = fraction_new(self.iboth_completeness) / degrees_per_bin

  def show(self):
    print("Warning: Use of the .show() method is deprecated. Use print(object) instead.")
    print(str(self))

  def __str__(self):
    return "\n".join((
        "Max. completeness (non-anom): %.2f %%" % (100 * flex.max(self.ieither_completeness)),
        "Max. completeness (anom): %.2f %%" %(100 * flex.max(self.iboth_completeness)),
    ))

def plot_statistics(statistics, prefix='', degrees_per_bin=5,
                    cutoff_anom=None, cutoff_non_anom=None):

  range_width = 1
  range_min = flex.min(statistics.dose) - range_width
  range_max = flex.max(statistics.dose)
  n_steps = 2 + int((range_max - range_min) - range_width)
  x = flex.double_range(n_steps) * range_width + range_min
  x *= degrees_per_bin

  dpi = 300
  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot as plt
  try: plt.style.use('ggplot')
  except AttributeError: pass
  line1,  = plt.plot(x, statistics.ieither_completeness, label='Unique reflections')
  line2,  = plt.plot(x, statistics.iboth_completeness, label='Bijvoet pairs')
  if cutoff_non_anom is not None:
    plt.plot([cutoff_non_anom, cutoff_non_anom], plt.ylim(), c=line1.get_color(), linestyle='dashed')
  if cutoff_anom is not None:
    plt.plot([cutoff_anom, cutoff_anom], plt.ylim(), c=line2.get_color(), linestyle='dotted')
  plt.xlim(0, plt.xlim()[1])
  plt.xlabel('Scan angle (degrees)')
  plt.ylabel('Completeness (%)')
  plt.ylim(0, 1)
  plt.legend(loc='lower right', fontsize='small')
  plt.savefig('%scompleteness_vs_scan_angle.png' %prefix, dpi=dpi)
  plt.clf()

  line1, = plt.plot(x[1:], 100 * statistics.frac_new_ref, label='Unique reflections')
  line2, = plt.plot(x[1:], 100 * statistics.frac_new_pairs, label='Bijvoet pairs')
  ylim = plt.ylim()
  if cutoff_non_anom is not None:
    plt.plot([cutoff_non_anom, cutoff_non_anom], ylim, c=line1.get_color(), linestyle='dashed')
  if cutoff_anom is not None:
    plt.plot([cutoff_anom, cutoff_anom], ylim, c=line2.get_color(), linestyle='dotted')
  plt.ylim(ylim)
  plt.xlim(0, plt.xlim()[1])
  plt.xlabel('Scan angle (degrees)')
  plt.ylabel('% new reflections per degree')
  plt.legend(loc='upper right', fontsize='small')
  plt.savefig('%spercent_new_reflections_vs_scan_angle.png' %prefix, dpi=dpi)
  plt.clf()



if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
