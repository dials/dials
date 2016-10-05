# LIBTBX_SET_DISPATCHER_NAME dev.dials.simple_strategy

from __future__ import division
from dials.array_family import flex
import iotbx.phil

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
    strategy.plot()
    return


class Strategy(object):

  def __init__(self, experiment, d_min=None, unit_cell_scale=1, degrees_per_bin=5,
               min_frac_new=0.001):
    self.experiment = experiment
    self.unit_cell_scale = unit_cell_scale
    self.degrees_per_bin = degrees_per_bin
    self.min_frac_new = min_frac_new

    scan = self.experiment.scan
    detector = self.experiment.detector
    beam = self.experiment.beam
    crystal = self.experiment.crystal

    from cctbx import uctbx
    from scitbx import matrix
    s = self.unit_cell_scale
    assert s > 0
    uc = crystal.get_unit_cell()
    a, b, c, alpha, beta, gamma = uc.parameters()
    uc_mod = uctbx.unit_cell((a/s, b/s, c/s, alpha, beta, gamma))
    crystal.set_B(matrix.sqr(uc_mod.fractionalization_matrix()).transpose())

    self.d_min = d_min
    if self.d_min is None:
      self.d_min = detector.get_max_inscribed_resolution(beam.get_s0())

    import math
    # bragg's law
    theta_max = math.asin(0.5 * beam.get_wavelength()/self.d_min) * 180/math.pi

    # Section 2.9, Dauter Acta Cryst. (1999). D55, 1703-1717
    #max_rotation = 360 + 2 * theta_max
    max_rotation = 360

    image_range = scan.get_image_range()
    oscillation = scan.get_oscillation()
    scan.set_image_range((image_range[0],
                          image_range[0]+int(max_rotation/oscillation[1])))
    self.predict()
    x, y, z = self.predicted['xyzcal.px'].parts()
    self.predicted['dose'] = flex.size_t(
      list(flex.ceil(z * oscillation[1]/self.degrees_per_bin).iround()))
    self.compute_stats()
    self.determine_cutoffs(self.min_frac_new)
    self.show()

  def predict(self):
    # Populate the reflection table with predictions
    self.predicted = flex.reflection_table.from_predictions(
      self.experiment,
      force_static=True,
      dmin=self.d_min)
    self.predicted['id'] = flex.int(len(self.predicted), 0)

  def compute_stats(self, n_bins=8):
    from cctbx import crystal, miller

    cs = crystal.symmetry(unit_cell=self.experiment.crystal.get_unit_cell(),
                          space_group=self.experiment.crystal.get_space_group())
    ms = miller.set(cs, indices=self.predicted['miller_index'], anomalous_flag=True)
    ma = miller.array(ms, data=flex.double(ms.size(), 1),
                      sigmas=flex.double(ms.size(), 1))
    if 1:
      o = ma.merge_equivalents().array().as_mtz_dataset('I').mtz_object()
      o.write('predicted.mtz')

    d_star_sq = ma.d_star_sq().data()

    binner = ma.setup_binner_d_star_sq_step(
      d_star_sq_step=(flex.max(d_star_sq)-flex.min(d_star_sq)+1e-8)/n_bins)

    dose = self.predicted['dose']
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

    self.ieither_completeness = chef_stats.ieither_completeness()
    self.iboth_completeness = chef_stats.iboth_completeness()
    self.frac_new_ref = fraction_new(self.ieither_completeness) / self.degrees_per_bin
    self.frac_new_pairs = fraction_new(self.iboth_completeness) / self.degrees_per_bin

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
    print "Max. completeness (non-anom): %.2f %%" %(100 * flex.max(self.ieither_completeness))
    print "Max. completeness (anom): %.2f %%" %(100 * flex.max(self.iboth_completeness))
    print "Suggested cutoff (non-anom): %.2f degrees" %self.cutoff_non_anom
    print "  (completeness: %.2f %%)" %(
      100 * self.ieither_completeness[int(self.cutoff_non_anom/self.degrees_per_bin)])
    print "Suggested cutoff (anom): %.2f degrees" %self.cutoff_anom
    print "  (completeness: %.2f %%)" %(
      100 * self.iboth_completeness[int(self.cutoff_anom/self.degrees_per_bin)])

  def plot(self, prefix=''):

    dose = self.predicted['dose']
    range_width = 1
    range_min = flex.min(dose) - range_width
    range_max = flex.max(dose)
    n_steps = 2 + int((range_max - range_min) - range_width)
    x = flex.double_range(n_steps) * range_width + range_min
    x *= self.degrees_per_bin

    dpi = 300
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    plt.style.use('ggplot')
    #plt.style.use('seaborn-colorblind')
    #plt.style.use('fivethirtyeight')
    line,  = plt.plot(x, self.ieither_completeness, label='Unique reflections')
    plt.plot([self.cutoff_non_anom, self.cutoff_non_anom], plt.ylim(), c=line.get_color(), linestyle='dashed')
    line,  = plt.plot(x, self.iboth_completeness, label='Bijvoet pairs')
    plt.plot([self.cutoff_anom, self.cutoff_anom], plt.ylim(), c=line.get_color(), linestyle='dotted')
    plt.xlim(0, plt.xlim()[1])
    plt.xlabel('Scan angle (degrees)')
    plt.ylabel('Completeness (%)')
    plt.ylim(0, 1)
    plt.legend(loc='lower right', fontsize='small')
    plt.savefig('%scompleteness_vs_scan_angle.png' %prefix, dpi=dpi)
    plt.clf()

    scale = self.unit_cell_scale**3
    line1, = plt.plot(
      x[1:], scale * self.frac_new_ref*flex.sum(self.n_complete),
      label='Unique reflections')
    line2, = plt.plot(
      x[1:], scale * self.frac_new_pairs*flex.sum(self.n_complete),
      label='Bijvoet pairs')
    plt.plot([self.cutoff_non_anom, self.cutoff_non_anom], plt.ylim(), c=line1.get_color(), linestyle='dashed')
    plt.plot([self.cutoff_anom, self.cutoff_anom], plt.ylim(), c=line2.get_color(), linestyle='dotted')
    plt.xlim(0, plt.xlim()[1])
    plt.xlabel('Scan angle (degrees)')
    plt.ylabel('# new reflections')
    plt.legend(loc='upper right', fontsize='small')
    plt.savefig('%sn_new_reflections_vs_scan_angle.png' %prefix, dpi=dpi)
    plt.clf()

    line1, = plt.plot(x[1:], 100 * self.frac_new_ref, label='Unique reflections')
    line2, = plt.plot(x[1:], 100 * self.frac_new_pairs, label='Bijvoet pairs')
    ylim = plt.ylim()
    plt.plot([self.cutoff_non_anom, self.cutoff_non_anom], ylim, c=line1.get_color(), linestyle='dashed')
    plt.plot([self.cutoff_anom, self.cutoff_anom], ylim, c=line2.get_color(), linestyle='dotted')
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
