# LIBTBX_SET_DISPATCHER_NAME dev.dials.estimate_resolution_limit

from __future__ import absolute_import, division, print_function

from libtbx.phil import command_line
from dials.util import Sorry
import iotbx.phil
# from dials.util.command_line import Importer
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_experiments
from dials.array_family import flex

help_message = '''

'''

phil_scope = iotbx.phil.parse("""
i_sigi_cutoff = 0.1
  .type = float(value_min=0.000001)
plot=False
  .type = bool
""", process_includes=True)


def run(args):
  import libtbx.load_env
  usage = "%s experiments.json indexed.pickle [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  if len(experiments) == 0:
    parser.print_help()
    return
  elif len(experiments) > 1:
    raise Sorry("More than one experiment present")

  experiment = experiments[0]
  assert(len(reflections) == 1)
  reflections = reflections[0]

  intensities = reflections['intensity.sum.value']
  variances = reflections['intensity.sum.variance']
  if 'intensity.prf.value' in reflections:
    intensities = reflections['intensity.prf.value']
    variances = reflections['intensity.prf.variance']
  sel = (variances > 0)
  intensities = intensities.select(sel)
  variances = variances.select(sel)
  sigmas = flex.sqrt(variances)
  indices = reflections['miller_index'].select(sel)

  from cctbx import crystal, miller
  crystal_symmetry = crystal.symmetry(
    space_group=experiment.crystal.get_space_group(),
    unit_cell=experiment.crystal.get_unit_cell())

  miller_set = miller.set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=True,
    indices=indices)
  miller_array = miller.array(
    miller_set=miller_set,
    data=intensities,
    sigmas=sigmas).set_observation_type_xray_intensity()

  #miller_array.setup_binner(n_bins=50, reflections_per_bin=100)
  miller_array.setup_binner(auto_binning=True, n_bins=20)
  result = miller_array.i_over_sig_i(use_binning=True)
  result.show()

  from cctbx import uctbx
  d_star_sq_centre = result.binner.bin_centers(2)
  i_over_sig_i = flex.double(
    [d if d is not None else 0 for d in result.data[1:-1]])
  sel = (i_over_sig_i > 0)
  d_star_sq_centre = d_star_sq_centre.select(sel)
  i_over_sig_i = i_over_sig_i.select(sel)
  log_i_over_sig_i = flex.log(i_over_sig_i)
  weights = result.binner.counts()[1:-1].as_double().select(sel)
  fit = flex.linear_regression(
    d_star_sq_centre, log_i_over_sig_i, weights=weights)

  m = fit.slope()
  c = fit.y_intercept()

  import math
  y_cutoff = math.log(params.i_sigi_cutoff)
  x_cutoff = (y_cutoff - c)/m

  estimated_d_min = uctbx.d_star_sq_as_d(x_cutoff)
  print("estimated d_min: %.2f" %estimated_d_min)

  if params.plot:
    from matplotlib import pyplot
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)

    ax.plot(
      list(d_star_sq_centre),
      list(log_i_over_sig_i),
      label=r"ln(I/sigI)")
    ax.plot(pyplot.xlim(), [(m * x + c) for x in pyplot.xlim()], color='red')
    ax.plot([x_cutoff, x_cutoff], pyplot.ylim(), color='grey', linestyle='dashed')
    ax.plot(pyplot.xlim(), [y_cutoff, y_cutoff], color='grey', linestyle='dashed')
    ax.set_xlabel("d_star_sq")
    ax.set_ylabel("ln(I/sigI)")

    ax_ = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
    xticks = ax.get_xticks()
    xlim = ax.get_xlim()
    xticks_d = [
      uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks ]
    xticks_ = [ds2/(xlim[1]-xlim[0]) for ds2 in xticks]
    ax_.set_xticks(xticks)
    ax_.set_xlim(ax.get_xlim())
    ax_.set_xlabel(r"Resolution ($\AA$)")
    ax_.set_xticklabels(["%.1f" %d for d in xticks_d])
    pyplot.savefig("estimate_resolution_limit.png")
    pyplot.clf()


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
