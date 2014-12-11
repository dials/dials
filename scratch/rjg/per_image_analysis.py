from __future__ import division

from libtbx import group_args
from dials.array_family import flex

def map_to_reciprocal_space(reflections, imageset):
  detector = imageset.get_detector()
  scan = imageset.get_scan()
  goniometer = imageset.get_goniometer()
  beam = imageset.get_beam()

  from dials.algorithms.indexing import indexer
  reflections = indexer.indexer_base.map_spots_pixel_to_mm_rad(
    reflections, detector, scan)
  indexer.indexer_base.map_centroids_to_reciprocal_space(
  reflections, detector, beam, goniometer)

  return reflections

def get_histogram(d_star_sq, target_n_per_bin=10, max_slots=20, min_slots=5):
  n_slots = len(d_star_sq)//target_n_per_bin
  n_slots = min(n_slots, max_slots)
  n_slots = max(n_slots, min_slots)
  return flex.histogram(d_star_sq, n_slots=n_slots)


def estimate_resolution_limit(reflections, imageset, plot_filename=None):
  reflections = map_to_reciprocal_space(reflections, imageset)

  d_star_sq = flex.pow2(reflections['rlp'].norms())

  from cctbx import uctbx
  d_spacings = 1/reflections['rlp'].norms()
  assert d_star_sq == uctbx.d_as_d_star_sq(d_spacings)

  intensities = reflections['intensity.sum.value']
  variances = reflections['intensity.sum.variance']

  sel = variances > 0
  intensities = intensities.select(sel)
  variances = intensities.select(sel)

  i_over_sigi = intensities/flex.sqrt(variances)
  log_i_over_sigi = flex.log(i_over_sigi)

  fit = flex.linear_regression(d_star_sq, log_i_over_sigi)
  m = fit.slope()
  c = fit.y_intercept()

  import math

  log_i_sigi_lower = flex.double()
  d_star_sq_lower = flex.double()
  log_i_sigi_upper = flex.double()
  d_star_sq_upper = flex.double()
  weights = flex.double()

  hist = get_histogram(d_star_sq)

  i_slot_max = flex.max_index(hist.slots())

  for i_slot, slot in enumerate(hist.slot_infos()):
    sel = (d_star_sq > slot.low_cutoff) & (d_star_sq < slot.high_cutoff)
    if sel.count(True) == 0:
      if i_slot > i_slot_max:
        hist = get_histogram(d_star_sq.select(d_star_sq < slot.low_cutoff))
        break

  low_percentile_limit = 0.05
  upper_percentile_limit = 1-low_percentile_limit
  for i_slot, slot in enumerate(hist.slot_infos()):
    sel = (d_star_sq > slot.low_cutoff) & (d_star_sq < slot.high_cutoff)
    if sel.count(True) == 0:
      if i_slot > i_slot_max:
        break
      else:
        continue
    log_i_over_sigi_sel = log_i_over_sigi.select(sel)
    d_star_sq_sel = d_star_sq.select(sel)
    perm = flex.sort_permutation(log_i_over_sigi_sel)
    i_lower = perm[int(math.floor(low_percentile_limit * len(perm)))]
    i_upper = perm[int(math.floor(upper_percentile_limit * len(perm)))]
    log_i_sigi_lower.append(log_i_over_sigi_sel[i_lower])
    log_i_sigi_upper.append(log_i_over_sigi_sel[i_upper])
    d_star_sq_upper.append(d_star_sq_sel[i_lower])
    d_star_sq_lower.append(d_star_sq_sel[i_upper])
    weights.append(slot.n)

  fit_upper = flex.linear_regression(
    d_star_sq_upper, log_i_sigi_upper, weights=weights)
  m_upper = fit_upper.slope()
  c_upper = fit_upper.y_intercept()
  fit_lower = flex.linear_regression(
    d_star_sq_lower, log_i_sigi_lower, weights=weights)
  m_lower = fit_lower.slope()
  c_lower = fit_lower.y_intercept()

  #fit_upper.show_summary()
  #fit_lower.show_summary()

  if m_upper == m_lower:
    intersection = (-1,-1)
    resolution_estimate = -1

  else:
    # http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_the_equations_of_the_lines
    intersection = (
      (c_lower-c_upper)/(m_upper-m_lower),
      (m_upper*c_lower-m_lower*c_upper)/(m_upper-m_lower))

    a = m_upper
    c_ = c_upper
    b = m_lower
    d = c_lower
    assert intersection == ((d-c_)/(a-b), (a*d-b*c_)/(a-b))

    d_star_sq_estimate = intersection[0]
    resolution_estimate = uctbx.d_star_sq_as_d(d_star_sq_estimate)

  resolution_estimate = max(resolution_estimate, flex.min(d_spacings))

  if plot_filename is not None:
    from cctbx import uctbx
    from matplotlib import pyplot
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(d_star_sq, log_i_over_sigi, marker='+')
    ax.scatter(d_star_sq_upper, log_i_sigi_upper, marker='+', color='red')
    ax.scatter(d_star_sq_lower, log_i_sigi_lower, marker='+', color='red')
    ax.scatter([intersection[0]], [intersection[1]], marker='x', s=50, color='b')
    #ax.hexbin(d_star_sq, log_i_over_sigi, gridsize=30)
    ax.plot(pyplot.xlim(), [(m * x + c) for x in pyplot.xlim()])
    ax.plot(pyplot.xlim(), [(m_upper * x + c_upper) for x in pyplot.xlim()], color='red')
    ax.plot(pyplot.xlim(), [(m_lower * x + c_lower) for x in pyplot.xlim()], color='red')
    ax.set_xlabel('d_star_sq')
    ax.set_ylabel('ln(I/sigI)')
    ax.set_xlim((max(-ax.get_xlim()[1], -0.05), ax.get_xlim()[1]))
    ax.set_ylim((0, ax.get_ylim()[1]))
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
    #pyplot.show()
    pyplot.savefig(plot_filename)
    pyplot.clf()

  return resolution_estimate


def resolution_histogram(reflections, imageset, plot_filename=None):
  reflections = map_to_reciprocal_space(reflections, imageset)
  d_star_sq = flex.pow2(reflections['rlp'].norms())
  hist = get_histogram(d_star_sq)

  if plot_filename is not None:
    from matplotlib import pyplot
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.bar(hist.slot_centers()-0.5*hist.slot_width(), hist.slots(),
               width=hist.slot_width())
    ax.set_xlabel("d_star_sq")
    ax.set_ylabel("Frequency")

    ax_ = ax.twiny() # ax2 is responsible for "top" axis and "right" axis
    xticks = ax.get_xticks()
    xlim = ax.get_xlim()
    from cctbx import uctbx
    xticks_d = [
      uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks ]
    xticks_ = [ds2/(xlim[1]-xlim[0]) for ds2 in xticks]
    ax_.set_xticks(xticks)
    ax_.set_xlim(ax.get_xlim())
    ax_.set_xlabel(r"Resolution ($\AA$)")
    ax_.set_xticklabels(["%.1f" %d for d in xticks_d])
    #pyplot.show()
    pyplot.savefig(plot_filename)
    pyplot.clf()

def plot_ordered_d_star_sq(reflections, imageset):
  reflections = map_to_reciprocal_space(reflections, imageset)
  d_star_sq = flex.pow2(reflections['rlp'].norms())

  from matplotlib import pyplot
  perm = flex.sort_permutation(d_star_sq)
  pyplot.scatter(list(range(len(perm))), list(d_star_sq.select(perm)), marker='+')
  pyplot.show()


def stats_single_image(imageset, reflections, i=None):
  if i is not None:
    filename = "i_over_sigi_vs_resolution_%d.png" %i
  #plot_ordered_d_star_sq(reflections, imageset)
  n_spots_total = len(reflections)
  estimated_d_min = estimate_resolution_limit(
    reflections, imageset, plot_filename=filename)

  return group_args(n_spots_total=n_spots_total,
                    estimated_d_min=estimated_d_min)

def stats_imageset(imageset, reflections):
  n_spots_total = []
  estimated_d_min = []

  image_number = reflections['xyzobs.px.value'].parts()[2]
  image_number = flex.floor(image_number)

  for i in range(len(imageset)):
    stats = stats_single_image(
      imageset[i:i+1], reflections.select(image_number==i), i=i)
    n_spots_total.append(stats.n_spots_total)
    estimated_d_min.append(stats.estimated_d_min)

  return group_args(n_spots_total=n_spots_total,
                    estimated_d_min=estimated_d_min)


def table(stats):
  n_spots_total = stats.n_spots_total
  estimated_d_min = stats.estimated_d_min
  rows = [("image", "#spots")]
  for i_image in range(len(n_spots_total)):
    rows.append((str(int(i_image)),
                 str(n_spots_total[i_image]),
                 "%.2f" %estimated_d_min[i_image]))
  return rows

def print_table(stats, out=None):
  if out is None:
    import sys
    out = sys.stdout
  from libtbx import table_utils

  rows = table(stats)
  print >> out, table_utils.format(
    rows, has_header=True, prefix="|", postfix="|")


if __name__ == '__main__':
  from dials.util.options import OptionParser
  from dials.util.options import flatten_reflections, flatten_datablocks

  parser = OptionParser(
    read_reflections=True,
    read_datablocks=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=False)
  reflections = flatten_reflections(params.input.reflections)
  datablocks = flatten_datablocks(params.input.datablock)

  assert len(reflections) == 1
  assert len(datablocks) == 1

  reflections = reflections[0]
  imageset = datablocks[0].extract_imagesets()[0]

  stats = stats_imageset(imageset, reflections)
  print_table(stats)

  #resolution_histogram(
  #  reflections, imageset, plot_filename="spot_count_vs_resolution.png")
  #estimated_d_min = estimate_resolution_limit(
  #  reflections, imageset, plot_filename="i_over_sigi_vs_resolution.png")
  #print "Estimated d_min: %.2f" %estimated_d_min
