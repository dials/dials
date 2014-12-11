from __future__ import division

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

def estimate_resolution_limit(reflections, imageset):
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

  if m_upper == m_lower:
    intersection = (-1,-1)
    resolution_estimate = -1

  else:
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

  return resolution_estimate


def resolution_histogram(reflections, imageset):
  reflections = map_to_reciprocal_space(reflections, imageset)
  d_star_sq = flex.pow2(reflections['rlp'].norms())
  hist = get_histogram(d_star_sq)

  # FIXME analyse here

if __name__ == '__main__':
  import sys

  # filter command parameters from filenames
  cl = []
  filenames = []
  for arg in sys.argv[1:]:
    if '=' in arg:
      cl.append(arg)
    else:
      filenames.append(arg)

  from dials.command_line.find_spots import phil_scope as params
  from dxtbx.datablock import DataBlockFactory
  from dials.array_family import flex
  interp = params.command_line_argument_interpreter()
  for cla in cl:
    params = params.fetch(interp.process(cla))
  params_extract = params.extract()
  for filename in filenames:
    datablock = DataBlockFactory.from_filenames([filename])[0]
    reflections = flex.reflection_table.from_observations(
      datablock, params_extract)

    imageset = datablock.extract_imagesets()[0]

    resolution_histogram(reflections, imageset)
    estimated_d_min = estimate_resolution_limit(reflections, imageset)
    print "%s: %.2f" % (filename, estimated_d_min)
