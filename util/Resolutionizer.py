from __future__ import absolute_import, division, print_function

import copy
import math
import sys
import time

import iotbx.phil
from cctbx.array_family import flex
from libtbx.utils import Sorry
from scitbx import lbfgs

def nint(a):
  return int(round(a))

start_time = time.time()
def stamp(message):
  if False:
    print("[%7.3f] %s" % (time.time() - start_time, message))
  return

def poly_residual(xp, y, params):
  '''Compute the residual between the observations y[i] and sum_j
  params[j] x[i]^j. For efficiency, x[i]^j are pre-calculated in xp.'''

  r = 0.0

  n = len(params)
  c = len(y)

  e = flex.double([flex.sum(xp[j] * params) for j in range(c)])

  return flex.sum(flex.pow2(y - e))

def poly_gradients(xp, y, params):
  '''Compute the gradient of the residual w.r.t. the parameters, N.B.
  will be performed using a finite difference method. N.B. this should
  be trivial to do algebraicly.'''

  eps = 1.0e-6

  g = flex.double()

  n = len(params)

  for j in range(n):
    rs = []
    for signed_eps in [- eps, eps]:
      params_eps = params[:]
      params_eps[j] += signed_eps
      rs.append(poly_residual(xp, y, params_eps))
    g.append((rs[1] - rs[0]) / (2 * eps))

  return g

class poly_fitter(object):
  '''A class to do the polynomial fit. This will fit observations y
  at points x with a polynomial of order n.'''

  def __init__(self, points, values, order):
    self.x = flex.double([1.0 for j in range(order)])
    self._x = flex.double(points)
    self._y = flex.double(values)

    # precalculate x[j]^[0-(n - 1)] values

    self._xp = [flex.double([math.pow(x, j) for j in range(order)])
                for x in self._x]

    return

  def refine(self):
    '''Actually perform the parameter refinement.'''

    tp = lbfgs.termination_parameters(max_iterations=1000)
    r = lbfgs.run(target_evaluator = self, termination_params=tp)
    return r

  def compute_functional_and_gradients(self):

    return poly_residual(self._xp, self._y, self.x), \
           poly_gradients(self._xp, self._y, self.x)

  def get_parameters(self):
    return list(self.x)

  def evaluate(self, x):
    '''Evaluate the resulting fit at point x.'''

    return sum([math.pow(x, k) * self.x[k] for k in range(len(self.x))])

def fit(x, y, order):
  '''Fit the values y(x) then return this fit. x, y should
  be iterables containing floats of the same size. The order is the order
  of polynomial to use for this fit. This will be useful for e.g. I/sigma.'''

  stamp("fitter: %s %s %s" % (x, y, order))
  pf = poly_fitter(x, y, order)
  stamp("fitter: refine")
  pf.refine()
  stamp("fitter: done")

  return [pf.evaluate(_x) for _x in x]

def tanh_fit(x, y, iqr_multiplier=None):

  from scitbx.math import curve_fitting

  tf = curve_fitting.tanh_fit(x, y)
  f = curve_fitting.tanh(*tf.params)

  if iqr_multiplier is not None:
    assert iqr_multiplier > 0
    yc = f(x)
    dy = y - yc

    from scitbx.math import five_number_summary
    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(dy)
    iqr_x = q3_x - q1_x
    cut_x = iqr_multiplier * iqr_x
    outliers = (dy > q3_x + cut_x) | (dy < q1_x - cut_x)
    if outliers.count(True) > 0:
      xo = x.select(~outliers)
      yo = y.select(~outliers)
      tf = curve_fitting.tanh_fit(xo, yo)
      f = curve_fitting.tanh(*tf.params)

  return f(x)

def log_fit(x, y, order):
  '''Fit the values log(y(x)) then return exp() to this fit. x, y should
  be iterables containing floats of the same size. The order is the order
  of polynomial to use for this fit. This will be useful for e.g. I/sigma.'''

  ly = [math.log(_y) for _y in y]

  pf = poly_fitter(x, ly, order)
  pf.refine()

  return [math.exp(pf.evaluate(_x)) for _x in x]

def log_inv_fit(x, y, order):
  '''Fit the values log(1 / y(x)) then return the inverse of this fit.
  x, y should be iterables, the order of the polynomial for the transformed
  fit needs to be specified. This will be useful for e.g. Rmerge.'''

  ly = [math.log(1.0 / _y) for _y in y]

  pf = poly_fitter(x, ly, order)
  pf.refine()

  return [(1.0 / math.exp(pf.evaluate(_x))) for _x in x]

def interpolate_value(x, y, t):
  '''Find the value of x: y(x) = t.'''

  if t > max(y) or t < min(y):
    raise RuntimeError('t outside of [%f, %f]' % (min(y), max(y)))

  for j in range(1, len(x)):
    x0 = x[j - 1]
    y0 = y[j - 1]

    x1 = x[j]
    y1 = y[j]

    if (y0 - t) * (y1 - t) < 0:
      return x0 + (t - y0) * (x1 - x0) / (y1 - y0)

phil_str = '''
  rmerge = None
    .type = float(value_min=0)
    .help = "Maximum value of Rmerge in the outer resolution shell"
    .short_caption = "Outer shell Rmerge"
    .expert_level = 1
  completeness = None
    .type = float(value_min=0)
    .help = "Minimum completeness in the outer resolution shell"
    .short_caption = "Outer shell completeness"
    .expert_level = 1
  cc_ref = 0.1
    .type = float(value_min=0)
    .help = "Minimum value of CC vs reference dataset in the outer resolution shell"
    .short_caption = "Outer shell CCref"
    .expert_level = 1
  cc_half = 0.5
    .type = float(value_min=0)
    .help = "Minimum value of CC1/2 in the outer resolution shell"
    .short_caption = "Outer shell CC1/2"
    .expert_level = 1
  cc_half_method = *half_dataset sigma_tau
    .type = choice
  cc_half_significance_level = None
    .type = float(value_min=0, value_max=1)
    .expert_level = 1
  cc_half_fit = polynomial *tanh
    .type = choice
    .expert_level = 1
  isigma = 0.25
    .type = float(value_min=0)
    .help = "Minimum value of the unmerged <I/sigI> in the outer resolution shell"
    .short_caption = "Outer shell unmerged <I/sigI>"
    .expert_level = 1
  misigma = 1.0
    .type = float(value_min=0)
    .help = "Minimum value of the merged <I/sigI> in the outer resolution shell"
    .short_caption = "Outer shell merged <I/sigI>"
    .expert_level = 1
  i_mean_over_sigma_mean = None
    .type = float(value_min=0)
    .help = "Minimum value of the unmerged <I>/<sigI> in the outer resolution shell"
    .short_caption = "Outer shell unmerged <I>/<sigI>"
    .expert_level = 2
  nbins = 100
    .type = int
    .help = "Number of resolution bins to use for estimation of resolution limit."
    .short_caption = "Number of resolution bins."
    .expert_level = 1
  binning_method = *counting_sorted volume
    .type = choice
    .help = "Use equal-volume bins or bins with approximately equal numbers of reflections per bin."
    .short_caption = "Equal-volume or equal #ref binning."
    .expert_level = 1
  anomalous = False
    .type = bool
    .short_caption = "Keep anomalous pairs separate in merging statistics"
    .expert_level = 1
  labels = None
    .type = strings
  space_group = None
    .type = space_group
    .expert_level = 1
  reference = None
    .type = path
'''


phil_defaults = iotbx.phil.parse('''
resolutionizer {
%s
  batch_range = None
    .type = ints(size=2, value_min=0)
  plot = False
    .type = bool
    .expert_level = 2
}
''' %phil_str)


class resolution_plot(object):
  def __init__(self, ylabel):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
    self.ylabel = ylabel
    self.fig = pyplot.figure()
    self.ax = self.fig.add_subplot(111)

  def plot(self, d_star_sq, values, label):
    self.ax.plot(d_star_sq, values, label=label)
    if label.startswith('CC'):
      ylim = self.ax.get_ylim()
      self.ax.set_ylim(0, max(ylim[1], 1.05))

  def plot_resolution_limit(self, d):
    from cctbx import uctbx
    d_star_sq = uctbx.d_as_d_star_sq(d)
    self.ax.plot([d_star_sq, d_star_sq], self.ax.get_ylim(), linestyle='--')

  def savefig(self, filename):
    from cctbx import uctbx
    xticks = self.ax.get_xticks()
    xticks_d = [
      '%.2f' %uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks]
    self.ax.set_xticklabels(xticks_d)
    self.ax.set_xlabel('Resolution (A)')
    self.ax.set_ylabel(self.ylabel)
    self.ax.legend(loc='best')
    self.fig.savefig(filename)


class resolutionizer(object):
  '''A class to calculate things from merging reflections.'''

  def __init__(self, i_obs, batches, params, reference=None):

    self._params = params
    self._reference = reference

    if self._reference is not None:
      self._reference = self._reference.merge_equivalents(
        use_internal_variance=False).array()

    i_obs = i_obs.customized_copy(anomalous_flag=params.anomalous, info=i_obs.info())

    if self._params.batch_range is not None:
      batch_min, batch_max = self._params.batch_range
      assert batches is not None
      sel = (batches.data() >= batch_min) & (batches.data() <= batch_max)
      batches = batches.select(sel).set_info(batches.info())
      i_obs = i_obs.select(sel).set_info(i_obs.info())

    if self._params.space_group is not None:
      i_obs = i_obs.customized_copy(space_group_info=self._params.space_group,
                                    info=i_obs.info())

    self._intensities = i_obs

    import iotbx.merging_statistics
    self._merging_statistics = iotbx.merging_statistics.dataset_statistics(
      i_obs=i_obs,
      n_bins=self._params.nbins,
      cc_one_half_significance_level=self._params.cc_half_significance_level,
      cc_one_half_method=self._params.cc_half_method,
      binning_method=self._params.binning_method,
      anomalous=params.anomalous,
      use_internal_variance=False,
      eliminate_sys_absent=False,
      assert_is_not_unique_set_under_symmetry=False,
    )

  @classmethod
  def from_unmerged_mtz(cls, scaled_unmerged, params):

    def miller_array_from_mtz(unmerged_mtz):
      from iotbx import reflection_file_reader
      hkl_in = reflection_file_reader.any_reflection_file(scaled_unmerged)
      miller_arrays = hkl_in.as_miller_arrays(merge_equivalents=False)
      i_obs = None
      batches = None
      all_i_obs = []
      for array in miller_arrays :
        labels = array.info().label_string()
        if (array.is_xray_intensity_array()) :
          all_i_obs.append(array)
        if (labels == 'BATCH'):
          assert batches is None
          batches = array
      if (i_obs is None) :
        if (len(all_i_obs) == 0) :
          raise Sorry("No intensities found in %s." % file_name)
        elif (len(all_i_obs) > 1) :
          if params.labels is not None:
            from iotbx.reflection_file_utils import label_table
            lab_tab = label_table(all_i_obs)
            i_obs = lab_tab.select_array(label=params.labels[0], command_line_switch='labels')
          if i_obs is None:
            raise Sorry("Multiple intensity arrays - please specify one:\n%s" %
              "\n".join(["  labels=%s"%a.info().label_string() for a in all_i_obs]))
        else :
          i_obs = all_i_obs[0]
      if hkl_in.file_type() == 'ccp4_mtz':
        # need original miller indices otherwise we don't get correct anomalous
        # merging statistics
        mtz_object = hkl_in.file_content()
        if "M_ISYM" in mtz_object.column_labels():
          indices = mtz_object.extract_original_index_miller_indices()
          i_obs = i_obs.customized_copy(indices=indices, info=i_obs.info())
      return i_obs, batches

    i_obs, batches = miller_array_from_mtz(scaled_unmerged)
    if params.reference is not None:
      reference, _ = miller_array_from_mtz(params.reference)
    else:
      reference = None

    return cls(i_obs, batches, params, reference=reference)

  def resolution_auto(self):
    '''Compute resolution limits based on the current self._params set.'''

    if self._params.rmerge:
      stamp("ra: rmerge")
      print('Resolution rmerge:       %.2f' % self.resolution_rmerge())

    if self._params.completeness:
      stamp("ra: comp")
      print('Resolution completeness: %.2f' % self.resolution_completeness())

    if self._params.cc_half:
      stamp("ra: cc")
      print('Resolution cc_half     : %.2f' % self.resolution_cc_half())

    if self._params.cc_ref and self._reference is not None:
      stamp("ra: cc")
      print('Resolution cc_ref      : %.2f' % self.resolution_cc_ref())

    if self._params.isigma:
      stamp("ra: isig")
      print('Resolution I/sig:        %.2f' % self.resolution_unmerged_isigma())

    if self._params.misigma:
      stamp("ra: mnisig")
      print('Resolution Mn(I/sig):    %.2f' % self.resolution_merged_isigma())

    if self._params.i_mean_over_sigma_mean:
      print('Resolution Mn(I)/Mn(sig):    %.2f' %
            self.resolution_i_mean_over_sigma_mean())

  def resolution_rmerge(self, limit=None, log=None):
    '''Compute a resolution limit where either rmerge = 1.0 (limit if
    set) or the full extent of the data. N.B. this fit is only meaningful
    for positive values.'''

    if limit is None:
      limit = self._params.rmerge

    rmerge_s = flex.double(
      [b.r_merge for b in self._merging_statistics.bins]).reversed()
    s_s = flex.double(
      [1/b.d_min**2 for b in self._merging_statistics.bins]).reversed()

    sel = rmerge_s > 0
    rmerge_s = rmerge_s.select(sel)
    s_s = s_s.select(sel)

    if limit == 0.0:
      r_rmerge = 1.0 / math.sqrt(flex.max(s_s))
      rmerge_f = None

    elif limit > flex.max(rmerge_s):
      r_rmerge = 1.0 / math.sqrt(flex.max(s_s))
      rmerge_f = None

    else:
      rmerge_f = log_inv_fit(s_s, rmerge_s, 6)

      if log:
        fout = open(log, 'w')
        for j, s in enumerate(s_s):
          d = 1.0 / math.sqrt(s)
          o = rmerge_s[j]
          m = rmerge_f[j]
          fout.write('%f %f %f %f\n' % (s, d, o, m))
        fout.close()

      try:
        r_rmerge = 1.0 / math.sqrt(interpolate_value(s_s, rmerge_f, limit))
      except Exception:
        r_rmerge = 1.0 / math.sqrt(flex.max(s_s))

    if self._params.plot:
      plot = resolution_plot(ylabel='Rmerge')
      if rmerge_f is not None:
        plot.plot(s_s, rmerge_f, label='fit')
      plot.plot(s_s, rmerge_s, label='Rmerge')
      plot.plot_resolution_limit(r_rmerge)
      plot.savefig('rmerge.png')

    return r_rmerge

  def resolution_i_mean_over_sigma_mean(self, limit=None, log=None):
    '''Compute a resolution limit where either <I>/<sigma> = 1.0 (limit if
    set) or the full extent of the data.'''

    if limit is None:
      limit = self._params.i_mean_over_sigma_mean

    isigma_s = flex.double(
      [b.i_mean_over_sigi_mean for b in self._merging_statistics.bins]).reversed()
    s_s = flex.double(
      [1/b.d_min**2 for b in self._merging_statistics.bins]).reversed()

    sel = isigma_s > 0
    isigma_s = isigma_s.select(sel)
    s_s = s_s.select(sel)

    if flex.min(isigma_s) > limit:
      r_isigma = 1.0 / math.sqrt(flex.max(s_s))
      isigma_f = None

    else:
      isigma_f = log_fit(s_s, isigma_s, 6)

      if log:
        fout = open(log, 'w')
        for j, s in enumerate(s_s):
          d = 1.0 / math.sqrt(s)
          o = isigma_s[j]
          m = isigma_f[j]
          fout.write('%f %f %f %f\n' % (s, d, o, m))
        fout.close()

      try:
        r_isigma = 1.0 / math.sqrt(interpolate_value(s_s, isigma_f, limit))
      except Exception:
        r_isigma = 1.0 / math.sqrt(flex.max(s_s))

    if self._params.plot:
      plot = resolution_plot(ylabel='Unmerged <I>/<sigma>')
      if isigma_f is not None:
        plot.plot(s_s, isigma_f, label='fit')
      plot.plot(s_s, isigma_s, label='Unmerged <I>/<sigma>')
      plot.plot_resolution_limit(r_isigma)
      plot.savefig('i_mean_over_sigma_mean.png')

    return r_isigma

  def resolution_unmerged_isigma(self, limit=None, log=None):
    '''Compute a resolution limit where either I/sigma = 1.0 (limit if
    set) or the full extent of the data.'''

    if limit is None:
      limit = self._params.isigma

    isigma_s = flex.double(
      [b.unmerged_i_over_sigma_mean for b in self._merging_statistics.bins]).reversed()
    s_s = flex.double(
      [1/b.d_min**2 for b in self._merging_statistics.bins]).reversed()

    sel = isigma_s > 0
    isigma_s = isigma_s.select(sel)
    s_s = s_s.select(sel)

    if flex.min(isigma_s) > limit:
      r_isigma = 1.0 / math.sqrt(flex.max(s_s))
      isigma_f = None

    else:
      isigma_f = log_fit(s_s, isigma_s, 6)

      if log:
        fout = open(log, 'w')
        for j, s in enumerate(s_s):
          d = 1.0 / math.sqrt(s)
          o = isigma_s[j]
          m = isigma_f[j]
          fout.write('%f %f %f %f\n' % (s, d, o, m))
        fout.close()

      try:
        r_isigma = 1.0 / math.sqrt(interpolate_value(s_s, isigma_f, limit))
      except Exception:
        r_isigma = 1.0 / math.sqrt(flex.max(s_s))

    if self._params.plot:
      plot = resolution_plot(ylabel='Unmerged I/sigma')
      if isigma_f is not None:
        plot.plot(s_s, isigma_f, label='fit')
      plot.plot(s_s, isigma_s, label='Unmerged I/sigma')
      plot.plot_resolution_limit(r_isigma)
      plot.savefig('isigma.png')

    return r_isigma

  def resolution_merged_isigma(self, limit=None, log=None):
    '''Compute a resolution limit where either Mn(I/sigma) = 1.0 (limit if
    set) or the full extent of the data.'''

    if limit is None:
      limit = self._params.misigma

    misigma_s = flex.double(
      [b.i_over_sigma_mean for b in self._merging_statistics.bins]).reversed()
    s_s = flex.double(
      [1/b.d_min**2 for b in self._merging_statistics.bins]).reversed()

    sel = misigma_s > 0
    misigma_s = misigma_s.select(sel)
    s_s = s_s.select(sel)

    if flex.min(misigma_s) > limit:
      r_misigma = 1.0 / math.sqrt(flex.max(s_s))
      misigma_f = None

    else:
      misigma_f = log_fit(s_s, misigma_s, 6)

      if log:
        fout = open(log, 'w')
        for j, s in enumerate(s_s):
          d = 1.0 / math.sqrt(s)
          o = misigma_s[j]
          m = misigma_f[j]
          fout.write('%f %f %f %f\n' % (s, d, o, m))
        fout.close()

      try:
        r_misigma = 1.0 / math.sqrt(
            interpolate_value(s_s, misigma_f, limit))
      except Exception:
        r_misigma = 1.0 / math.sqrt(flex.max(s_s))

    if self._params.plot:
      plot = resolution_plot(ylabel='Merged I/sigma')
      if misigma_f is not None:
        plot.plot(s_s, misigma_f, label='fit')
      plot.plot(s_s, misigma_s, label='Merged I/sigma')
      plot.plot_resolution_limit(r_misigma)
      plot.savefig('misigma.png')

    return r_misigma

  def resolution_completeness(self, limit=None, log=None):
    '''Compute a resolution limit where completeness < 0.5 (limit if
    set) or the full extent of the data. N.B. this completeness is
    with respect to the *maximum* completeness in a shell, to reflect
    triclinic cases.'''

    if limit is None:
      limit = self._params.completeness

    comp_s = flex.double(
      [b.completeness for b in self._merging_statistics.bins]).reversed()
    s_s = flex.double(
      [1/b.d_min**2 for b in self._merging_statistics.bins]).reversed()

    if flex.min(comp_s) > limit:
      r_comp = 1.0 / math.sqrt(flex.max(s_s))
      comp_f = None

    else:
      comp_f = fit(s_s, comp_s, 6)

      rlimit = limit * max(comp_s)

      if log:
        fout = open(log, 'w')
        for j, s in enumerate(s_s):
          d = 1.0 / math.sqrt(s)
          o = comp_s[j]
          m = comp_f[j]
          fout.write('%f %f %f %f\n' % (s, d, o, m))
        fout.close()

      try:
        r_comp = 1.0 / math.sqrt(
            interpolate_value(s_s, comp_f, rlimit))
      except Exception:
        r_comp = 1.0 / math.sqrt(flex.max(s_s))

    if self._params.plot:
      plot = resolution_plot(ylabel='Completeness')
      if comp_f is not None:
        plot.plot(s_s, comp_f, label='fit')
      plot.plot(s_s, comp_s, label='Completeness')
      plot.plot_resolution_limit(r_comp)
      plot.savefig('completeness.png')

    return r_comp

  def resolution_cc_half(self, limit=None, log=None):
    '''Compute a resolution limit where cc_half < 0.5 (limit if
    set) or the full extent of the data.'''

    if limit is None:
      limit = self._params.cc_half

    if self._params.cc_half_method == 'sigma_tau':
      cc_s = flex.double(
        [b.cc_one_half_sigma_tau for b in self._merging_statistics.bins]).reversed()
    else:
      cc_s = flex.double(
        [b.cc_one_half for b in self._merging_statistics.bins]).reversed()
    s_s = flex.double(
      [1/b.d_min**2 for b in self._merging_statistics.bins]).reversed()

    p = self._params.cc_half_significance_level
    if p is not None:
      if self._params.cc_half_method == 'sigma_tau':
        significance = flex.bool(
          [b.cc_one_half_sigma_tau_significance for b in self._merging_statistics.bins]).reversed()
        cc_half_critical_value = flex.double(
          [b.cc_one_half_sigma_tau_critical_value for b in self._merging_statistics.bins]).reversed()
      else:
        significance = flex.bool(
          [b.cc_one_half_significance for b in self._merging_statistics.bins]).reversed()
        cc_half_critical_value = flex.double(
          [b.cc_one_half_critical_value for b in self._merging_statistics.bins]).reversed()
      # index of last insignificant bin
      i = flex.last_index(significance, False)
      if i is None or i == len(significance) - 1:
        i = 0
      else:
        i += 1
    else:
      i = 0
    if self._params.cc_half_fit == 'tanh':
      cc_f = tanh_fit(s_s[i:], cc_s[i:], iqr_multiplier=4)
    else:
      cc_f = fit(s_s[i:], cc_s[i:], 6)

    stamp("rch: fits")
    rlimit = limit * max(cc_s)

    if log:
      fout = open(log, 'w')
      for j, s in enumerate(s_s):
        d = 1.0 / math.sqrt(s)
        o = cc_s[j]
        m = cc_f[j]
        fout.write('%f %f %f %f\n' % (s, d, o, m))
      fout.close()

    try:
      r_cc = 1.0 / math.sqrt(
          interpolate_value(s_s[i:], cc_f, rlimit))
    except Exception:
      r_cc = 1.0 / math.sqrt(max(s_s[i:]))
    stamp("rch: done : %s" % r_cc)

    if self._params.plot:
      plot = resolution_plot('CC1/2')
      plot.plot(s_s[i:], cc_f, label='fit')
      plot.plot(s_s, cc_s, label='CC1/2')
      if p is not None:
        plot.plot(
          s_s, cc_half_critical_value, label='Confidence limit (p=%g)' %p)
      plot.plot_resolution_limit(r_cc)
      plot.savefig('cc_half.png')

    return r_cc

  def resolution_cc_ref(self, limit=None, log=None):
    '''Compute a resolution limit where cc_ref < 0.5 (limit if
    set) or the full extent of the data.'''

    if limit is None:
      limit = self._params.cc_ref

    intensities = self._intensities.merge_equivalents(
      use_internal_variance=False).array()
    cc_s = flex.double()
    for b in self._merging_statistics.bins:
      sel = intensities.resolution_filter_selection(
        d_min=b.d_min, d_max=b.d_max)
      sel_ref = self._reference.resolution_filter_selection(
        d_min=b.d_min, d_max=b.d_max)
      d = intensities.select(sel)
      dref = self._reference.select(sel_ref)
      cc = d.correlation(dref, assert_is_similar_symmetry=False)
      cc_s.append(cc.coefficient())
    cc_s = cc_s.reversed()

    s_s = flex.double(
      [1/b.d_min**2 for b in self._merging_statistics.bins]).reversed()

    if self._params.cc_half_fit == 'tanh':
      cc_f = tanh_fit(s_s, cc_s, iqr_multiplier=4)
    else:
      cc_f = fit(s_s, cc_s, 6)

    stamp("rch: fits")
    rlimit = limit * max(cc_s)

    if log:
      fout = open(log, 'w')
      for j, s in enumerate(s_s):
        d = 1.0 / math.sqrt(s)
        o = cc_s[j]
        m = cc_f[j]
        fout.write('%f %f %f %f\n' % (s, d, o, m))
      fout.close()

    try:
      r_cc = 1.0 / math.sqrt(
          interpolate_value(s_s, cc_f, rlimit))
    except Exception:
      r_cc = 1.0 / math.sqrt(max(s_s))
    stamp("rch: done : %s" % r_cc)

    if self._params.plot:
      plot = resolution_plot('CCref')
      plot.plot(s_s, cc_f, label='fit')
      plot.plot(s_s, cc_s, label='CCref')
      plot.plot_resolution_limit(r_cc)
      plot.savefig('cc_ref.png')

    return r_cc

def run(args):
  working_phil = phil_defaults
  interp = working_phil.command_line_argument_interpreter(
    home_scope='resolutionizer')
  params, unhandled = interp.process_and_fetch(
    args, custom_processor='collect_remaining')
  params = params.extract().resolutionizer
  if len(unhandled) == 0:
    working_phil.show()
    exit()

  assert len(unhandled) == 1
  scaled_unmerged = unhandled[0]

  stamp("Resolutionizer.py starting")
  m = resolutionizer.from_unmerged_mtz(scaled_unmerged, params)
  stamp("instantiated")
  m.resolution_auto()
  stamp("the end.")


if __name__ == '__main__':
  run(sys.argv[1:])
