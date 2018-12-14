from __future__ import division, absolute_import, print_function

import logging
logger = logging.getLogger(__name__)

import math

import scipy.stats

import libtbx
from libtbx import table_utils
from scitbx.array_family import flex
from cctbx import sgtbx

from dials.algorithms.symmetry import symmetry_base

class determine_space_group(symmetry_base):

  def __init__(self, intensities,
               normalisation='ml_aniso',
               lattice_symmetry_max_delta=2.0,
               d_min=libtbx.Auto,
               min_i_mean_over_sigma_mean=4,
               min_cc_half=0.6,
               relative_length_tolerance=None,
               absolute_angle_tolerance=None):

    super(determine_space_group, self).__init__(
      intensities,
      normalisation=normalisation,
      lattice_symmetry_max_delta=lattice_symmetry_max_delta,
      d_min=d_min,
      min_i_mean_over_sigma_mean=min_i_mean_over_sigma_mean,
      min_cc_half=min_cc_half,
      relative_length_tolerance=relative_length_tolerance,
      absolute_angle_tolerance=absolute_angle_tolerance)

    self.estimate_cc_sig_fac()
    self.estimate_cc_true()
    self.score_symmetry_elements()
    self.score_laue_groups()

  def estimate_cc_sig_fac(self):

    # A1.1. Estimation of sigma(CC) as a function of sample size.

    max_bins = 500
    reflections_per_bin = max(
      200, int(math.ceil(self.intensities.size()/max_bins)))
    binner = self.intensities.setup_binner_counting_sorted(
      reflections_per_bin=reflections_per_bin)

    a = flex.double()
    b = flex.double()
    for i in range(binner.n_bins_all()):
      count = binner.counts()[i]
      if count == 0:
        continue
      bin_isel = binner.array_indices(i)
      p = flex.random_permutation(count)
      p = p[:2 * (count // 2)] # ensure even count
      a.extend(self.intensities.data().select(bin_isel.select(p[:count//2])))
      b.extend(self.intensities.data().select(bin_isel.select(p[count//2:])))

    perm = flex.random_selection(a.size(), min(20000, a.size()))
    a = a.select(perm)
    b = b.select(perm)

    self.corr_unrelated = CorrelationCoefficientAccumulator(a, b)

    n_pairs = a.size()
    min_num_groups = 10 # minimum number of groups
    max_n_group = int(min(n_pairs/min_num_groups, 200)) # maximum number in group
    min_n_group = int(min(5, max_n_group)) # minimum number in group

    mean_ccs = flex.double()
    rms_ccs = flex.double()
    ns = flex.double()
    for n in range(min_n_group, max_n_group):
      ns.append(n)
      ccs = flex.double()
      for i in range(200):
        isel = flex.random_selection(a.size(), n)
        corr = CorrelationCoefficientAccumulator(a.select(isel), b.select(isel))
        ccs.append(corr.coefficient())

      mean_ccs.append(flex.mean(ccs))
      rms_ccs.append(flex.mean(flex.pow2(ccs))**0.5)

    x = 1/flex.pow(ns, 0.5)
    y = rms_ccs
    fit = flex.linear_regression(x, y)

    assert fit.is_well_defined()
    self.cc_sig_fac = fit.slope()

    if 0:
      from matplotlib import pyplot as plt
      plt.plot(x, y)
      plt.plot(
        plt.xlim(), [fit.slope() * x_ + fit.y_intercept() for x_ in plt.xlim()])
      plt.show()

  def estimate_cc_true(self):

    # A1.2. Estimation of E(CC; S).

    # (i)

    var_intensities = flex.mean_and_variance(
      self.intensities.data()).unweighted_sample_variance()
    var_sigmas = flex.mean_and_variance(
      flex.pow2(self.intensities.sigmas())).mean()
    self.E_cc_true = var_intensities/(var_intensities + var_sigmas)

    # (ii)

    reindexed_intensities = self.intensities.change_basis(
      sgtbx.change_of_basis_op('-x,-y,-z')).map_to_asu()
    x, y = self.intensities.common_sets(
      reindexed_intensities, assert_is_similar_symmetry=False)
    self.cc_identity = CorrelationCoefficientAccumulator(x.data(), y.data())

    min_sd = 0.05
    min_sample = 10
    sigma_1 = max(min_sd, self.cc_sig_fac/200**0.5)
    w1 = 0
    w2 = 0
    if sigma_1 > 0.0001:
      w1 = 1/sigma_1**2
    if self.cc_identity.n() > min_sample:
      sigma_2 = max(min_sd, self.cc_sig_fac/self.cc_identity.n()**0.5)
      w2 = 1/sigma_2**2

    assert (w1 + w2) > 0
    self.cc_true = (w1 * self.E_cc_true + w2 * self.cc_identity.coefficient())/(w1 + w2)

    logger.debug('cc_true = w1 * E_cc_true + w2 * cc_identity)/(w1 + w2)')
    logger.debug('w1: %g', w1)
    logger.debug('w2: %g', w2)
    logger.debug('E_cc_true: %g', self.E_cc_true)
    logger.debug('cc_identity: %g', self.cc_identity.coefficient())
    logger.debug('cc_true: %g', self.cc_true)

  def score_symmetry_elements(self):
    self.sym_op_scores = []
    for smx in self.lattice_group.smx():
      if smx.r().info().sense() < 0: continue
      self.sym_op_scores.append(
        ScoreSymmetryElement(self.intensities, smx, self.cc_true, self.cc_sig_fac))

  def score_laue_groups(self):
    subgroup_scores = [
      ScoreSubGroup(subgrp, self.sym_op_scores) for subgrp in self.subgroups.result_groups]
    total_likelihood = sum(score.likelihood for score in subgroup_scores)
    sort_order = flex.sort_permutation(
      flex.double(score.likelihood for score in subgroup_scores), reverse=True)
    self.subgroup_scores = [subgroup_scores[i] for i in sort_order]
    for score in self.subgroup_scores:
      score.likelihood /= total_likelihood

    # The 'confidence' scores are derived from the total probability of the best
    # solution p_best and that for the next best solution p_next:
    #   confidence = [p_best * (p_best - p_next)]^1/2.

    confidence = flex.double(len(self.subgroup_scores), 0)
    for i, score in enumerate(self.subgroup_scores[:-1]):
      next_score = subgroup_scores[i+1]
      if score.likelihood > 0 and next_score.likelihood > 0:
        lgc = score.likelihood * (score.likelihood - next_score.likelihood)
        confidence = abs(lgc)**0.5
        if lgc < 0: confidence = -confidence
        score.confidence = confidence

    self.best_solution = self.subgroup_scores[0]

  def __str__(self):
    output = []
    output.append('Input crystal symmetry:')
    output.append(str(self.input_intensities[0].space_group_info()))
    output.append(str(self.median_unit_cell))
    output.append('Change of basis op to minimum cell: %s' % self.cb_op_inp_min)
    output.append('Crystal symmetry in minimum cell:')
    output.append(str(self.intensities.space_group_info()))
    output.append(str(self.intensities.unit_cell()))
    output.append('Lattice point group: %s' % self.lattice_group.info())
    output.append('Overall CC for %i unrelated pairs: %.3f' % (
      self.corr_unrelated.n(), self.corr_unrelated.coefficient()))
    output.append(
      'Estimated expectation value of true correlation coefficient E(CC) = %.3f'
      % self.E_cc_true)
    output.append('Estimated sd(CC) = %.3f / sqrt(N)' % self.cc_sig_fac)
    output.append(
      'Estimated E(CC) of true correlation coefficient from identity = %.3f'
      % self.cc_true)

    header = ('likelihood', 'Z-CC', 'CC', 'N', '', 'Operator')
    rows = [header]
    for score in self.sym_op_scores:
      if score.likelihood > 0.9: stars = '***'
      elif score.likelihood > 0.7: stars = '**'
      elif score.likelihood > 0.5: stars = '*'
      else: stars = ''
      rows.append((
        '%.3f' % score.likelihood,
        '%.2f' % score.z_cc,
        '%.2f' % score.cc.coefficient(),
        '%i' % score.n_refs,
        stars,
        '%s' % score.sym_op.r().info()))
    output.append('Scoring individual symmetry elements')
    output.append(table_utils.format(rows, has_header=True, delim='  '))

    header = ('Patterson group', '', 'Likelihood', 'NetZcc', 'Zcc+', 'Zcc-',
              'CC', 'CC-', 'delta', 'Reindex operator')
    rows = [header]
    for score in self.subgroup_scores:
      if score.likelihood > 0.8: stars = '***'
      elif score.likelihood > 0.6: stars = '**'
      elif score.likelihood > 0.4: stars = '*'
      else: stars = ''
      rows.append((
        '%s' % score.subgroup['best_subsym'].space_group_info(),
        stars,
        '%.3f' % score.likelihood,
        '% .2f' % score.z_cc_net,
        '% .2f' % score.z_cc_for,
        '% .2f' % score.z_cc_against,
        '% .2f' % score.cc_for.coefficient(),
        '% .2f' % score.cc_against.coefficient(),
        '%.1f' % score.subgroup['max_angular_difference'],
        '%s' % (score.subgroup['cb_op_inp_best'] * self.cb_op_inp_min)))
    output.append('Scoring all possible sub-groups')
    output.append(table_utils.format(rows, has_header=True, delim='  '))

    output.append(
      'Best solution: %s' % self.best_solution.subgroup['best_subsym'].space_group_info())
    output.append(
      'Unit cell: %s' % self.best_solution.subgroup['best_subsym'].unit_cell())
    output.append(
      'Reindex operator: %s' % (
        self.best_solution.subgroup['cb_op_inp_best'] * self.cb_op_inp_min))
    output.append('Laue group probability: %.3f' % self.best_solution.likelihood)
    output.append('Laue group confidence: %.3f' % self.best_solution.confidence)
    return "\n".join(output)

  def as_dict(self):
    d = {
      'input_symmetry': {
        'hall_symbol': self.input_intensities[0].space_group().type().hall_symbol(),
        'unit_cell': self.median_unit_cell.parameters(),
      },
      'cb_op_inp_min': self.cb_op_inp_min.as_xyz(),
      'min_cell_symmetry': {
        'hall_symbol': self.intensities.space_group().type().hall_symbol(),
        'unit_cell': self.intensities.unit_cell().parameters(),
      },
      'lattice_point_group': self.lattice_group.type().hall_symbol(),
      'cc_unrelated_pairs': self.corr_unrelated.coefficient(),
      'n_unrelated_pairs': self.corr_unrelated.n(),
      'E_cc_true': self.E_cc_true,
      'cc_sig_fac': self.cc_sig_fac,
      'cc_true': self.cc_true,
    }

    d['sym_op_scores'] = [score.as_dict() for score in self.sym_op_scores]
    d['subgroup_scores'] = [score.as_dict() for score in self.subgroup_scores]
    return d

  def as_json(self, filename=None, indent=2):
    d = self.as_dict()
    import json
    if filename is not None:
      return json.dump(d, open(filename, 'wb'), indent=indent)
    else:
      return json.dumps(d, indent=indent)

class ScoreSymmetryElement(object):

  def __init__(self, intensities, sym_op, cc_true, cc_sig_fac):

    self.sym_op = sym_op
    assert self.sym_op.r().info().sense() >= 0
    self.cc = CorrelationCoefficientAccumulator()
    cb_op = sgtbx.change_of_basis_op(self.sym_op)
    cb_ops = [cb_op]
    if self.sym_op.r().order() > 2:
      # include inverse symmetry operation
      cb_ops.append(cb_op.inverse())
    for cb_op in cb_ops:
      if cb_op.is_identity_op():
        cb_op = sgtbx.change_of_basis_op('-x,-y,-z')
      reindexed_intensities = intensities.change_basis(cb_op).map_to_asu()
      x, y = intensities.common_sets(
        reindexed_intensities, assert_is_similar_symmetry=False)
      tmp_sg = sgtbx.space_group().expand_smx(self.sym_op) \
        .info().as_reference_setting().group()
      sel = ~(tmp_sg.is_sys_absent(x.indices()))
      sel &= tmp_sg.epsilon(x.indices()) == 1
      x = x.select(sel)
      y = y.select(sel)
      self.cc += CorrelationCoefficientAccumulator(x.data(), y.data())

    self.n_refs = self.cc.n()
    if self.n_refs <= 0:
      self.likelihood = 0
      self.z_cc = 0
      return

    self.sigma_cc = max(0.1, cc_sig_fac/self.n_refs**0.5)
    self.z_cc = self.cc.coefficient()/self.sigma_cc

    # Probability of observing this CC if the sym op is present
    # Modelled by a Cauchy distribution centred on cc_true and width gamma = sigma_cc
    # p(CC; S)
    self.p_cc_given_s = trunccauchy_pdf(self.cc.coefficient(), -1, 1, cc_true, self.sigma_cc)

    def DM_power_pdf(m, k=2):
      return (1.0-pow(m,k))**(1.0/k)

    def numerator(x, loc, scale, k):
      return trunccauchy_pdf(x, -1, 1, loc=loc, scale=scale) * DM_power_pdf(x, k=k)

    def denominator(x, k):
      return DM_power_pdf(x, k=k)

    # Probability of observing this CC if the sym op is NOT present
    # p(CC; !S)

    k = 2
    sump = scipy.integrate.quad(
      numerator, 0, 1, args=(self.cc.coefficient(), self.sigma_cc, k))[0]
    sumw = scipy.integrate.quad(denominator, 0, 1, args=(k,))[0]
    self.p_cc_given_not_s = sump/sumw

    # Likelihood of symmetry element being present
    # p(S; CC) = p(CC; S) / (p(CC; S) + p(CC; !S))
    self.likelihood = self.p_cc_given_s/(self.p_cc_given_s + self.p_cc_given_not_s)

  def __str__(self):
    return '%.3f %.2f %.2f %i %s' % (self.likelihood,
                                     self.z_cc,
                                     self.cc.coefficient(),
                                     self.n_refs,
                                     self.sym_op.r().info())

  def as_dict(self):
    return {
      'likelihood': self.likelihood,
      'z_cc': self.z_cc,
      'cc': self.cc.coefficient(),
      'n_ref': self.n_refs,
      'operator': self.sym_op.as_xyz(),
    }

class ScoreSubGroup(object):
  def __init__(self, subgroup, sym_op_scores):

    # Combined correlation coefficients for symmetry operations
    # present/absent from subgroup
    self.subgroup = subgroup
    cb_op_inp_best = subgroup['cb_op_inp_best']
    patterson_group = subgroup['best_subsym'].space_group().change_basis(
      cb_op_inp_best.inverse())
    self.cc_for = CorrelationCoefficientAccumulator()
    self.cc_against = CorrelationCoefficientAccumulator()
    for score in sym_op_scores:
      if score.sym_op in patterson_group:
        self.cc_for += score.cc
      else:
        self.cc_against += score.cc

    # Overall Zcc scores for symmetry elements present/absent from subgroup
    self.z_cc_for = 0
    self.z_cc_against = 0
    n_for = 0
    n_against = 0
    PL_for = 0
    PL_against = 0
    power = 2
    for score in sym_op_scores:
      if score.n_refs <= 2:
        continue
      if score.sym_op in patterson_group:
        self.z_cc_for += score.z_cc**power
        n_for += 1
        PL_for += math.log(score.p_cc_given_s)
      else:
        self.z_cc_against += score.z_cc**power
        n_against += 1
        PL_against += math.log(score.p_cc_given_not_s)

    # Overall likelihood for this subgroup
    self.likelihood = math.exp(PL_for + PL_against)

    if n_against > 0:
      self.z_cc_against = (self.z_cc_against/n_against)**(1/power)
    if n_for > 0:
      self.z_cc_for = (self.z_cc_for/n_for)**(1/power)
    self.z_cc_net = self.z_cc_for - self.z_cc_against
    self.confidence = 0

  def __str__(self):
    return '%s %.3f %.2f %.2f %.2f %.2f %.2f' % (
      self.subgroup['best_subsym'].space_group_info(),
      self.likelihood,
      self.z_cc_net,
      self.z_cc_for,
      self.z_cc_against,
      self.cc_for.coefficient(),
      self.cc_against.coefficient())

  def as_dict(self):
    return {
      'patterson_group': self.subgroup['best_subsym'].space_group().type().hall_symbol(),
      'likelihood': self.likelihood,
      'confidence': self.confidence,
      'z_cc_net': '% .2f' % self.z_cc_net,
      'z_cc_for': '% .2f' % self.z_cc_for,
      'z_cc_against': '% .2f' % self.z_cc_against,
      'cc_for': '% .2f' % self.cc_for.coefficient(),
      'cc_against': '% .2f' % self.cc_against.coefficient(),
      'max_angular_difference': '%.1f' % self.subgroup['max_angular_difference'],
      'cb_op': '%s' % (self.subgroup['cb_op_inp_best'])
    }


# Single-pass formula for Pearson correlation coefficient
# https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample
class CorrelationCoefficientAccumulator(object):
  def __init__(self, x=None, y=None):
    self._n = 0
    self._sum_x = 0
    self._sum_y = 0
    self._sum_xy = 0
    self._sum_x_sq = 0
    self._sum_y_sq = 0
    if x is not None and y is not None:
      self.accumulate(x, y)

  def accumulate(self, x, y):
    assert x.size() == y.size()
    self._n += x.size()
    self._sum_x += flex.sum(x)
    self._sum_y += flex.sum(y)
    self._sum_xy += flex.sum(x * y)
    self._sum_x_sq += flex.sum(flex.pow2(x))
    self._sum_y_sq += flex.sum(flex.pow2(y))

  def coefficient(self):
    if self._n == 0: return 0
    return self.numerator()/self.denominator()

  def n(self):
    return self._n

  def numerator(self):
    return self._n * self._sum_xy - self._sum_x * self._sum_y

  def denominator(self):
    return (math.sqrt(self._n * self._sum_x_sq - self._sum_x**2) *
            math.sqrt(self._n * self._sum_y_sq - self._sum_y**2))

  def __iadd__(self, other):
    self._n += other._n
    self._sum_x += other._sum_x
    self._sum_y += other._sum_y
    self._sum_xy += other._sum_xy
    self._sum_x_sq += other._sum_x_sq
    self._sum_y_sq += other._sum_y_sq
    return self


def trunccauchy_pdf(x, a, b, loc=0, scale=1):
  assert b > a
  rv = scipy.stats.cauchy(loc=loc, scale=scale)
  return rv.pdf(x) / (rv.cdf(b) - rv.cdf(a))

